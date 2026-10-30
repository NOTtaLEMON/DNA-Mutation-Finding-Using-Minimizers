# Stage 4: Bloom filter (C)

## 1. Where the files go

Create this folder in the repository:

```text
src/
  stage1_load/
    seq.c
    seq.h
  stage2_minimizers/
    minimizer.c
    minimizer.h
  stage3_index/
    index.c
    index.h
  stage4_bloom/              <-- new folder
    bloom.c                  <-- new
    bloom.h                  <-- new
    test_bloom.c             <-- new tests
```

Do not overwrite Stages 1, 2 or 3.

## 2. What the Bloom filter does

The Bloom filter is a compact bit array that quickly answers a membership question:

- **Definitely absent (0):** at least one required bit is zero. The queried key was not inserted, so the Stage 3 hash-table lookup can safely be skipped.
- **Possibly present (1):** all required bits are one. The key might be present, so call Stage 3 `index_lookup()` to confirm and retrieve reference positions.
- **Error (-1):** the Bloom-filter object is invalid.

A Bloom filter can have false positives (it says “possibly present” for an absent key). It must not have false negatives for keys inserted and queried consistently. The Bloom filter does not store minimizer strings or reference positions. Stage 3's hash table remains the authoritative structure for retrieving positions.

## 3. What is inserted into the filter?

Each key is the pair:

```text
(Stage 2 minimizer hash, Stage 3 packed word_id)
```

This matches Stage 3's identity check, which requires both `hash` and `word_id` to match. It is safer than storing only the Stage 2 hash: two different k-mers could theoretically have the same 64-bit hash. A Bloom false positive is harmless because Stage 3 performs the exact check afterwards.

The Bloom filter is built from the reference minimizers, in parallel conceptually with building the Stage 3 index. Repeated minimizer occurrences may insert the same key multiple times; that only sets already-set bits. The filter is sized using the occurrence count, so duplicates can make it somewhat larger than necessary.

The Bloom filter uses its own two derived 64-bit hash streams to calculate several bit-array positions. These are **not** the same thing as the Stage 2 minimizer hash. Stage 2's hash helps rank/index minimizers; the Bloom hash streams choose bits in the Bloom array.

## 4. Public API

```c
int bloom_build_from_minimizers(const Seq *seq,
                                const MinimizerSet *reference_mins,
                                size_t k,
                                double false_positive_rate,
                                BloomFilter *out);

int bloom_maybe_contains(const BloomFilter *filter,
                         uint64_t hash,
                         uint64_t word_id);

void bloom_free(BloomFilter *filter);
```

Use a target false-positive rate such as `0.01` (1%) as a starting example. It is a sizing target, not a promise that every arbitrary query set will have exactly that observed rate.

`k` must be in the range 1..32 because Stage 3 stores up to 32 two-bit bases in a `uint64_t` word ID.

The output filter must be zero-initialized or already freed before calling `bloom_build_from_minimizers()`:

```c
BloomFilter reference_bloom = {0};
```

## 5. Build the two reference-side structures

After Stage 2 has generated minimizers for the reference, build the Stage 3 index and Stage 4 Bloom filter from the same reference and minimizer set:

```c
const size_t k = 15;
const size_t window_kmers = 10;
const double target_fp_rate = 0.01;

MinimizerSet reference_mins = {0};
Index reference_index = {0};
BloomFilter reference_bloom = {0};

if (minimizers_generate(&reference, k, window_kmers,
                        &reference_mins) != 0) {
    /* Handle Stage 2 failure. */
}

if (index_build(&reference, &reference_mins, k, 0,
                &reference_index) != 0) {
    /* Handle Stage 3 failure. */
}

if (bloom_build_from_minimizers(&reference, &reference_mins, k,
                                target_fp_rate, &reference_bloom) != 0) {
    /* Handle Stage 4 failure. */
}
```

The two structures are conceptually parallel products of the reference minimizers:

```text
                         Reference minimizers
                           /              \
                          /                \
                 Stage 3 Index         Stage 4 Bloom filter
                 stores positions      stores membership bits
```

For an empty reference minimizer set, the builder creates a valid empty filter; every query returns definitely absent.

## 6. Query a sample read

Stage 2 first generates minimizers for one sample read. For each one, compute the packed word ID in exactly the same way Stage 3 does, then ask the Bloom filter:

```c
MinimizerSet read_mins = {0};
size_t i;

if (minimizers_generate(&read, k, window_kmers, &read_mins) != 0) {
    /* Handle Stage 2 failure. */
}

for (i = 0; i < read_mins.count; i++) {
{
    uint64_t hash = read_mins.items[i].hash;
    uint64_t word_id = word_id_at(&read, read_mins.items[i].pos, k);
    int answer = bloom_maybe_contains(&reference_bloom, hash, word_id);

    if (answer == 0) {
        /* DEFINITELY ABSENT: skip this minimizer's hash-table lookup. */
        continue;
    }
    if (answer < 0) {
        /* Invalid filter: handle as an error; do NOT treat as absent. */
        /* Break/return according to the calling function's error policy. */
        break;
    }

    /* POSSIBLY PRESENT: the hash table confirms membership and gives positions. */
    {
        size_t probes = 0;
        const IndexEntry *entry = index_lookup(&reference_index,
                                               hash,
                                               word_id,
                                               &probes);
        if (entry != NULL) {
            const PosNode *p;
            for (p = entry->pos_head; p != NULL; p = p->next) {
                /* Use p->pos as a reference occurrence in Stage 5 voting. */
            }
        }
        /* entry == NULL is possible because Bloom filters allow false positives. */
    }
}

minimizers_free(&read_mins);
```

The sample read's `word_id` is calculated from the read sequence, while the Bloom filter was built using the corresponding reference minimizer's word ID. Equal DNA k-mers, generated with the same k and encoding, have the same word ID.

## 7. How the Bloom filter chooses bits

For each key pair, the implementation derives two 64-bit hash values and uses double hashing:

`bit_i = (h1 + i * h2) % bit_count`

for `i = 0 .. hash_count-1`.

Insertion sets each selected bit to 1. Querying checks those same bits:

- If any is 0, the key is definitely absent.
- If all are 1, the key is possibly present.

The bit array is compact: eight membership bits are stored in each byte.

The filter dimensions use the standard Bloom-filter formulas:

`m = -n * ln(p) / (ln(2)^2)`  (bit count)

`h = (m / n) * ln(2)`         (number of hash functions)

where `n` is the number of reference minimizer occurrences used for sizing and `p` is the target false-positive rate. The implementation rounds the bit count upward and the hash-function count to a positive integer.

## 8. Compile and run tests

From the repository root, after completing the Stage 1 TODOs (`base_to_code()` and `seq_get()`):

```sh
cc -std=c11 -Wall -Wextra -Wpedantic \
  -Isrc/stage1_load -Isrc/stage2_minimizers \
  -Isrc/stage3_index -Isrc/stage4_bloom \
  src/stage1_load/seq.c \
  src/stage2_minimizers/minimizer.c \
  src/stage3_index/index.c \
  src/stage4_bloom/bloom.c \
  src/stage4_bloom/test_bloom.c \
  -lm -o test_bloom

./test_bloom
```

`-lm` links the C math library for `log()`, `ceil()`, `floor()` and `isfinite()` on toolchains where it is a separate library (including common GCC/Linux setups).

Expected success message:

```text
all Bloom filter checks passed
```

The tests cover inserted keys (no false negatives), a definitely-absent query, an empty filter, and invalid inputs. The absent-key test searches for a definitely-absent example because a Bloom filter is allowed to return a false positive for any particular absent key.

## 9. Important assumptions and items to revisit

1. **Stage 2 / Stage 3 compatibility:** Both stages must use the same k-mer length, minimizer hash function, and position convention. The Bloom filter uses the exact `(hash, word_id)` pair expected by Stage 3.
2. **Unknown bases:** Stage 2 must exclude every k-mer overlapping an N/unknown position. This Stage 4 builder trusts that contract.
3. **Hash collisions:** Bloom membership is not proof. Always call `index_lookup()` after a “possibly present” result.
4. **False-positive rate:** `0.01` is an initial example. We should benchmark memory use and the observed false-positive rate on the actual reference dataset before choosing a final value.
5. **The Bloom filter is an optimization:** It must never be used to produce reference positions or mutation calls. It only decides whether a Stage 3 lookup can safely be skipped.
6. **Finimizers:** This code uses the existing fixed-k minimizer output. It does not implement the frequency-bounded finimizer method from the 2025 paper.
7. **Output ownership:** Pass a zero-initialized/empty filter to the build function. Call `bloom_free()` before rebuilding into an already-used filter.
