# Stage 2: Minimizer generation (C)

## Where these files go in the project

Copy this folder into the repository:

```text
src/
  stage1_load/
    seq.c
    seq.h
  stage2_minimizers/
    minimizer.c       <-- new
    minimizer.h       <-- new
    test_minimizer.c  <-- new tests for this stage
```

## Important design choice

This implementation uses the common read-mapping convention:

- `k` = number of bases in each k-mer.
- `window_kmers` = number of consecutive k-mers in one window.
- Hash every valid k-mer and choose the smallest hash in each window.

This matches the project's candidate-start formula:

`reference_kmer_start - read_kmer_start`

It is **not** the finimizer method from the 2025 paper. The paper's formal definition uses an m-mer selected inside a k-mer, and finimizers introduce frequency-dependent variable lengths. Do not mix these definitions. The project specification/source code should decide which convention is ultimately required.

## Stage 1 prerequisite: finish the two TODO functions

In `src/stage1_load/seq.c`, implement the existing TODO functions as follows:

```c
int base_to_code(char c)
{
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

uint8_t seq_get(const Seq *s, size_t pos)
{
    size_t byte_index;
    unsigned int shift;

    /* The public contract says pos < s->len. The guard prevents accidental
       invalid access if a caller violates that contract. */
    if (s == NULL || s->cells == NULL || pos >= s->len)
        return 0;

    byte_index = pos / 4;
    shift = (unsigned int)(6 - 2 * (pos % 4));
    return (uint8_t)(((unsigned int)s->cells[byte_index] >> shift) & 0x03u);
}
```

Add `#include <stddef.h>` only if it is not already indirectly available; `seq.h` already includes it.

## API

```c
int minimizers_generate(const Seq *seq,
                        size_t k,
                        size_t window_kmers,
                        MinimizerSet *out);
void minimizers_free(MinimizerSet *set);
```

`out` must be initialized as an empty set:

```c
MinimizerSet read_minimizers = {0};
```

For each result:

- `items[i].hash` is the 64-bit FNV-1a hash used for ordering.
- `items[i].pos` is the zero-based start position of the selected k-mer in the input sequence.

The same `k` and `window_kmers` settings and same hash/order definition must be used when generating minimizers for the reference and for reads.

Example call (the reference and each read use the same function):

```c
MinimizerSet reference_minimizers = {0};
MinimizerSet read_minimizers = {0};

if (minimizers_generate(&reference, 15, 10, &reference_minimizers) != 0) {
    /* Handle invalid parameters or allocation failure. */
}

if (minimizers_generate(&one_read, 15, 10, &read_minimizers) != 0) {
    /* Handle invalid parameters or allocation failure. */
}

/* Later stages can use items[i].hash and items[i].pos. */

minimizers_free(&reference_minimizers);
minimizers_free(&read_minimizers);
```

The values `15` and `10` are illustrative only, not recommended final parameters. A later benchmark should select values appropriate to the read lengths and dataset.

## Unknown (`N`) bases

The Stage 1 loader turns non-ACGT characters into a packed placeholder A and stores their positions in `Seq.runs`. This module checks `Seq.runs` and treats every k-mer overlapping an unknown position as invalid. It also breaks the stream at unknowns, so a minimizer window never bridges an N.

## Window behaviour

A window contains `window_kmers` consecutive valid k-mers, not `window_kmers` DNA bases. If a valid segment is shorter than the requested window, this implementation emits the minimum from that short segment. This is an explicit, documented policy; confirm it against the project specification before finalizing.

Repeated selection of the same k-mer occurrence in adjacent windows is suppressed. Equal hash values retain the leftmost occurrence.

## Compile and run the tests

From the repository root (adjust paths if you put `test_minimizer.c` elsewhere):

```sh
cc -std=c11 -Wall -Wextra -Wpedantic \
  -Isrc/stage1_load -Isrc/stage2_minimizers \
  src/stage1_load/seq.c \
  src/stage2_minimizers/minimizer.c \
  src/stage2_minimizers/test_minimizer.c \
  -o test_minimizer

./test_minimizer
```

The tests need Stage 1's `base_to_code()` and `seq_get()` to be implemented first.

## Important limitations to revisit before Stage 3

1. **Reverse-complement reads:** Step 2 currently hashes the forward orientation only. If reads may be reverse-complemented relative to the reference, Stage 5 must handle strand/orientation or Step 2 must be extended to canonical minimizers with strand metadata.
2. **Hash collisions:** a 64-bit hash is used to rank k-mers. A later index must not treat a hash match as unquestionable proof of sequence equality; verify the actual k-mer when needed.
3. **Performance:** this first correct implementation hashes each k-mer by reading its `k` bases, so hashing costs O(k) per k-mer. A rolling hash can optimize this later after correctness is established.
4. **Parameter definition:** check whether your project report specifies window size in bases or in number of k-mers. This code uses number of k-mers.
5. **Finimizers:** not implemented here. Finimizers are a separate variable-length, frequency-bounded indexing approach and would require a compatible index design.
