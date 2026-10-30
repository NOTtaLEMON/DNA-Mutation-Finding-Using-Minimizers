#ifndef BLOOM_H
#define BLOOM_H

#include <stddef.h>
#include <stdint.h>

#include "seq.h"
#include "minimizer.h"

/*
 * A Bloom filter is a compact bit array used only to answer membership
 * questions. It does NOT store reference positions; Stage 3's Index does that.
 *
 * `bits`       Packed bit array. A bit set to 1 means one or more inserted
 *              keys mapped to that bit.
 * `bit_count`  Number of usable bits (not bytes).
 * `byte_count` Number of allocated bytes in `bits`.
 * `hash_count` Number of bit positions checked for each key.
 * `expected_items` Number of minimizer occurrences used to size the filter.
 *                 This can be larger than the number of distinct words,
 *                 because repeated minimizer occurrences may be inserted.
 * `inserted_items` Number of occurrences inserted (duplicates are harmless).
 */
typedef struct {
    uint8_t *bits;
    size_t bit_count;
    size_t byte_count;
    size_t hash_count;
    size_t expected_items;
    size_t inserted_items;
} BloomFilter;

/*
 * Build a Bloom filter from the reference minimizers.
 *
 * Each key is the pair (minimizer hash, packed k-mer word ID). This matches
 * Stage 3's exact identity check: Stage 3 checks BOTH the hash and word_id.
 * Using the pair avoids treating a 64-bit minimizer-hash collision as proof
 * that two different DNA words are equal. Bloom-filter collisions are still
 * possible, but those only create false positives, which Stage 3 verifies.
 *
 * `seq`          Reference sequence loaded by Stage 1.
 * `reference_mins` Minimizers generated from that same reference by Stage 2.
 * `k`            K-mer length used by Stage 2; must be 1..32 because Stage 3
 *                stores a packed word ID in uint64_t.
 * `false_positive_rate` Desired sizing target, strictly between 0 and 1
 *                (for example, 0.01 means 1%). It is a theoretical target,
 *                not an absolute guarantee for every finite query set.
 * `out`          Must be zero-initialized or previously passed to bloom_free.
 *
 * Returns 0 on success, -1 for invalid input, size overflow, or allocation
 * failure. On failure, `out` is left empty.
 *
 * DESIGN NOTE: The filter is built from the reference minimizer set, as shown
 * in the project's Step 4 reference branch. Repeated occurrences may insert
 * the same key more than once. Setting an already-set Bloom bit is harmless;
 * sizing by occurrence count may simply allocate somewhat more bits than
 * would be needed for the number of distinct minimizer words.
 *
 * N POLICY: Stage 2 must already have excluded k-mers overlapping N/unknown
 * positions. This function checks sequence bounds but relies on that Stage 2
 * contract for N handling, so it does not index artificial placeholder As.
 */
int bloom_build_from_minimizers(const Seq *seq,
                                const MinimizerSet *reference_mins,
                                size_t k,
                                double false_positive_rate,
                                BloomFilter *out);

/*
 * Ask whether (hash, word_id) might be in the reference filter.
 *
 * Returns:
 *   1  = POSSIBLY PRESENT: all required bits are set. Call index_lookup().
 *   0  = DEFINITELY ABSENT: at least one required bit is zero. Skip lookup.
 *  -1  = invalid/corrupted filter argument.
 *
 * A return value of 1 is NOT proof of membership; Bloom filters can have
 * false positives. A correctly built filter has no false negatives for keys
 * that were inserted and queried using the same pair and hash scheme.
 */
int bloom_maybe_contains(const BloomFilter *filter,
                         uint64_t hash,
                         uint64_t word_id);

/* Frees the bit array and resets the filter to empty. Safe on an empty filter. */
void bloom_free(BloomFilter *filter);

#endif /* BLOOM_H */
