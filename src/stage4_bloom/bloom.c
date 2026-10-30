#include "bloom.h"

#include "index.h"  /* Reuse Stage 3's exact packed-word identity function. */

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Natural logarithm of 2, used in the standard Bloom-filter sizing formula. */
#define BLOOM_LN2 0.69314718055994530942

/*
 * SplitMix64 finalizer/mixer.
 *
 * A Bloom filter needs several well-distributed bit positions for each key.
 * We derive two 64-bit values from the Stage 2 hash and Stage 3 word ID, then
 * use double hashing to generate the remaining positions efficiently.
 */
static uint64_t mix64(uint64_t x)
{
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

/* Rotate left without relying on compiler-specific intrinsics. */
static uint64_t rotl64(uint64_t x, unsigned int r)
{
    return (x << r) | (x >> (64U - r));
}

/*
 * Turn the full key pair (Stage 2 hash, packed DNA word ID) into two hash
 * streams. The second stream is made odd; this is a common double-hashing
 * technique and helps visit different positions in the bit array.
 */
static void make_hash_pair(uint64_t hash,
                           uint64_t word_id,
                           uint64_t *h1,
                           uint64_t *h2)
{
    *h1 = mix64(hash ^ rotl64(word_id, 23) ^ UINT64_C(0x243f6a8885a308d3));
    *h2 = mix64(word_id ^ rotl64(hash, 41) ^ UINT64_C(0x13198a2e03707344));
    *h2 |= UINT64_C(1);
}

/*
 * Calculate Bloom-filter dimensions from the expected number of inserted
 * items and target false-positive rate.
 *
 * Standard formula:
 *   m = -n * ln(p) / (ln(2)^2)   (number of bits)
 *   h = (m/n) * ln(2)            (number of hash functions)
 *
 * We round m upward and h to the nearest positive integer. We explicitly
 * check conversions/allocations to avoid overflowing size_t.
 */
static int bloom_dimensions(size_t n,
                            double p,
                            size_t *bit_count_out,
                            size_t *byte_count_out,
                            size_t *hash_count_out)
{
    double m_real;
    double h_real;
    size_t bit_count;
    size_t byte_count;
    size_t hash_count;

    if (n == 0 || !isfinite(p) || p <= 0.0 || p >= 1.0)
        return -1;

    m_real = -((double)n * log(p)) / (BLOOM_LN2 * BLOOM_LN2);
    if (!isfinite(m_real) || m_real < 1.0 || m_real > (double)SIZE_MAX - 8.0)
        return -1;

    bit_count = (size_t)ceil(m_real);
    if (bit_count == 0)
        bit_count = 1;

    /* Compute ceil(bit_count / 8) without overflowing bit_count + 7. */
    byte_count = bit_count / 8;
    if (bit_count % 8 != 0) {
        if (byte_count == SIZE_MAX)
            return -1;
        byte_count++;
    }

    h_real = ((double)bit_count / (double)n) * BLOOM_LN2;
    if (!isfinite(h_real) || h_real > (double)SIZE_MAX - 1.0)
        return -1;

    hash_count = (size_t)floor(h_real + 0.5);
    if (hash_count == 0)
        hash_count = 1;

    *bit_count_out = bit_count;
    *byte_count_out = byte_count;
    *hash_count_out = hash_count;
    return 0;
}

/* Return the bit-array index for one of the double-hashed positions. */
static size_t bit_index_for(uint64_t h1,
                            uint64_t h2,
                            size_t round,
                            size_t bit_count)
{
    /* Unsigned arithmetic wrapping is defined by C and is intentional here. */
    uint64_t combined = h1 + (uint64_t)round * h2;
    return (size_t)(combined % (uint64_t)bit_count);
}

/* Set all Bloom bits for one key. Called only after successful allocation. */
static void bloom_insert_key(BloomFilter *filter,
                             uint64_t hash,
                             uint64_t word_id)
{
    uint64_t h1;
    uint64_t h2;
    size_t i;

    make_hash_pair(hash, word_id, &h1, &h2);

    for (i = 0; i < filter->hash_count; i++) {
        size_t bit = bit_index_for(h1, h2, i, filter->bit_count);
        size_t byte = bit / 8;
        unsigned int offset = (unsigned int)(bit % 8);

        /* Bit numbering within a byte is least-significant-bit first. */
        filter->bits[byte] |= (uint8_t)(1U << offset);
    }

    filter->inserted_items++;
}

int bloom_build_from_minimizers(const Seq *seq,
                                const MinimizerSet *reference_mins,
                                size_t k,
                                double false_positive_rate,
                                BloomFilter *out)
{
    BloomFilter filter;
    size_t i;

    if (out == NULL)
        return -1;

    /* Contract: caller supplies an empty/zero-initialized filter. */
    memset(out, 0, sizeof *out);

    if (seq == NULL || reference_mins == NULL || k == 0 || k > 32)
        return -1;
    if (reference_mins->count > 0 && reference_mins->items == NULL)
        return -1;
    if (reference_mins->count > 0 && (seq->cells == NULL || k > seq->len))
        return -1;
    if (!isfinite(false_positive_rate) ||
        false_positive_rate <= 0.0 || false_positive_rate >= 1.0)
        return -1;

    memset(&filter, 0, sizeof filter);
    filter.expected_items = reference_mins->count;

    /* An empty reference minimizer set is valid: every query is definitely absent. */
    if (reference_mins->count == 0) {
        *out = filter;
        return 0;
    }

    if (bloom_dimensions(reference_mins->count,
                         false_positive_rate,
                         &filter.bit_count,
                         &filter.byte_count,
                         &filter.hash_count) != 0)
        return -1;

    filter.bits = (uint8_t *)calloc(filter.byte_count, sizeof *filter.bits);
    if (filter.bits == NULL)
        return -1;

    for (i = 0; i < reference_mins->count; i++) {
        size_t pos = reference_mins->items[i].pos;
        uint64_t word_id;

        /* Check the complete k-mer is inside the reference before reading it. */
        if (pos > seq->len - k) {
            free(filter.bits);
            return -1;
        }

        /*
         * Stage 3 uses this same function to define exact word identity.
         * Keeping the same definition prevents Bloom/index disagreement.
         */
        word_id = word_id_at(seq, pos, k);
        bloom_insert_key(&filter, reference_mins->items[i].hash, word_id);
    }

    *out = filter;
    return 0;
}

int bloom_maybe_contains(const BloomFilter *filter,
                         uint64_t hash,
                         uint64_t word_id)
{
    uint64_t h1;
    uint64_t h2;
    size_t i;

    if (filter == NULL)
        return -1;

    /* A deliberately empty filter represents an empty reference minimizer set. */
    if (filter->bit_count == 0 && filter->byte_count == 0 &&
        filter->hash_count == 0 && filter->bits == NULL &&
        filter->expected_items == 0)
        return 0;

    /* Catch partially initialized/corrupted filters instead of reading invalid memory. */
    if (filter->bits == NULL || filter->bit_count == 0 ||
        filter->byte_count == 0 || filter->hash_count == 0 ||
        filter->byte_count < filter->bit_count / 8 + (filter->bit_count % 8 != 0))
        return -1;

    make_hash_pair(hash, word_id, &h1, &h2);

    for (i = 0; i < filter->hash_count; i++) {
        size_t bit = bit_index_for(h1, h2, i, filter->bit_count);
        size_t byte = bit / 8;
        unsigned int offset = (unsigned int)(bit % 8);

        if ((filter->bits[byte] & (uint8_t)(1U << offset)) == 0)
            return 0; /* One zero bit proves this key was not inserted. */
    }

    return 1; /* All bits are set: possible membership, not proof. */
}

void bloom_free(BloomFilter *filter)
{
    if (filter == NULL)
        return;

    free(filter->bits);
    memset(filter, 0, sizeof *filter);
}
