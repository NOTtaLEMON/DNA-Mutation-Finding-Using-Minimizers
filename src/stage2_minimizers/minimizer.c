#include "minimizer.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* A candidate stored in the monotonic deque used for sliding-window minima. */
typedef struct {
    size_t pos;       /* start position of this k-mer in the DNA sequence */
    uint64_t hash;    /* hash used to rank this k-mer */
} DequeEntry;

/*
 * Hash one k-mer using 64-bit FNV-1a.
 *
 * The hash is deterministic across runs and platforms that use uint64_t.
 * The DNA is read through seq_get(), so this function relies on Stage 1's
 * seq_get() correctly decoding the two-bit packed representation.
 *
 * IMPORTANT: A 64-bit hash can theoretically collide. Later index lookup
 * should use the hash to find candidates and, where correctness requires it,
 * compare the actual k-mer sequence before treating two keys as identical.
 */
static uint64_t hash_kmer(const Seq *seq, size_t pos, size_t k)
{
    const uint64_t fnv_offset_basis = UINT64_C(14695981039346656037);
    const uint64_t fnv_prime = UINT64_C(1099511628211);
    uint64_t hash = fnv_offset_basis;
    size_t i;

    for (i = 0; i < k; i++) {
        /* seq_get returns A=0, C=1, G=2, T=3. */
        uint8_t code = seq_get(seq, pos + i);
        hash ^= (uint64_t)code;
        hash *= fnv_prime; /* unsigned overflow is defined modulo 2^64 */
    }

    return hash;
}

/*
 * Return non-zero if the k-mer [pos, pos + k) overlaps any unknown-base run.
 *
 * Seq.runs is sorted by start and runs do not overlap. Binary search finds
 * the first run whose exclusive end is after `pos`; only that run can be the
 * first one overlapping this k-mer.
 */
static int kmer_overlaps_unknown(const Seq *seq, size_t pos, size_t k)
{
    size_t low = 0;
    size_t high = seq->nruns;
    const size_t kmer_end = pos + k; /* safe: caller only passes valid k-mer starts */

    while (low < high) {
        size_t mid = low + (high - low) / 2;
        size_t run_end = seq->runs[mid].start + seq->runs[mid].length;

        if (run_end <= pos)
            low = mid + 1;
        else
            high = mid;
    }

    return low < seq->nruns && seq->runs[low].start < kmer_end;
}

/* Append one result, growing the output array when necessary. */
static int append_minimizer(MinimizerSet *set, uint64_t hash, size_t pos)
{
    if (set->count > 0 && set->items[set->count - 1].pos == pos) {
        /* The same occurrence won consecutive windows; store it only once. */
        return 0;
    }

    if (set->count == set->capacity) {
        size_t new_capacity;
        Minimizer *new_items;

        if (set->capacity == 0) {
            new_capacity = 16;
        } else {
            if (set->capacity > SIZE_MAX / 2)
                return -1;
            new_capacity = set->capacity * 2;
        }

        if (new_capacity > SIZE_MAX / sizeof *new_items)
            return -1;

        new_items = (Minimizer *)realloc(set->items,
                                         new_capacity * sizeof *new_items);
        if (new_items == NULL)
            return -1;

        set->items = new_items;
        set->capacity = new_capacity;
    }

    set->items[set->count].hash = hash;
    set->items[set->count].pos = pos;
    set->count++;
    return 0;
}

/*
 * Process one contiguous run of valid k-mers [run_start, run_end).
 * All positions are k-mer start positions in the original sequence.
 *
 * The deque is a monotonic queue: hashes increase from its front to its back.
 * Therefore the front always holds the minimum hash in the current window.
 */
static int process_valid_run(const Seq *seq,
                             size_t k,
                             size_t window_kmers,
                             size_t run_start,
                             size_t run_end,
                             DequeEntry *deque,
                             size_t deque_capacity,
                             MinimizerSet *out)
{
    const size_t run_length = run_end - run_start;
    size_t pos;

    if (run_length == 0)
        return 0;

    /*
     * Short-segment policy: if there is no complete window, retain the
     * single minimum from the available valid k-mers. See header comment.
     */
    if (run_length < window_kmers) {
        size_t best_pos = run_start;
        uint64_t best_hash = hash_kmer(seq, run_start, k);

        for (pos = run_start + 1; pos < run_end; pos++) {
            uint64_t hash = hash_kmer(seq, pos, k);
            /* Strict > keeps the leftmost occurrence in a tie. */
            if (hash < best_hash) {
                best_hash = hash;
                best_pos = pos;
            }
        }

        return append_minimizer(out, best_hash, best_pos);
    }

    /* A full window is possible, so the deque must hold up to window_kmers. */
    if (deque == NULL || deque_capacity < window_kmers)
        return -1;

    {
        size_t head = 0;
        size_t count = 0;

        for (pos = run_start; pos < run_end; pos++) {
            uint64_t hash = hash_kmer(seq, pos, k);

            /*
             * 1. Remove entries that are now left of the current window.
             * A window ending at pos begins at pos - window_kmers + 1.
             */
            while (count > 0 && pos - deque[head].pos >= window_kmers) {
                head = (head + 1) % deque_capacity;
                count--;
            }

            /*
             * 2. Remove larger hashes from the back. For equal hashes we
             *    keep the older/leftmost entry, giving deterministic ties.
             */
            while (count > 0) {
                size_t back = (head + count - 1) % deque_capacity;
                if (deque[back].hash <= hash)
                    break;
                count--;
            }

            /* 3. Add this k-mer at the back of the deque. */
            {
                size_t tail = (head + count) % deque_capacity;
                deque[tail].pos = pos;
                deque[tail].hash = hash;
                count++;
            }

            /*
             * 4. Once a complete window has been seen, the front is its
             *    minimum. Append it; append_minimizer removes repeated
             *    selections of the same occurrence in adjacent windows.
             */
            if (pos - run_start + 1 >= window_kmers) {
                if (append_minimizer(out, deque[head].hash, deque[head].pos) != 0)
                    return -1;
            }
        }
    }

    return 0;
}

int minimizers_generate(const Seq *seq,
                        size_t k,
                        size_t window_kmers,
                        MinimizerSet *out)
{
    MinimizerSet result = {0};
    size_t number_of_kmers;
    size_t pos;
    size_t deque_capacity;
    DequeEntry *deque = NULL;

    if (seq == NULL || out == NULL || k == 0 || window_kmers == 0)
        return -1;

    /* Validate the relevant parts of the Seq structure before using them. */
    if ((seq->len > 0 && seq->cells == NULL) ||
        (seq->nruns > 0 && seq->runs == NULL))
        return -1;

    /*
     * The caller must pass an empty output set. This avoids silently leaking
     * an earlier allocation if a non-empty set is accidentally reused.
     */
    if (out->items != NULL || out->count != 0 || out->capacity != 0)
        return -1;

    if (k > seq->len) {
        *out = result; /* valid request, but there are no complete k-mers */
        return 0;
    }

    number_of_kmers = seq->len - k + 1;

    /* If window_kmers is larger than all available k-mers, no full window
       exists. The short-segment rule will still emit one minimum per valid run. */
    deque_capacity = window_kmers < number_of_kmers
                         ? window_kmers
                         : number_of_kmers;

    if (window_kmers <= number_of_kmers) {
        if (deque_capacity > SIZE_MAX / sizeof *deque)
            return -1;
        deque = (DequeEntry *)malloc(deque_capacity * sizeof *deque);
        if (deque == NULL)
            return -1;
    }

    /*
     * Walk all k-mer starts, splitting them into contiguous valid runs.
     * A k-mer overlapping any Stage-1 N-run is not allowed to become a
     * minimizer. No window is allowed to bridge across an unknown base.
     */
    pos = 0;
    while (pos < number_of_kmers) {
        size_t run_start;

        /* Skip invalid k-mers until a valid run begins. */
        while (pos < number_of_kmers && kmer_overlaps_unknown(seq, pos, k))
            pos++;

        if (pos == number_of_kmers)
            break;

        run_start = pos;

        /* Find the end of this contiguous valid run. */
        while (pos < number_of_kmers && !kmer_overlaps_unknown(seq, pos, k))
            pos++;

        if (process_valid_run(seq,
                              k,
                              window_kmers,
                              run_start,
                              pos,
                              deque,
                              deque_capacity,
                              &result) != 0) {
            free(deque);
            minimizers_free(&result);
            return -1;
        }
    }

    free(deque);
    *out = result;
    return 0;
}

void minimizers_free(MinimizerSet *set)
{
    if (set == NULL)
        return;

    free(set->items);
    memset(set, 0, sizeof *set);
}
