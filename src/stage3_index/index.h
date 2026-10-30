#ifndef INDEX_H
#define INDEX_H

#include <stddef.h>
#include <stdint.h>

#include "seq.h"
#include "minimizer.h"

/* One reference position where a word was kept. Positions form a singly linked list. */
typedef struct PosNode {
    size_t pos;
    struct PosNode *next;
} PosNode;

/*
 * One distinct kept word.
 *
 * word_id: the word's 2-bit codes packed into one integer (first letter in the highest bits).
 *          Two words are equal exactly when their IDs are equal, so this is the real identity check.
 * hash:    the stage 2 hash. Picks the bucket and is a quick first filter, but is NOT proof of identity.
 */
typedef struct IndexEntry {
    uint64_t word_id;
    uint64_t hash;
    PosNode *pos_head;           /* positions in the order they were added (increasing) */
    PosNode *pos_tail;           /* last position node, so appending is O(1) */
    size_t npos;
    struct IndexEntry *next;     /* next entry in the same bucket's chain */
} IndexEntry;

/* Hash table with chaining: buckets[b] is the head of a linked list of entries. */
typedef struct {
    IndexEntry **buckets;
    size_t nbuckets;
    size_t nentries;             /* number of distinct words */
    size_t k;                    /* word length the index was built with */
} Index;

/*
 * Packs the k letters starting at pos into one integer, 2 bits per letter.
 * Example: AAGC -> 00 00 10 01 -> 9. The caller guarantees k <= 32 and pos + k <= seq->len.
 */
uint64_t word_id_at(const Seq *seq, size_t pos, size_t k);

/*
 * Builds the index from the reference's kept words.
 *
 *   seq       the reference from stage 1 (used to read each word's letters)
 *   mins      the reference's minimizers from stage 2
 *   k         word length used in stage 2; must be 1..32
 *   nbuckets  0 = size from the input (about mins->count / 0.75); any other value is used as given
 *   out       filled in on success
 *
 * Words already in the table get the new position appended; new words are added at the end of their
 * chain, so each chain is in first-seen order. Returns 0 on success, -1 on bad arguments (including a
 * minimizer position too close to the end of seq) or if memory runs out; out is left empty on failure.
 */
int index_build(const Seq *seq, const MinimizerSet *mins, size_t k, size_t nbuckets, Index *out);

/*
 * Finds the entry for a word, given its stage 2 hash and its word ID.
 * Returns the entry, or NULL if the word is not in the index.
 * If probes is not NULL, it receives the number of entries checked (for benchmarking chain length).
 */
const IndexEntry *index_lookup(const Index *idx, uint64_t hash, uint64_t word_id, size_t *probes);

/* Frees every entry, every position node and the bucket array, then resets idx to empty. */
void index_free(Index *idx);

#endif /* INDEX_H */
