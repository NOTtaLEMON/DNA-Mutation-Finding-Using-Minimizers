#ifndef VOTE_H
#define VOTE_H

#include <stddef.h>
#include <stdint.h>

#include "seq.h"
#include "minimizer.h"
#include "index.h"
#include "bloom.h"

/*
 * One node of the vote tally: a binary search tree keyed by candidate start.
 *
 * start: candidate start of the read on the reference = p - read_pos, where p is the reference
 *        position of a matching word and read_pos is that word's position inside the read.
 *        Signed, because a chance match can give a negative start (stage 7 rejects those).
 * count: number of read words that voted for this start.
 */
typedef struct VoteNode {
    long start;
    size_t count;
    struct VoteNode *left;     /* starts smaller than this one */
    struct VoteNode *right;    /* starts larger than this one */
} VoteNode;

/* The tally for one read. Empty when root is NULL. */
typedef struct {
    VoteNode *root;
    size_t nodes;              /* number of distinct starts */
} VoteTree;

/* Counters for one read, for the benchmark. */
typedef struct {
    size_t words;              /* read words (minimizers) considered */
    size_t bloom_skipped;      /* words the Bloom filter reported definitely absent */
    size_t lookups;            /* words looked up in the index */
    size_t hits;               /* reference positions returned, i.e. votes cast */
} VoteStats;

/* Called once per node by vote_tree_walk, in increasing order of start. Stage 6's insert fits this. */
typedef void (*VoteVisitFn)(long start, size_t count, void *ctx);

/*
 * Adds one vote for start: walks down from the root (smaller -> left, larger -> right).
 * If start is already in the tree its count goes up by 1; otherwise a new node with count 1 is added
 * at the empty spot reached. Returns 0 on success, -1 if memory ran out.
 */
int vote_tree_add(VoteTree *tree, long start);

/* In-order walk: calls visit(start, count, ctx) for every node, smallest start first. */
void vote_tree_walk(const VoteTree *tree, VoteVisitFn visit, void *ctx);

/* Frees every node, children before parent, and leaves the tree empty for the next read. */
void vote_tree_clear(VoteTree *tree);

/* Number of nodes on the longest root-to-leaf path (0 for an empty tree). */
size_t vote_tree_height(const VoteTree *tree);

/*
 * Stage 5 for one read: for each read minimizer, ask the Bloom filter; if it says "maybe", look the
 * word up in the index; for every reference position p returned, add a vote for p - read_pos.
 *
 *   read       the read, loaded like any sequence
 *   read_mins  the read's minimizers from stage 2 (same k and window as the reference)
 *   idx        the reference index from stage 3 (built with the same k)
 *   bloom      the reference Bloom filter from stage 4 (required)
 *   k          word length
 *   tree       tally to add to; normally empty (call vote_tree_clear between reads)
 *   stats      if not NULL, filled with this read's counters
 *
 * Returns 0 on success, -1 on bad arguments, a broken Bloom filter, or if memory ran out.
 */
int vote_read(const Seq *read, const MinimizerSet *read_mins, const Index *idx,
              const BloomFilter *bloom, size_t k, VoteTree *tree, VoteStats *stats);

#endif /* VOTE_H */
