#ifndef TOP3_H
#define TOP3_H

#include <stddef.h>

#include "vote.h"

/*
 * Stage 6 keeps at most this many candidate read-start positions.
 * The value 3 comes from the project design/handwritten Stage 6 notes.
 */
#define TOP3_CAPACITY 3

/*
 * One candidate starting position and the number of votes it received in
 * Stage 5. `start` is signed because Stage 5 intentionally allows negative
 * candidate starts; the later alignment-verification stage decides whether
 * a candidate is a valid reference coordinate.
 */
typedef struct {
    long start;
    size_t votes;
} TopCandidate;

/*
 * A fixed-capacity MIN-heap of candidates.
 *
 * IMPORTANT INVARIANT:
 *   heap[0] is always the WORST candidate currently retained:
 *   - fewer votes is worse;
 *   - if vote counts tie, the larger start coordinate is worse.
 *
 * The heap never stores more than TOP3_CAPACITY candidates, so no dynamic
 * allocation is needed for Stage 6.
 */
typedef struct {
    TopCandidate heap[TOP3_CAPACITY];
    size_t count;
} Top3Heap;

/* Reset a heap to the empty state. Safe for a non-NULL heap. */
void top3_init(Top3Heap *top);

/*
 * Offer one already-aggregated candidate to the top-three heap.
 *
 * If fewer than three candidates are retained, the candidate is inserted.
 * Once full, it replaces heap[0] only when it is strictly better than the
 * current weakest retained candidate.
 *
 * Ranking / tie rule (a project choice because the Stage 6 sketch does not
 * specify tie handling):
 *   1. More votes is better.
 *   2. If votes tie, the smaller reference start is better.
 *
 * `votes` must be greater than zero. The input should contain one item per
 * distinct start, with its FINAL Stage 5 vote count (normally obtained by
 * walking the VoteTree). Do not call this once per individual vote.
 *
 * Returns 0 on success, -1 for a NULL heap, corrupt heap count, or zero votes.
 */
int top3_offer(Top3Heap *top, long start, size_t votes);

/*
 * Build a top-three heap from all nodes in one Stage 5 VoteTree.
 * This uses vote_tree_walk(), which visits each distinct start and its
 * accumulated vote count. The resulting heap contains at most three starts.
 *
 * `out` is reset before collecting. An empty VoteTree is a successful input
 * and produces an empty heap. Returns 0 on success, -1 for invalid arguments
 * or an invalid (zero-vote) tree node.
 */
int top3_from_vote_tree(const VoteTree *tree, Top3Heap *out);

/*
 * Copy the retained candidates into `out` in BEST-FIRST order:
 *   highest vote count first; ties ordered by smaller start coordinate.
 *
 * The heap itself is not modified. `out` must have room for
 * TOP3_CAPACITY TopCandidate values. `out_count` receives the number copied.
 * Returns 0 on success, -1 for invalid arguments or a corrupt heap count.
 */
int top3_get_sorted(const Top3Heap *top,
                    TopCandidate out[TOP3_CAPACITY],
                    size_t *out_count);

#endif /* TOP3_H */
