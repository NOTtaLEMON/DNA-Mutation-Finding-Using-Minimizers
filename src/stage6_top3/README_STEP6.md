# Stage 6: Keep the top three candidate read starts (C)

## Where the files go

Create this folder in the repository and copy the three C files into it:

```text
src/
  stage1_load/
  stage2_minimizers/
  stage3_index/
  stage4_bloom/
  stage5_vote/
  stage6_top3/              <-- new
    top3.c
    top3.h
    test_top3.c
```

This README can sit in `src/stage6_top3/` or in the repository's documentation/material folder. It is not compiled.

## What Stage 6 does

Stage 5 (`vote_read`) processes one read and builds a `VoteTree`. Each tree node represents one distinct candidate reference start and its accumulated vote count.

Stage 6 examines every `(start, votes)` pair in that tree and retains only the three best candidates using a fixed-capacity min-heap.

The min-heap contains at most three items. Its root (`heap[0]`) is deliberately the weakest candidate currently retained:

- fewer votes means a weaker candidate;
- when vote counts tie, the larger start coordinate is treated as weaker.

When a new candidate arrives:

1. If the heap has fewer than three items, insert it and restore the heap property.
2. If the heap already has three items, compare the new candidate with the root.
3. If the new candidate is better than the root, replace the root and restore the heap property.
4. Otherwise discard the new candidate.

The tie rule (higher vote count first, then smaller start coordinate first) is an explicit deterministic project choice. The handwritten Stage 6 notes specify the top-three/min-heap behavior but do not define tie handling.

## Files

- `top3.h`: top-three data structures and function declarations.
- `top3.c`: min-heap implementation and conversion from the Stage 5 vote tree.
- `test_top3.c`: tests for top-three retention, replacement, ties, fewer than three candidates, invalid zero-vote candidates, and reuse.

## How Stage 6 connects to Stage 5

The Stage 5 header provides:

```c
void vote_tree_walk(const VoteTree *tree, VoteVisitFn visit, void *ctx);
```

Stage 6 uses this traversal. It does not need to know how the BST is internally balanced or shaped. It receives each distinct candidate start and its final vote count.

Example integration after Stage 5 has called `vote_read()` successfully:

```c
#include "vote.h"
#include "top3.h"
#include <stdio.h>  /* printf */

VoteTree votes = {0};
VoteStats stats;
Top3Heap top = {0};
TopCandidate ordered[TOP3_CAPACITY];
size_t candidate_count = 0;
size_t i;

/* Stage 5: add votes for this one read. */
if (vote_read(&read, &read_mins, &reference_index,
              &reference_bloom, k, &votes, &stats) != 0) {
    /* Handle an error. */
}

/* Stage 6: retain at most the three candidates with the most votes. */
if (top3_from_vote_tree(&votes, &top) != 0) {
    /* Handle an error. */
}

/* Get the retained candidates in best-first order for Stage 7. */
if (top3_get_sorted(&top, ordered, &candidate_count) != 0) {
    /* Handle an error. */
}

for (i = 0; i < candidate_count; i++) {
    printf("candidate start=%ld, votes=%zu\n",
           ordered[i].start, ordered[i].votes);
}

/* The VoteTree and Top3Heap are per-read state. Clear/reset before next read. */
vote_tree_clear(&votes);
top3_init(&top);
```

The above is an integration example, not a separate file required for Stage 6. The existing `vote_read()` does not clear its input tree, so the tree must be empty before processing a new read. `vote_tree_clear()` is suitable for this.

## API

```c
#define TOP3_CAPACITY 3

void top3_init(Top3Heap *top);
int top3_offer(Top3Heap *top, long start, size_t votes);
int top3_from_vote_tree(const VoteTree *tree, Top3Heap *out);
int top3_get_sorted(const Top3Heap *top,
                    TopCandidate out[TOP3_CAPACITY],
                    size_t *out_count);
```

`top3_from_vote_tree()` resets `out`, walks every vote-tree node, and builds the bounded heap. An empty vote tree is valid and produces zero candidates.

`top3_get_sorted()` copies the retained items into best-first order without changing the heap. This is useful because a heap is not generally stored in fully sorted order; Stage 7 should receive an ordered list.

## Compile and run the unit tests

From the repository root, after the previous stage source files are present:

```sh
cc -std=c11 -Wall -Wextra -Wpedantic \
  -Isrc/stage1_load -Isrc/stage2_minimizers \
  -Isrc/stage3_index -Isrc/stage4_bloom \
  -Isrc/stage5_vote -Isrc/stage6_top3 \
  src/stage1_load/seq.c \
  src/stage2_minimizers/minimizer.c \
  src/stage3_index/index.c \
  src/stage4_bloom/bloom.c \
  src/stage5_vote/vote.c \
  src/stage6_top3/top3.c \
  src/stage6_top3/test_top3.c \
  -lm -o test_top3

./test_top3
```

Adjust the Stage 5 folder name if your repository uses a different name. The Stage 6 unit tests directly exercise the heap. The build includes Stage 1-5 because `top3_from_vote_tree()` uses the actual `vote_tree_walk()` implementation declared in `vote.h`.

## Assumptions and limitations to keep visible

1. **Top-three capacity:** the value 3 comes directly from the project diagram.
2. **Tie handling:** the diagram does not specify ties. This implementation prefers the smaller reference start when vote counts are equal. Change this only if the project rubric specifies another rule.
3. **Vote counts:** Stage 5's VoteTree has one node per distinct start, with the final count at that node. Stage 6 consumes each node once; it must not receive one callback per individual vote.
4. **Negative starts:** they are retained or discarded only according to vote count and tie order. Stage 5 explicitly allows negative starts; the later alignment stage is expected to reject invalid negative reference starts.
5. **Heap order vs output order:** the internal heap is not fully sorted. Use `top3_get_sorted()` to obtain candidates in best-first order.
6. **Per-read state:** create/reset the top-three heap for every read. Clear the Stage 5 VoteTree between reads too.
7. **No edit-distance verification here:** Stage 6 only ranks candidates. It does not align reads, reject bad alignments, or call mutations; those belong to later stages.
8. **BST performance:** Stage 5 currently uses an ordinary, unbalanced BST. Stage 6 visits every vote-tree node, so its total collection work is O(number of distinct candidate starts); the heap work per candidate is O(log 3), effectively constant. A sorted insertion pattern can make the Stage 5 BST tall, but Stage 6 does not change that behavior.
