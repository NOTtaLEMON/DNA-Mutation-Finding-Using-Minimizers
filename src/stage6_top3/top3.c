#include "top3.h"

#include <string.h>

/*
 * Return non-zero when `a` should rank BELOW `b` (i.e. `a` is worse).
 * This ordering makes the weakest retained candidate the min-heap root.
 *
 * More votes is better. For equal votes, a smaller start is preferred.
 */
static int candidate_is_worse(const TopCandidate *a, const TopCandidate *b)
{
    if (a->votes != b->votes)
        return a->votes < b->votes;

    return a->start > b->start;
}

/* The inverse ranking: `a` is better than `b`. */
static int candidate_is_better(const TopCandidate *a, const TopCandidate *b)
{
    return candidate_is_worse(b, a);
}

static void swap_candidates(TopCandidate *a, TopCandidate *b)
{
    TopCandidate temp = *a;
    *a = *b;
    *b = temp;
}

/*
 * Restore the min-heap property after appending a new item.
 * A worse child must move toward the root, because the root is the weakest
 * candidate currently retained.
 */
static void sift_up(Top3Heap *top, size_t child_index)
{
    while (child_index > 0) {
        size_t parent_index = (child_index - 1) / 2;

        if (!candidate_is_worse(&top->heap[child_index],
                                &top->heap[parent_index]))
            break;

        swap_candidates(&top->heap[child_index], &top->heap[parent_index]);
        child_index = parent_index;
    }
}

/*
 * Restore the min-heap property after replacing the root.
 * Move the worse of the two children upward until the weakest candidate is
 * again at heap[0].
 */
static void sift_down(Top3Heap *top, size_t parent_index)
{
    for (;;) {
        size_t left = parent_index * 2 + 1;
        size_t right = left + 1;
        size_t worst = parent_index;

        if (left < top->count &&
            candidate_is_worse(&top->heap[left], &top->heap[worst]))
            worst = left;

        if (right < top->count &&
            candidate_is_worse(&top->heap[right], &top->heap[worst]))
            worst = right;

        if (worst == parent_index)
            break;

        swap_candidates(&top->heap[parent_index], &top->heap[worst]);
        parent_index = worst;
    }
}

void top3_init(Top3Heap *top)
{
    if (top == NULL)
        return;

    memset(top, 0, sizeof *top);
}

int top3_offer(Top3Heap *top, long start, size_t votes)
{
    TopCandidate candidate;

    if (top == NULL || top->count > TOP3_CAPACITY || votes == 0)
        return -1;

    candidate.start = start;
    candidate.votes = votes;

    /* Until three candidates have been collected, insert and heapify upward. */
    if (top->count < TOP3_CAPACITY) {
        size_t inserted_at = top->count;
        top->heap[inserted_at] = candidate;
        top->count++;
        sift_up(top, inserted_at);
        return 0;
    }

    /*
     * The heap is full. heap[0] is its weakest retained candidate.
     * A new candidate matters only if it ranks strictly better than that root.
     */
    if (candidate_is_better(&candidate, &top->heap[0])) {
        top->heap[0] = candidate;
        sift_down(top, 0);
    }

    /* If it did not beat the root, it is safely discarded. */
    return 0;
}

/* Context passed through VoteVisitFn while converting a VoteTree into a heap. */
typedef struct {
    Top3Heap *top;
    int status;
} CollectContext;

/*
 * Callback invoked once for each VoteTree node.
 * A VoteNode already aggregates all votes for one distinct start, so this
 * callback offers that final (start, count) pair exactly once.
 */
static void collect_vote_node(long start, size_t count, void *ctx)
{
    CollectContext *collect = (CollectContext *)ctx;

    if (collect == NULL || collect->top == NULL || collect->status != 0)
        return;

    if (top3_offer(collect->top, start, count) != 0)
        collect->status = -1;
}

int top3_from_vote_tree(const VoteTree *tree, Top3Heap *out)
{
    CollectContext collect;

    if (out == NULL)
        return -1;

    top3_init(out);

    if (tree == NULL)
        return -1;

    collect.top = out;
    collect.status = 0;

    /* Stage 5 provides this in-order traversal callback interface. */
    vote_tree_walk(tree, collect_vote_node, &collect);

    if (collect.status != 0) {
        top3_init(out);
        return -1;
    }

    return 0;
}

int top3_get_sorted(const Top3Heap *top,
                    TopCandidate out[TOP3_CAPACITY],
                    size_t *out_count)
{
    size_t i;

    if (top == NULL || out == NULL || out_count == NULL ||
        top->count > TOP3_CAPACITY)
        return -1;

    *out_count = top->count;

    /* Copy first: sorting the presentation order must not damage the heap. */
    for (i = 0; i < top->count; i++)
        out[i] = top->heap[i];

    /*
     * Insertion-sort at most three items into best-first order.
     * This tiny final sort is simpler and safer than changing the heap layout.
     */
    for (i = 1; i < *out_count; i++) {
        TopCandidate value = out[i];
        size_t j = i;

        while (j > 0 && candidate_is_better(&value, &out[j - 1])) {
            out[j] = out[j - 1];
            j--;
        }
        out[j] = value;
    }

    return 0;
}
