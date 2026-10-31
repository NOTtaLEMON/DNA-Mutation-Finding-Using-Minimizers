#include "top3.h"

#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

#define CHECK(condition)                                                       \
    do {                                                                       \
        if (!(condition)) {                                                    \
            printf("FAIL line %d: %s\n", __LINE__, #condition);               \
            failures++;                                                       \
        }                                                                      \
    } while (0)

static void check_sorted_item(const TopCandidate *items, size_t index,
                              long expected_start, size_t expected_votes)
{
    CHECK(items[index].start == expected_start);
    CHECK(items[index].votes == expected_votes);
}

/* More than three candidates are offered; only the three strongest survive. */
static void test_keeps_top_three(void)
{
    Top3Heap top = {0};
    TopCandidate sorted[TOP3_CAPACITY];
    size_t count = 0;

    CHECK(top3_offer(&top, 10, 2) == 0);
    CHECK(top3_offer(&top, 20, 8) == 0);
    CHECK(top3_offer(&top, 30, 5) == 0);
    CHECK(top3_offer(&top, 40, 1) == 0);
    CHECK(top3_offer(&top, 50, 7) == 0);

    CHECK(top.count == 3);
    CHECK(top3_get_sorted(&top, sorted, &count) == 0);
    CHECK(count == 3);
    check_sorted_item(sorted, 0, 20, 8);
    check_sorted_item(sorted, 1, 50, 7);
    check_sorted_item(sorted, 2, 30, 5);
}

/*
 * The root is the weakest retained candidate. A stronger candidate replaces
 * it; a weaker candidate is discarded without changing the selected set.
 */
static void test_root_replacement_and_discard(void)
{
    Top3Heap top = {0};
    TopCandidate sorted[TOP3_CAPACITY];
    size_t count = 0;

    CHECK(top3_offer(&top, 100, 3) == 0);
    CHECK(top3_offer(&top, 200, 5) == 0);
    CHECK(top3_offer(&top, 300, 4) == 0);
    CHECK(top.count == 3);

    /* New candidate has 6 votes, so it replaces the weakest (3 votes). */
    CHECK(top3_offer(&top, 400, 6) == 0);
    /* New candidate has only 2 votes, so it must not enter the top three. */
    CHECK(top3_offer(&top, 500, 2) == 0);

    CHECK(top3_get_sorted(&top, sorted, &count) == 0);
    CHECK(count == 3);
    check_sorted_item(sorted, 0, 400, 6);
    check_sorted_item(sorted, 1, 200, 5);
    check_sorted_item(sorted, 2, 300, 4);
}

/* Equal vote counts use the documented tie rule: smaller start wins. */
static void test_tie_breaking(void)
{
    Top3Heap top = {0};
    TopCandidate sorted[TOP3_CAPACITY];
    size_t count = 0;

    CHECK(top3_offer(&top, 30, 5) == 0);
    CHECK(top3_offer(&top, 20, 5) == 0);
    CHECK(top3_offer(&top, 40, 5) == 0);
    CHECK(top3_offer(&top, 10, 5) == 0);

    CHECK(top3_get_sorted(&top, sorted, &count) == 0);
    CHECK(count == 3);
    check_sorted_item(sorted, 0, 10, 5);
    check_sorted_item(sorted, 1, 20, 5);
    check_sorted_item(sorted, 2, 30, 5);
}

/* Fewer than three candidates are returned; zero-vote candidates are invalid. */
static void test_fewer_than_three_and_invalid_votes(void)
{
    Top3Heap top = {0};
    TopCandidate sorted[TOP3_CAPACITY];
    size_t count = 99;

    CHECK(top3_offer(&top, -5, 2) == 0);  /* Negative starts are allowed in Stage 5. */
    CHECK(top3_offer(&top, 12, 1) == 0);
    CHECK(top3_offer(&top, 99, 0) == -1);

    CHECK(top3_get_sorted(&top, sorted, &count) == 0);
    CHECK(count == 2);
    check_sorted_item(sorted, 0, -5, 2);
    check_sorted_item(sorted, 1, 12, 1);
}

/* Resetting the heap allows it to be reused for a different read. */
static void test_init_and_reuse(void)
{
    Top3Heap top = {0};
    TopCandidate sorted[TOP3_CAPACITY];
    size_t count = 0;

    CHECK(top3_offer(&top, 10, 3) == 0);
    top3_init(&top);
    CHECK(top.count == 0);

    CHECK(top3_offer(&top, 77, 9) == 0);
    CHECK(top3_get_sorted(&top, sorted, &count) == 0);
    CHECK(count == 1);
    check_sorted_item(sorted, 0, 77, 9);
}


/* Stage 5 VoteTree conversion: Stage 6 should inspect every node and retain the best three. */
static void test_build_from_vote_tree(void)
{
    VoteNode n10 = {10, 9, NULL, NULL};
    VoteNode n25 = {25, 8, NULL, NULL};
    VoteNode n40 = {40, 1, NULL, NULL};
    VoteNode n30 = {30, 7, &n25, &n40};
    VoteNode n20 = {20, 4, &n10, &n30};
    VoteTree tree = {&n20, 6};
    Top3Heap top = {0};
    TopCandidate sorted[TOP3_CAPACITY];
    size_t count = 0;

    CHECK(top3_from_vote_tree(&tree, &top) == 0);
    CHECK(top.count == 3);
    CHECK(top3_get_sorted(&top, sorted, &count) == 0);
    CHECK(count == 3);
    check_sorted_item(sorted, 0, 10, 9);
    check_sorted_item(sorted, 1, 25, 8);
    check_sorted_item(sorted, 2, 30, 7);
}

static void test_empty_vote_tree(void)
{
    VoteTree tree = {NULL, 0};
    Top3Heap top = {{ {0, 0}, {0, 0}, {0, 0} }, 2};
    TopCandidate sorted[TOP3_CAPACITY];
    size_t count = 99;

    CHECK(top3_from_vote_tree(&tree, &top) == 0);
    CHECK(top.count == 0);
    CHECK(top3_get_sorted(&top, sorted, &count) == 0);
    CHECK(count == 0);
}

int main(void)
{
    test_keeps_top_three();
    test_root_replacement_and_discard();
    test_tie_breaking();
    test_fewer_than_three_and_invalid_votes();
    test_init_and_reuse();
    test_build_from_vote_tree();
    test_empty_vote_tree();

    if (failures == 0)
        printf("all Stage 6 top-three heap checks passed\n");
    else
        printf("%d Stage 6 check(s) failed\n", failures);

    return failures != 0;
}
