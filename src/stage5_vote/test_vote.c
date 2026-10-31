#include "vote.h"

#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("FAIL line %d: %s\n", __LINE__, #cond);                   \
            failures++;                                                      \
        }                                                                    \
    } while (0)

static const char LETTER[] = "ACGT";

/* Pairs collected by a walk, so tests can check the order and counts. */
typedef struct {
    long start[64];
    size_t count[64];
    size_t n;
} Collected;

/* Walk callback for tests: records each (start, count) in a Collected. */
static void collect(long start, size_t count, void *ctx)
{
    Collected *c = (Collected *)ctx;
    if (c->n < 64) {
        c->start[c->n] = start;
        c->count[c->n] = count;
        c->n++;
    }
}

/* Walk callback that remembers the start with the most votes (first one wins a tie). */
static void track_best(long start, size_t count, void *ctx)
{
    long *best = (long *)ctx;          /* best[0] = start, best[1] = count */
    if ((size_t)best[1] < count) {
        best[0] = start;
        best[1] = (long)count;
    }
}

/* The toy read #4 votes 3, 3, 26, 3, 3, 15, 3, 38: root 3 (count 5), right 26, whose children are 15 and 38. */
static void test_toy_votes(void)
{
    static const long votes[] = {3, 3, 26, 3, 3, 15, 3, 38};
    VoteTree tree = {0};
    Collected c = {{0}, {0}, 0};
    size_t i;

    for (i = 0; i < sizeof votes / sizeof votes[0]; i++)
        CHECK(vote_tree_add(&tree, votes[i]) == 0);

    CHECK(tree.nodes == 4);
    CHECK(tree.root->start == 3 && tree.root->count == 5);
    CHECK(tree.root->left == NULL);
    CHECK(tree.root->right->start == 26);
    CHECK(tree.root->right->left->start == 15);
    CHECK(tree.root->right->right->start == 38);
    CHECK(vote_tree_height(&tree) == 3);

    vote_tree_walk(&tree, collect, &c);       /* in-order: sorted by start */
    CHECK(c.n == 4);
    CHECK(c.start[0] == 3 && c.count[0] == 5);
    CHECK(c.start[1] == 15 && c.count[1] == 1);
    CHECK(c.start[2] == 26 && c.count[2] == 1);
    CHECK(c.start[3] == 38 && c.count[3] == 1);

    vote_tree_clear(&tree);
    CHECK(tree.root == NULL && tree.nodes == 0 && vote_tree_height(&tree) == 0);
}

/* Negative starts are stored like any other key, and a cleared tree can be reused. */
static void test_negative_and_reuse(void)
{
    VoteTree tree = {0};
    Collected c = {{0}, {0}, 0};

    CHECK(vote_tree_add(&tree, 5) == 0);
    CHECK(vote_tree_add(&tree, -6) == 0);
    CHECK(vote_tree_add(&tree, -6) == 0);
    vote_tree_walk(&tree, collect, &c);
    CHECK(c.n == 2 && c.start[0] == -6 && c.count[0] == 2 && c.start[1] == 5);
    vote_tree_clear(&tree);

    CHECK(vote_tree_add(&tree, 42) == 0);
    CHECK(tree.nodes == 1 && tree.root->start == 42 && tree.root->count == 1);
    vote_tree_clear(&tree);
}

/* Sorted starts turn a plain BST into a chain: height equals the number of nodes. */
static void test_sorted_degenerates(void)
{
    VoteTree tree = {0};
    long s;

    for (s = 0; s < 20; s++)
        CHECK(vote_tree_add(&tree, s) == 0);
    CHECK(vote_tree_height(&tree) == 20);
    vote_tree_clear(&tree);
}

/*
 * Cuts len letters at 'from' out of ref, optionally changes the letter at offset 'change_at'
 * (pass len or more for no change), writes them as a one-record FASTA and loads them as a read.
 */
static int make_read(const Seq *ref, size_t from, size_t len, size_t change_at, Seq *read)
{
    FILE *f = fopen("test_vote_read.fasta", "wb");
    size_t i;

    if (f == NULL)
        return -1;
    fputs(">read\n", f);
    for (i = 0; i < len; i++) {
        uint8_t code = seq_get(ref, from + i);
        if (i == change_at)
            code = (uint8_t)((code + 1) % 4);      /* a different letter: a planted mutation */
        fputc(LETTER[code], f);
    }
    fputc('\n', f);
    fclose(f);
    return seq_load_fasta("test_vote_read.fasta", read);
}

/* Votes one read and prints the winner; returns the start with the most votes. */
static long vote_and_report(const char *label, const Seq *read, const Index *idx,
                            const BloomFilter *bloom, size_t k, size_t w)
{
    MinimizerSet mins = {0};
    VoteTree tree = {0};
    VoteStats st;
    long best[2] = {0, 0};

    CHECK(minimizers_generate(read, k, w, &mins) == 0);
    CHECK(vote_read(read, &mins, idx, bloom, k, &tree, &st) == 0);
    vote_tree_walk(&tree, track_best, best);

    printf("  %s: %lu words, %lu skipped by Bloom, %lu lookups, %lu votes, %lu distinct starts, "
           "tree height %lu, winner start %ld with %ld votes\n",
           label, (unsigned long)st.words, (unsigned long)st.bloom_skipped,
           (unsigned long)st.lookups, (unsigned long)st.hits, (unsigned long)tree.nodes,
           (unsigned long)vote_tree_height(&tree), best[0], best[1]);

    vote_tree_clear(&tree);
    minimizers_free(&mins);
    return best[0];
}

/* Stages 1 to 5 on the real reference, with reads cut from position 1000. */
static void test_real_reference(const char *path)
{
    const size_t k = 15, w = 10, from = 1000, len = 150;
    Seq ref, read;
    MinimizerSet ref_mins = {0};
    Index idx;
    BloomFilter bloom = {0};

    if (seq_load_fasta(path, &ref) != 0) {
        printf("(could not open %s, skipping)\n", path);
        return;
    }
    CHECK(minimizers_generate(&ref, k, w, &ref_mins) == 0);
    CHECK(index_build(&ref, &ref_mins, k, 0, &idx) == 0);
    CHECK(bloom_build_from_minimizers(&ref, &ref_mins, k, 0.01, &bloom) == 0);
    printf("%s, k=%lu, w=%lu, read of %lu letters cut from position %lu:\n", path,
           (unsigned long)k, (unsigned long)w, (unsigned long)len, (unsigned long)from);

    CHECK(make_read(&ref, from, len, len, &read) == 0);
    CHECK(vote_and_report("exact copy", &read, &idx, &bloom, k, w) == (long)from);
    seq_free(&read);

    CHECK(make_read(&ref, from, len, 70, &read) == 0);       /* one changed letter in the middle */
    CHECK(vote_and_report("one change at offset 70", &read, &idx, &bloom, k, w) == (long)from);
    seq_free(&read);
    remove("test_vote_read.fasta");

    bloom_free(&bloom);
    index_free(&idx);
    minimizers_free(&ref_mins);
    seq_free(&ref);
}

/* Runs all checks, then stages 1 to 5 on argv[1] (default sequence.fasta). */
int main(int argc, char **argv)
{
    test_toy_votes();
    test_negative_and_reuse();
    test_sorted_degenerates();
    test_real_reference(argc > 1 ? argv[1] : "sequence.fasta");

    if (failures == 0)
        printf("all checks passed\n");
    else
        printf("%d check(s) failed\n", failures);
    return failures != 0;
}
