#include "pileup.h"
#include "verify.h"

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

/* Writes a one-record FASTA with the given letters and loads it through stage 1. */
static int load_letters(const char *letters, Seq *seq)
{
    FILE *f = fopen("test_pileup_tmp.fasta", "wb");
    int rc;

    if (f == NULL)
        return -1;
    fprintf(f, ">s\n%s\n", letters);
    fclose(f);
    rc = seq_load_fasta("test_pileup_tmp.fasta", seq);
    remove("test_pileup_tmp.fasta");
    return rc;
}

/* Loads letters as a read and adds it to the pile at start; returns pileup_add_read's result. */
static int add(Pileup *pile, const char *letters, size_t start)
{
    Seq read;
    int rc;

    if (load_letters(letters, &read) != 0)
        return -2;
    rc = pileup_add_read(pile, &read, start);
    seq_free(&read);
    return rc;
}

/* 3 of 4 reads carry G->T at position 2: one call; the other positions stay quiet. */
static void test_toy_call(void)
{
    Seq ref;
    Pileup pile;
    CallList calls = {0};

    CHECK(load_letters("ACGTACGTAC", &ref) == 0);
    CHECK(pileup_init(&pile, ref.len) == 0);
    CHECK(add(&pile, "ACTT", 0) == 0);
    CHECK(add(&pile, "ACTTA", 0) == 0);
    CHECK(add(&pile, "CTTAC", 1) == 0);
    CHECK(add(&pile, "ACGT", 0) == 0);

    CHECK(pileup_call(&pile, &ref, PILEUP_DEFAULT_MIN_DEPTH,
                      PILEUP_DEFAULT_MIN_ALT_PERCENT, &calls) == 0);
    CHECK(calls.count == 1);
    CHECK(calls.items[0].pos == 2 && calls.items[0].ref == 2 && calls.items[0].alt == 3);
    CHECK(calls.items[0].alt_count == 3 && calls.items[0].depth == 4);

    calls_free(&calls);
    pileup_free(&pile);
    seq_free(&ref);
}

/* One read in four disagrees (a sequencing error): below half, so no call. */
static void test_single_error_not_called(void)
{
    Seq ref;
    Pileup pile;
    CallList calls = {0};

    CHECK(load_letters("ACGTACGTAC", &ref) == 0);
    CHECK(pileup_init(&pile, ref.len) == 0);
    CHECK(add(&pile, "ACGT", 0) == 0);
    CHECK(add(&pile, "ACGT", 0) == 0);
    CHECK(add(&pile, "ACGT", 0) == 0);
    CHECK(add(&pile, "AAGT", 0) == 0);
    CHECK(pileup_call(&pile, &ref, 3, 50, &calls) == 0);
    CHECK(calls.count == 0);

    calls_free(&calls);
    pileup_free(&pile);
    seq_free(&ref);
}

/* Two agreeing reads fail min depth 3 but pass with min depth 2. */
static void test_min_depth(void)
{
    Seq ref;
    Pileup pile;
    CallList calls = {0};

    CHECK(load_letters("ACGTACGTAC", &ref) == 0);
    CHECK(pileup_init(&pile, ref.len) == 0);
    CHECK(add(&pile, "ACTT", 0) == 0);
    CHECK(add(&pile, "ACTT", 0) == 0);
    CHECK(pileup_call(&pile, &ref, 3, 50, &calls) == 0);
    CHECK(calls.count == 0);
    CHECK(pileup_call(&pile, &ref, 2, 50, &calls) == 0);
    CHECK(calls.count == 1 && calls.items[0].pos == 2);

    calls_free(&calls);
    pileup_free(&pile);
    seq_free(&ref);
}

/* A read past the end is clipped, and depth from the difference array equals A+C+G+T everywhere. */
static void test_clip_and_depth(void)
{
    Seq ref;
    Pileup pile;
    uint32_t depth[10];
    size_t p;

    CHECK(load_letters("ACGTACGTAC", &ref) == 0);
    CHECK(pileup_init(&pile, ref.len) == 0);
    CHECK(add(&pile, "ACGTA", 0) == 0);
    CHECK(add(&pile, "GTAC", 2) == 0);
    CHECK(add(&pile, "ACGG", 8) == 0);          /* only "AC" lands on positions 8 and 9 */
    CHECK(pileup_depths(&pile, depth) == 0);

    CHECK(depth[0] == 1 && depth[2] == 2 && depth[4] == 2 && depth[5] == 1);
    CHECK(depth[6] == 0 && depth[8] == 1 && depth[9] == 1);
    for (p = 0; p < 10; p++) {
        const uint32_t *c = pile.cells[p].count;
        CHECK(depth[p] == c[0] + c[1] + c[2] + c[3]);
    }

    pileup_free(&pile);
    seq_free(&ref);
}

/* Bad starts and reads with N are refused; reference N positions are never called. */
static void test_bad_input_and_reference_n(void)
{
    Seq ref;
    Pileup pile;
    CallList calls = {0};

    CHECK(load_letters("ACNNACGTAC", &ref) == 0);
    CHECK(pileup_init(&pile, ref.len) == 0);
    CHECK(add(&pile, "ACGT", 10) == -1);
    CHECK(add(&pile, "ACNT", 0) == -1);
    CHECK(add(&pile, "ACTT", 0) == 0);          /* T at an N position: 3 times */
    CHECK(add(&pile, "ACTT", 0) == 0);
    CHECK(add(&pile, "ACTT", 0) == 0);
    CHECK(pileup_call(&pile, &ref, 3, 50, &calls) == 0);
    CHECK(calls.count == 0);
    CHECK(pile.reads == 3);

    calls_free(&calls);
    pileup_free(&pile);
    seq_free(&ref);
}

/*
 * Cuts len letters at from out of ref as a read. The letter at reference
 * position mut_pos (if inside the read) is changed to simulate the sample's
 * mutation; error_at (read offset, or len or more for none) adds one
 * sequencing error. Loads the result as a Seq.
 */
static int make_read(const Seq *ref, size_t from, size_t len, size_t mut_pos,
                     size_t error_at, Seq *read)
{
    char buf[512];
    size_t i;

    if (len >= sizeof buf)
        return -1;
    for (i = 0; i < len; i++) {
        uint8_t code = seq_get(ref, from + i);
        if (from + i == mut_pos)
            code = (uint8_t)((code + 1) % 4);
        if (i == error_at)
            code = (uint8_t)((code + 2) % 4);
        buf[i] = LETTER[code];
    }
    buf[len] = '\0';
    return load_letters(buf, read);
}

/* Stages 2 to 7 for one read; returns the accepted start, or -1 if the read is not placed. */
static long place_read(const Seq *ref, const Seq *read, const Index *idx,
                       const BloomFilter *bloom, size_t k, size_t w)
{
    MinimizerSet mins = {0};
    VoteTree tree = {0};
    Top3Heap top;
    TopCandidate sorted[TOP3_CAPACITY];
    size_t n = 0;
    VerifyResult vr;
    long placed = -1;

    if (minimizers_generate(read, k, w, &mins) == 0 &&
        vote_read(read, &mins, idx, bloom, k, &tree, NULL) == 0 &&
        top3_from_vote_tree(&tree, &top) == 0 &&
        top3_get_sorted(&top, sorted, &n) == 0 &&
        verify_candidates(ref, read, sorted, n, VERIFY_DEFAULT_BAND,
                          VERIFY_DEFAULT_MAX_EDIT_DISTANCE, &vr) == 0 &&
        vr.mapped)
        placed = vr.start;

    vote_tree_clear(&tree);
    minimizers_free(&mins);
    return placed;
}

/*
 * Stages 1 to 8 on the real reference: 21 reads of 150 letters start every 10
 * letters from 900. The sample has one mutation at 1070; one read also has a
 * sequencing error. Expect exactly one call, at 1070.
 */
static void test_real_reference(const char *path)
{
    const size_t k = 15, w = 10, len = 150, mut_pos = 1070;
    Seq ref, read;
    MinimizerSet ref_mins = {0};
    Index idx;
    BloomFilter bloom = {0};
    Pileup pile;
    CallList calls = {0};
    size_t from, placed = 0, total = 0;

    if (seq_load_fasta(path, &ref) != 0) {
        printf("(could not open %s, skipping)\n", path);
        return;
    }
    CHECK(minimizers_generate(&ref, k, w, &ref_mins) == 0);
    CHECK(index_build(&ref, &ref_mins, k, 0, &idx) == 0);
    CHECK(bloom_build_from_minimizers(&ref, &ref_mins, k, 0.01, &bloom) == 0);
    CHECK(pileup_init(&pile, ref.len) == 0);

    for (from = 900; from <= 1100; from += 10) {
        size_t error_at = (from == 1000) ? 40 : len;    /* one read gets an error at 1040 */
        long start;

        CHECK(make_read(&ref, from, len, mut_pos, error_at, &read) == 0);
        start = place_read(&ref, &read, &idx, &bloom, k, w);
        total++;
        CHECK(start == (long)from);
        if (start >= 0 && pileup_add_read(&pile, &read, (size_t)start) == 0)
            placed++;
        seq_free(&read);
    }

    CHECK(pileup_call(&pile, &ref, PILEUP_DEFAULT_MIN_DEPTH,
                      PILEUP_DEFAULT_MIN_ALT_PERCENT, &calls) == 0);
    printf("%s: %lu of %lu reads placed, %lu call(s)\n", path,
           (unsigned long)placed, (unsigned long)total, (unsigned long)calls.count);
    if (calls.count > 0)
        printf("  first call: position %lu, %c -> %c, %lu of %lu reads\n",
               (unsigned long)calls.items[0].pos, LETTER[calls.items[0].ref],
               LETTER[calls.items[0].alt], (unsigned long)calls.items[0].alt_count,
               (unsigned long)calls.items[0].depth);
    CHECK(calls.count == 1 && calls.items[0].pos == mut_pos);

    calls_free(&calls);
    pileup_free(&pile);
    bloom_free(&bloom);
    index_free(&idx);
    minimizers_free(&ref_mins);
    seq_free(&ref);
}

/* Runs all checks, then stages 1 to 8 on argv[1] (default sequence.fasta). */
int main(int argc, char **argv)
{
    test_toy_call();
    test_single_error_not_called();
    test_min_depth();
    test_clip_and_depth();
    test_bad_input_and_reference_n();
    test_real_reference(argc > 1 ? argv[1] : "sequence.fasta");

    if (failures == 0)
        printf("all checks passed\n");
    else
        printf("%d check(s) failed\n", failures);
    return failures != 0;
}
