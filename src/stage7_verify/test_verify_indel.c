/*
 * INACTIVE: tests for verify_indel.c. Compiled only with -DVERIFY_INDEL_ENABLED.
 *
 * Build (from the repository root):
 *   gcc -std=c99 -Wall -Wextra -DVERIFY_INDEL_ENABLED \
 *     -Isrc/stage1_load -Isrc/stage2_minimizers -Isrc/stage3_index \
 *     -Isrc/stage4_bloom -Isrc/stage5_vote -Isrc/stage6_top3 -Isrc/stage7_verify \
 *     src/stage1_load/seq.c src/stage7_verify/verify.c \
 *     src/stage7_verify/verify_indel.c src/stage7_verify/test_verify_indel.c \
 *     -o test_verify_indel
 */
#ifdef VERIFY_INDEL_ENABLED

#include "verify_indel.h"

#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition)                                                        \
    do {                                                                        \
        if (!(condition)) {                                                     \
            printf("FAIL line %d: %s\n", __LINE__, #condition);                \
            failures++;                                                         \
        }                                                                       \
    } while (0)

/* The reference used by every test (15 letters). */
static const char REFERENCE[] = ">ref\nGATCCTAGGACTTGA\n";

/* Write a tiny FASTA record and load it through the real Stage 1 loader. */
static int load_text(const char *path, const char *text, Seq *seq)
{
    FILE *file = fopen(path, "wb");

    if (file == NULL)
        return -1;
    fputs(text, file);
    fclose(file);
    return seq_load_fasta(path, seq);
}

/* Count how many times op appears in an alignment path. */
static size_t count_ops(const IndelCandidateResult *r, char op)
{
    size_t i, n = 0;

    for (i = 0; i < r->nops; i++)
        n += r->ops[i] == op;
    return n;
}

/*
 * Run both Stage 7 versions on one read placed at start 0 and check the
 * indel version's distance, path and reference span. Also checks whether
 * the original verify.c accepts it, to show the difference.
 */
static void check_read(const char *label, const char *read_text,
                       size_t want_distance, size_t want_i, size_t want_d,
                       size_t want_span, int old_accepts)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[1] = {{0, 10}};
    VerifyResult old_result;
    IndelVerifyResult result;
    const IndelCandidateResult *best;

    CHECK(load_text("test_indel_ref.fa", REFERENCE, &reference) == 0);
    CHECK(load_text("test_indel_read.fa", read_text, &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 3, &old_result) == 0);
    CHECK(old_result.mapped == old_accepts);

    CHECK(verify_candidates_indel(&reference, &read, candidates, 1, 2, 3, &result) == 0);
    CHECK(result.mapped == 1);
    best = &result.candidates[result.best];
    CHECK(best->edit_distance == want_distance);
    CHECK(count_ops(best, ALIGN_INSERT) == want_i);
    CHECK(count_ops(best, ALIGN_DELETE) == want_d);
    CHECK(best->ref_span == want_span);
    CHECK(count_ops(best, ALIGN_MATCH) + count_ops(best, ALIGN_INSERT) == read.len);

    printf("  %-28s old: %s, distance %lu | indel: distance %lu, path %.*s\n",
           label, old_result.mapped ? "accepted" : "rejected",
           (unsigned long)old_result.candidates[0].edit_distance,
           (unsigned long)best->edit_distance, (int)best->nops, best->ops);

    verify_indel_result_free(&result);
    seq_free(&reference);
    seq_free(&read);
    remove("test_indel_ref.fa");
    remove("test_indel_read.fa");
}

/* A read that cannot fit within the threshold is rejected and has no path. */
static void test_rejected_read_has_no_path(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[1] = {{0, 10}};
    IndelVerifyResult result;

    CHECK(load_text("test_indel_ref.fa", REFERENCE, &reference) == 0);
    CHECK(load_text("test_indel_read.fa", ">read\nTTTTTTTTTTTT\n", &read) == 0);

    CHECK(verify_candidates_indel(&reference, &read, candidates, 1, 2, 3, &result) == 0);
    CHECK(result.mapped == 0);
    CHECK(result.candidates[0].accepted == 0);
    CHECK(result.candidates[0].ops == NULL);

    verify_indel_result_free(&result);
    seq_free(&reference);
    seq_free(&read);
    remove("test_indel_ref.fa");
    remove("test_indel_read.fa");
}

/* Runs every check and prints a side-by-side comparison with verify.c. */
int main(void)
{
    /* ref GATCCTAGGACTTGA; each read is 12 letters placed at start 0. */
    check_read("exact copy", ">read\nGATCCTAGGACT\n", 0, 0, 0, 12, 1);
    check_read("1 extra letter", ">read\nGATCCATAGGAC\n", 1, 1, 0, 11, 1);
    check_read("1 missing letter", ">read\nGATCCAGGACTT\n", 1, 0, 1, 13, 1);
    check_read("2 changes + 1 extra letter", ">read\nGTTCCATAGCAC\n", 3, 1, 0, 11, 0);
    test_rejected_read_has_no_path();

    if (failures == 0)
        printf("all indel Stage 7 tests passed\n");
    else
        printf("%d indel Stage 7 test(s) failed\n", failures);
    return failures != 0;
}

#endif /* VERIFY_INDEL_ENABLED */
