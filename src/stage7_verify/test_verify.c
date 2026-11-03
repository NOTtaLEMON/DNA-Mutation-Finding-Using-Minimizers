#include "verify.h"

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

/* Write a tiny FASTA record and load it through the real Stage 1 loader. */
static int load_text(const char *path, const char *text, Seq *seq)
{
    FILE *file = fopen(path, "wb");

    if (file == NULL)
        return -1;
    if (fwrite(text, 1, strlen(text), file) != strlen(text)) {
        fclose(file);
        remove(path);
        return -1;
    }
    if (fclose(file) != 0) {
        remove(path);
        return -1;
    }
    return seq_load_fasta(path, seq);
}

static void test_exact_match_and_best_distance(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {
        {4, 20}, /* ACGT vs ACGA: distance 1, more votes */
        {0, 5},  /* ACGT vs ACGT: distance 0, fewer votes */
        {0, 0}
    };
    VerifyResult result;

    CHECK(load_text("test_verify_ref_a.fa", ">ref\nACGTACGA\n", &reference) == 0);
    CHECK(load_text("test_verify_read_a.fa", ">read\nACGT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 2, 2, 3, &result) == 0);
    CHECK(result.candidate_count == 2);
    CHECK(result.candidates[0].edit_distance == 1);
    CHECK(result.candidates[0].accepted == 1);
    CHECK(result.candidates[1].edit_distance == 0);
    CHECK(result.candidates[1].accepted == 1);
    CHECK(result.mapped == 1);
    CHECK(result.start == 0); /* Smaller distance wins, not the larger vote count. */
    CHECK(result.votes == 5);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_a.fa");
    remove("test_verify_read_a.fa");
}

static void test_substitution_distance(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{0, 7}, {0, 0}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_b.fa", ">ref\nAGGT\n", &reference) == 0);
    CHECK(load_text("test_verify_read_b.fa", ">read\nACGT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 3, &result) == 0);
    CHECK(result.candidates[0].edit_distance == 1);
    CHECK(result.candidates[0].accepted == 1);
    CHECK(result.mapped == 1 && result.start == 0);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_b.fa");
    remove("test_verify_read_b.fa");
}

static void test_insertion_or_deletion_cost(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{0, 4}, {0, 0}, {0, 0}};
    VerifyResult result;

    /* Read has one extra G compared with the reference window. */
    CHECK(load_text("test_verify_ref_c.fa", ">ref\nACGT\n", &reference) == 0);
    CHECK(load_text("test_verify_read_c.fa", ">read\nACGGT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 3, &result) == 0);
    CHECK(result.candidates[0].edit_distance == 1);
    CHECK(result.candidates[0].accepted == 1);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_c.fa");
    remove("test_verify_read_c.fa");
}

static void test_threshold_rejection(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{0, 9}, {0, 0}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_d.fa", ">ref\nACGT\n", &reference) == 0);
    CHECK(load_text("test_verify_read_d.fa", ">read\nTTTT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 2, &result) == 0);
    CHECK(result.candidates[0].edit_distance == 3);
    CHECK(result.candidates[0].accepted == 0);
    CHECK(result.mapped == 0);
    CHECK(result.edit_distance == VERIFY_DISTANCE_INF);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_d.fa");
    remove("test_verify_read_d.fa");
}

static void test_band_length_difference(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{0, 3}, {0, 0}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_e.fa", ">ref\nA\n", &reference) == 0);
    CHECK(load_text("test_verify_read_e.fa", ">read\nACGT\n", &read) == 0);

    /* Length difference is 3, greater than band 2: no path can reach the end. */
    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 3, &result) == 0);
    CHECK(result.candidates[0].edit_distance == VERIFY_DISTANCE_INF);
    CHECK(result.candidates[0].status == VERIFY_CANDIDATE_OUTSIDE_BAND);
    CHECK(result.candidates[0].accepted == 0);
    CHECK(result.mapped == 0);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_e.fa");
    remove("test_verify_read_e.fa");
}

static void test_reference_right_edge_is_clipped(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{4, 6}, {0, 0}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_f.fa", ">ref\nACGTAC\n", &reference) == 0);
    CHECK(load_text("test_verify_read_f.fa", ">read\nTACG\n", &read) == 0);

    /* Start 4 leaves only "AC" in the reference; two deletions fit band 2. */
    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 2, &result) == 0);
    CHECK(result.candidates[0].edit_distance == 2);
    CHECK(result.candidates[0].accepted == 1);
    CHECK(result.mapped == 1 && result.start == 4);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_f.fa");
    remove("test_verify_read_f.fa");
}

static void test_unknown_base_is_not_treated_as_A(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{0, 5}, {0, 0}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_g.fa", ">ref\nACNT\n", &reference) == 0);
    CHECK(load_text("test_verify_read_g.fa", ">read\nACGT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 3, &result) == 0);
    CHECK(result.candidates[0].status == VERIFY_CANDIDATE_UNKNOWN_BASE);
    CHECK(result.candidates[0].edit_distance == VERIFY_DISTANCE_INF);
    CHECK(result.mapped == 0);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_g.fa");
    remove("test_verify_read_g.fa");
}

static void test_invalid_candidate_starts_and_empty_list(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{-2, 5}, {10, 4}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_h.fa", ">ref\nACGT\n", &reference) == 0);
    CHECK(load_text("test_verify_read_h.fa", ">read\nACGT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 2, 2, 3, &result) == 0);
    CHECK(result.candidates[0].status == VERIFY_CANDIDATE_NEGATIVE_START);
    CHECK(result.candidates[1].status == VERIFY_CANDIDATE_START_OUT_OF_RANGE);
    CHECK(result.mapped == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 0, 2, 3, &result) == 0);
    CHECK(result.candidate_count == 0 && result.mapped == 0);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_h.fa");
    remove("test_verify_read_h.fa");
}

static void test_equal_distance_uses_votes_then_start(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{0, 5}, {4, 10}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_i.fa", ">ref\nACGAACGG\n", &reference) == 0);
    CHECK(load_text("test_verify_read_i.fa", ">read\nACGT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 2, 2, 3, &result) == 0);
    CHECK(result.candidates[0].edit_distance == 1);
    CHECK(result.candidates[1].edit_distance == 1);
    CHECK(result.mapped == 1 && result.start == 4 && result.votes == 10);

    /* Equal distance AND equal votes: deterministic final tie goes to start 0. */
    candidates[0].votes = 10;
    CHECK(verify_candidates(&reference, &read, candidates, 2, 2, 3, &result) == 0);
    CHECK(result.mapped == 1 && result.start == 0);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_i.fa");
    remove("test_verify_read_i.fa");
}

static void test_bad_arguments(void)
{
    Seq reference = {0}, read = {0};
    TopCandidate candidates[TOP3_CAPACITY] = {{0, 1}, {0, 0}, {0, 0}};
    VerifyResult result;

    CHECK(load_text("test_verify_ref_j.fa", ">ref\nACGT\n", &reference) == 0);
    CHECK(load_text("test_verify_read_j.fa", ">read\nACGT\n", &read) == 0);

    CHECK(verify_candidates(&reference, &read, candidates, 4, 2, 3, &result) == -1);
    CHECK(verify_candidates(&reference, &read, NULL, 1, 2, 3, &result) == -1);
    CHECK(verify_candidates(&reference, &read, candidates, 1, 2, 3, NULL) == -1);

    seq_free(&reference);
    seq_free(&read);
    remove("test_verify_ref_j.fa");
    remove("test_verify_read_j.fa");
}

int main(void)
{
    test_exact_match_and_best_distance();
    test_substitution_distance();
    test_insertion_or_deletion_cost();
    test_threshold_rejection();
    test_band_length_difference();
    test_reference_right_edge_is_clipped();
    test_unknown_base_is_not_treated_as_A();
    test_invalid_candidate_starts_and_empty_list();
    test_equal_distance_uses_votes_then_start();
    test_bad_arguments();

    if (failures == 0)
        printf("all Stage 7 verification tests passed\n");
    else
        printf("%d Stage 7 test(s) failed\n", failures);

    return failures != 0;
}
