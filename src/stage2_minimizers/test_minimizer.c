#include "minimizer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(condition)                                                       \
    do {                                                                        \
        if (!(condition)) {                                                    \
            printf("FAIL line %d: %s\n", __LINE__, #condition);              \
            failures++;                                                        \
        }                                                                       \
    } while (0)

/*
 * Load a tiny FASTA string through Stage 1 so these tests exercise the real
 * packed-Seq representation and its N-run metadata.
 */
static int load_text(const char *path, const char *text, Seq *seq)
{
    FILE *file = fopen(path, "wb");
    size_t length = strlen(text);

    if (file == NULL)
        return -1;
    if (fwrite(text, 1, length, file) != length) {
        fclose(file);
        return -1;
    }
    if (fclose(file) != 0)
        return -1;

    return seq_load_fasta(path, seq);
}

static void test_repeated_kmers_emit_one_occurrence(void)
{
    /* All 3-mers are AAA and have the same hash. Leftmost tie-breaking means
       each 2-k-mer window selects its leftmost occurrence. As the window
       slides, that occurrence moves right, so all seven positions are kept. */
    Seq seq = {0};
    MinimizerSet mins = {0};

    CHECK(load_text("test_minimizer_repeat.fa", ">repeat\nAAAAAAAAAA\n", &seq) == 0);
    CHECK(minimizers_generate(&seq, 3, 2, &mins) == 0);
    CHECK(mins.count == 7);
    if (mins.count == 7) {
        size_t i;
        for (i = 0; i < mins.count; i++)
            CHECK(mins.items[i].pos == i);
    }

    minimizers_free(&mins);
    seq_free(&seq);
    remove("test_minimizer_repeat.fa");
}

static void test_unknown_bases_break_windows(void)
{
    /*
     * ACGTNNACGT, k=3 gives valid k-mer runs at starts [0,2) and [6,8).
     * The Ns invalidate all k-mers that overlap positions 4 or 5.
     * With window_kmers=2, each valid run contributes one minimizer.
     */
    Seq seq = {0};
    MinimizerSet mins = {0};
    size_t i;

    CHECK(load_text("test_minimizer_n.fa", ">n-test\nACGTNNACGT\n", &seq) == 0);
    CHECK(minimizers_generate(&seq, 3, 2, &mins) == 0);
    CHECK(mins.count == 2);

    for (i = 0; i < mins.count; i++) {
        /* Each output is a valid 3-mer and must not overlap the N-run [4,6). */
        CHECK(mins.items[i].pos + 3 <= 4 || mins.items[i].pos >= 6);
    }

    minimizers_free(&mins);
    seq_free(&seq);
    remove("test_minimizer_n.fa");
}

static void test_short_valid_segment(void)
{
    /* Only two 3-mers exist, fewer than the requested window size of five.
       Our documented policy emits the minimum of this short valid segment. */
    Seq seq = {0};
    MinimizerSet mins = {0};

    CHECK(load_text("test_minimizer_short.fa", ">short\nACGT\n", &seq) == 0);
    CHECK(minimizers_generate(&seq, 3, 5, &mins) == 0);
    CHECK(mins.count == 1);
    if (mins.count == 1)
        CHECK(mins.items[0].pos < 2);

    minimizers_free(&mins);
    seq_free(&seq);
    remove("test_minimizer_short.fa");
}

static void test_bad_parameters_and_empty_result(void)
{
    Seq seq = {0};
    MinimizerSet mins = {0};

    CHECK(load_text("test_minimizer_args.fa", ">args\nACGTAC\n", &seq) == 0);
    CHECK(minimizers_generate(&seq, 0, 2, &mins) == -1);
    CHECK(minimizers_generate(&seq, 3, 0, &mins) == -1);
    CHECK(minimizers_generate(&seq, 20, 2, &mins) == 0);
    CHECK(mins.count == 0);

    minimizers_free(&mins);
    seq_free(&seq);
    remove("test_minimizer_args.fa");
}

static void test_known_window_minima(void)
{
    /*
     * Expected positions were calculated from 64-bit FNV-1a hashes of each
     * 4-mer, using windows of three consecutive k-mers and leftmost tie-breaks.
     * This checks the actual sliding-window selection, not just output shape.
     */
    static const size_t expected_positions[] = {0, 2, 5, 7, 9, 12, 15, 18};
    Seq seq = {0};
    MinimizerSet mins = {0};
    size_t i;

    CHECK(load_text("test_minimizer_order.fa", ">order\nACGTTGCATGTCAGTACGATCG\n", &seq) == 0);
    CHECK(minimizers_generate(&seq, 4, 3, &mins) == 0);
    CHECK(mins.count == sizeof expected_positions / sizeof expected_positions[0]);

    if (mins.count == sizeof expected_positions / sizeof expected_positions[0]) {
        for (i = 0; i < mins.count; i++)
            CHECK(mins.items[i].pos == expected_positions[i]);
    }

    minimizers_free(&mins);
    seq_free(&seq);
    remove("test_minimizer_order.fa");
}

int main(void)
{
    test_repeated_kmers_emit_one_occurrence();
    test_unknown_bases_break_windows();
    test_short_valid_segment();
    test_bad_parameters_and_empty_result();
    test_known_window_minima();

    if (failures == 0)
        printf("All Step 2 minimizer tests passed.\n");
    else
        printf("%d Step 2 minimizer test(s) failed.\n", failures);

    return failures != 0;
}
