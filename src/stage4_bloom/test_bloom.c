#include "bloom.h"
#include "index.h"

#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond)                                                            \
    do {                                                                       \
        if (!(cond)) {                                                         \
            printf("FAIL line %d: %s\n", __LINE__, #cond);                   \
            failures++;                                                        \
        }                                                                      \
    } while (0)

/* Write a tiny FASTA file and load it through the real Stage 1 loader. */
static int load_text(const char *path, const char *text, Seq *seq)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL)
        return -1;
    if (fwrite(text, 1, strlen(text), f) != strlen(text)) {
        fclose(f);
        return -1;
    }
    fclose(f);
    return seq_load_fasta(path, seq);
}

/*
 * Test that every reference minimizer inserted into the Bloom filter is
 * reported as POSSIBLY PRESENT (there must be no false negatives).
 */
static void test_inserted_keys_are_present(void)
{
    Seq reference;
    BloomFilter filter = {0};
    Minimizer items[4];
    MinimizerSet mins;
    uint64_t id_acgt;
    uint64_t id_cgta;

    CHECK(load_text("test_bloom_ref.fasta", ">ref\nACGTACGT\n", &reference) == 0);

    /*
     * These are the reference minimizer occurrences for this test. The same
     * DNA word at positions 0 and 4 has the same Stage 2 hash and word ID.
     * Hash values are test fixtures, not outputs of a particular hash call.
     */
    items[0].hash = 111; items[0].pos = 0; /* ACGT */
    items[1].hash = 222; items[1].pos = 1; /* CGTA */
    items[2].hash = 111; items[2].pos = 4; /* ACGT repeated */
    items[3].hash = 333; items[3].pos = 2; /* GTAC */
    mins.items = items;
    mins.count = 4;
    mins.capacity = 4;

    CHECK(bloom_build_from_minimizers(&reference, &mins, 4, 0.01, &filter) == 0);
    CHECK(filter.bits != NULL);
    CHECK(filter.expected_items == 4);
    CHECK(filter.inserted_items == 4);
    CHECK(filter.bit_count > 0 && filter.hash_count > 0);

    id_acgt = word_id_at(&reference, 0, 4);
    id_cgta = word_id_at(&reference, 1, 4);

    CHECK(bloom_maybe_contains(&filter, 111, id_acgt) == 1);
    CHECK(bloom_maybe_contains(&filter, 222, id_cgta) == 1);
    CHECK(bloom_maybe_contains(&filter, 333, word_id_at(&reference, 2, 4)) == 1);

    bloom_free(&filter);
    seq_free(&reference);
    remove("test_bloom_ref.fasta");
}

/*
 * A Bloom filter may return a false positive, so do not assume one arbitrary
 * absent key must be rejected. Search a range and ensure the filter can give
 * the DEFINITELY ABSENT answer for keys that were not inserted.
 */
static void test_definitely_absent_result(void)
{
    Seq reference;
    BloomFilter filter = {0};
    Minimizer item;
    MinimizerSet mins;
    uint64_t inserted_id;
    uint64_t candidate;
    int found_absent = 0;

    CHECK(load_text("test_bloom_absent.fasta", ">ref\nACGTACGT\n", &reference) == 0);

    item.hash = 111;
    item.pos = 0;
    mins.items = &item;
    mins.count = 1;
    mins.capacity = 1;

    CHECK(bloom_build_from_minimizers(&reference, &mins, 4, 0.001, &filter) == 0);
    inserted_id = word_id_at(&reference, 0, 4);

    /* First verify the inserted key is never reported absent. */
    CHECK(bloom_maybe_contains(&filter, 111, inserted_id) == 1);

    /* Try different non-inserted pairs until one is certainly rejected. */
    for (candidate = 1; candidate < 10000; candidate++) {
        int answer = bloom_maybe_contains(&filter,
                                          UINT64_C(100000) + candidate,
                                          UINT64_C(900000) + candidate);
        if (answer == 0) {
            found_absent = 1;
            break;
        }
    }
    CHECK(found_absent == 1);

    bloom_free(&filter);
    seq_free(&reference);
    remove("test_bloom_absent.fasta");
}

/* An empty minimizer set gives an empty filter: every query is definitely absent. */
static void test_empty_reference_minimizers(void)
{
    Seq reference;
    BloomFilter filter = {0};
    MinimizerSet mins = {0};

    CHECK(load_text("test_bloom_empty.fasta", ">ref\nACGT\n", &reference) == 0);
    CHECK(bloom_build_from_minimizers(&reference, &mins, 4, 0.01, &filter) == 0);
    CHECK(filter.bits == NULL);
    CHECK(bloom_maybe_contains(&filter, 123, 456) == 0);

    bloom_free(&filter);
    seq_free(&reference);
    remove("test_bloom_empty.fasta");
}

/* Invalid parameters/positions must fail and leave the output filter empty. */
static void test_invalid_inputs(void)
{
    Seq reference;
    BloomFilter filter = {0};
    Minimizer item;
    MinimizerSet mins;

    CHECK(load_text("test_bloom_bad.fasta", ">ref\nACGTACGT\n", &reference) == 0);
    item.hash = 1;
    item.pos = 6; /* k=4 would run beyond the 8-base reference. */
    mins.items = &item;
    mins.count = 1;
    mins.capacity = 1;

    CHECK(bloom_build_from_minimizers(&reference, &mins, 4, 0.01, &filter) == -1);
    CHECK(filter.bits == NULL && filter.bit_count == 0);

    item.pos = 0;
    CHECK(bloom_build_from_minimizers(&reference, &mins, 0, 0.01, &filter) == -1);
    CHECK(bloom_build_from_minimizers(&reference, &mins, 4, 0.0, &filter) == -1);
    CHECK(bloom_build_from_minimizers(&reference, &mins, 4, 1.0, &filter) == -1);

    CHECK(bloom_maybe_contains(NULL, 1, 2) == -1);

    bloom_free(&filter);
    seq_free(&reference);
    remove("test_bloom_bad.fasta");
}

int main(void)
{
    test_inserted_keys_are_present();
    test_definitely_absent_result();
    test_empty_reference_minimizers();
    test_invalid_inputs();

    if (failures == 0)
        printf("all Bloom filter checks passed\n");
    else
        printf("%d Bloom filter check(s) failed\n", failures);

    return failures != 0;
}
