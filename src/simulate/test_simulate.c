#include "simulate.h"

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

/* Writes letters as a one-record FASTA and loads it through stage 1. */
static int load_letters(const char *letters, Seq *seq)
{
    FILE *f = fopen("test_simulate_tmp.fasta", "wb");
    int rc;

    if (f == NULL)
        return -1;
    fprintf(f, ">s\n%s\n", letters);
    fclose(f);
    rc = seq_load_fasta("test_simulate_tmp.fasta", seq);
    remove("test_simulate_tmp.fasta");
    return rc;
}

/* Same seed, same numbers; different seed, different numbers; ranges hold. */
static void test_rng(void)
{
    SimRng a, b, c;
    int i, same = 1, differ = 0;

    sim_rng_seed(&a, 42);
    sim_rng_seed(&b, 42);
    sim_rng_seed(&c, 43);
    for (i = 0; i < 1000; i++) {
        uint64_t x = sim_rng_next(&a), y = sim_rng_next(&b), z = sim_rng_next(&c);
        same &= x == y;
        differ |= x != z;
    }
    CHECK(same && differ);
    for (i = 0; i < 1000; i++) {
        double u = sim_rng_unit(&a);
        CHECK(sim_rng_below(&a, 7) < 7);
        CHECK(u >= 0.0 && u < 1.0);
    }
}

/* Planted positions are distinct, sorted, never in an N run, and the sample carries them. */
static void test_plant_and_sample(void)
{
    Seq ref, sample;
    AnswerKey key = {0};
    SimRng rng;
    size_t i, changed = 0;

    CHECK(load_letters("ACGTACGTNNNNACGTACGTACGT", &ref) == 0);   /* 24 letters, N at 8..11 */
    sim_rng_seed(&rng, 1);
    CHECK(sim_plant(&ref, 20, &rng, &key) == 0);                 /* all 20 usable positions */
    CHECK(key.count == 20);
    for (i = 0; i < key.count; i++) {
        CHECK(key.items[i].pos < 8 || key.items[i].pos >= 12);
        CHECK(key.items[i].ref == seq_get(&ref, key.items[i].pos));
        CHECK(key.items[i].alt != key.items[i].ref);
        if (i > 0)
            CHECK(key.items[i - 1].pos < key.items[i].pos);
    }
    CHECK(sim_plant(&ref, 21, &rng, &key) == -1);                /* more than usable */

    CHECK(sim_plant(&ref, 5, &rng, &key) == 0);
    CHECK(sim_make_sample(&ref, &key, &sample) == 0);
    CHECK(sample.len == ref.len && sample.nruns == 1);
    for (i = 0; i < ref.len; i++) {
        const KeyEntry *e = key_find(&key, i);
        if (e != NULL)
            CHECK(seq_get(&sample, i) == e->alt);
        else
            CHECK(seq_get(&sample, i) == seq_get(&ref, i));
        changed += seq_get(&sample, i) != seq_get(&ref, i);
    }
    CHECK(changed == 5);

    seq_free(&sample);
    key_free(&key);
    seq_free(&ref);
}

/* Error-free reads are exact copies at their true start and never touch the N run. */
static void test_reads(void)
{
    Seq ref;
    ReadSet set = {0};
    SimRng rng;
    size_t r, i, diffs = 0, letters = 0;

    CHECK(load_letters("ACGTACGTNNNNACGTACGTACGTTTGCA", &ref) == 0);   /* 29 letters */
    sim_rng_seed(&rng, 5);
    CHECK(sim_reads(&ref, 6, 10.0, 0.0, &rng, &set) == 0);
    CHECK(set.count == 48);                                      /* 10 * 29 / 6, rounded */
    for (r = 0; r < set.count; r++) {
        size_t s = set.true_start[r];
        CHECK(set.reads[r].len == 6 && set.reads[r].nruns == 0);
        CHECK(s + 6 <= 8 || s >= 12);
        for (i = 0; i < 6; i++)
            CHECK(seq_get(&set.reads[r], i) == seq_get(&ref, s + i));
    }
    sim_reads_free(&set);

    /* With a 10% error rate about one letter in ten differs. */
    CHECK(sim_reads(&ref, 6, 2000.0, 0.10, &rng, &set) == 0);
    for (r = 0; r < set.count; r++)
        for (i = 0; i < 6; i++, letters++)
            diffs += seq_get(&set.reads[r], i) != seq_get(&ref, set.true_start[r] + i);
    CHECK(diffs * 100 > letters * 8 && diffs * 100 < letters * 12);
    sim_reads_free(&set);
    CHECK(sim_reads(&ref, 0, 1.0, 0.0, &rng, &set) == -1);

    seq_free(&ref);
}

/* Runs all checks. */
int main(void)
{
    test_rng();
    test_plant_and_sample();
    test_reads();

    if (failures == 0)
        printf("all checks passed\n");
    else
        printf("%d check(s) failed\n", failures);
    return failures != 0;
}
