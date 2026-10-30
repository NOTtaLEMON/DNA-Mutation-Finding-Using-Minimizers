#include "index.h"

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

/* Writes text to a new file at path and loads it with stage 1. Returns 0 on success. */
static int load_text(const char *path, const char *text, Seq *s)
{
    FILE *f = fopen(path, "wb");
    if (f == NULL)
        return -1;
    fwrite(text, 1, strlen(text), f);
    fclose(f);
    return seq_load_fasta(path, s);
}

/* Returns 1 if pos appears in the entry's position list, 0 otherwise. */
static int entry_has_pos(const IndexEntry *e, size_t pos)
{
    const PosNode *p;
    for (p = e->pos_head; p != NULL; p = p->next)
        if (p->pos == pos)
            return 1;
    return 0;
}

/* Word IDs use the oracle's packing: first letter in the highest bits. */
static void test_word_id(void)
{
    Seq s;
    CHECK(load_text("test_idx_a.fasta", ">a\nAAGCACGT\n", &s) == 0);
    CHECK(word_id_at(&s, 0, 4) == 9);     /* AAGC = 00 00 10 01 */
    CHECK(word_id_at(&s, 4, 4) == 27);    /* ACGT = 00 01 10 11 */
    seq_free(&s);
    remove("test_idx_a.fasta");
}

/*
 * Hand-made kept words on ACGTACGTTTGCA (k = 4), with only 2 buckets so chains are shared:
 *   ACGT at 0 and 4 (same word twice), TTGC at 8 with the SAME hash as ACGT (a hash collision),
 *   CGTA at 1 with its own hash.
 */
static void test_repeats_and_collisions(void)
{
    Seq s;
    Index idx;
    Minimizer items[4];
    MinimizerSet mins;
    const IndexEntry *e;
    size_t probes;
    uint64_t id_acgt, id_ttgc, id_cgta;

    CHECK(load_text("test_idx_b.fasta", ">b\nACGTACGTTTGCA\n", &s) == 0);
    id_acgt = word_id_at(&s, 0, 4);
    id_ttgc = word_id_at(&s, 8, 4);
    id_cgta = word_id_at(&s, 1, 4);

    items[0].hash = 111; items[0].pos = 0;   /* ACGT */
    items[1].hash = 222; items[1].pos = 1;   /* CGTA */
    items[2].hash = 111; items[2].pos = 4;   /* ACGT again */
    items[3].hash = 111; items[3].pos = 8;   /* TTGC: same hash, different word */
    mins.items = items; mins.count = 4; mins.capacity = 4;

    CHECK(index_build(&s, &mins, 4, 2, &idx) == 0);
    CHECK(idx.nentries == 3);

    e = index_lookup(&idx, 111, id_acgt, &probes);
    CHECK(e != NULL);
    if (e != NULL) {
        CHECK(e->npos == 2);
        CHECK(e->pos_head->pos == 0 && e->pos_tail->pos == 4);   /* positions kept in order */
    }

    e = index_lookup(&idx, 111, id_ttgc, &probes);
    CHECK(e != NULL && e->npos == 1 && e->pos_head->pos == 8);   /* not merged with ACGT */
    CHECK(probes == 2);                                            /* ACGT checked first, then TTGC */

    e = index_lookup(&idx, 222, id_cgta, NULL);
    CHECK(e != NULL && e->npos == 1 && e->pos_head->pos == 1);

    CHECK(index_lookup(&idx, 111, 0, NULL) == NULL);         /* AAAA: never kept */
    CHECK(index_lookup(&idx, 999, id_acgt, NULL) == NULL);   /* right word, wrong hash */

    index_free(&idx);
    seq_free(&s);
    remove("test_idx_b.fasta");
}

/* Bad arguments are rejected and leave the output empty. */
static void test_bad_arguments(void)
{
    Seq s;
    Index idx;
    Minimizer item;
    MinimizerSet mins;

    CHECK(load_text("test_idx_c.fasta", ">c\nACGTACGT\n", &s) == 0);
    item.hash = 1; item.pos = 5;             /* 5 + 4 > 8: runs past the end */
    mins.items = &item; mins.count = 1; mins.capacity = 1;

    CHECK(index_build(&s, &mins, 33, 0, &idx) == -1);   /* k too large for a 64-bit word ID */
    CHECK(index_build(&s, &mins, 4, 0, &idx) == -1);    /* position out of range */
    CHECK(idx.buckets == NULL && idx.nentries == 0);

    seq_free(&s);
    remove("test_idx_c.fasta");
}

/*
 * Full stages 1 to 3 on the real reference: every kept word must be found again by lookup,
 * with its position in the entry's list. Prints table statistics.
 */
static void test_real_reference(const char *path)
{
    const size_t k = 15, w = 10;
    Seq s;
    MinimizerSet mins = {0};
    Index idx;
    size_t i, b, longest = 0, used = 0, total_probes = 0;

    if (seq_load_fasta(path, &s) != 0) {
        printf("(could not open %s, skipping)\n", path);
        return;
    }
    CHECK(minimizers_generate(&s, k, w, &mins) == 0);
    CHECK(index_build(&s, &mins, k, 0, &idx) == 0);

    for (i = 0; i < mins.count; i++) {
        size_t probes;
        uint64_t id = word_id_at(&s, mins.items[i].pos, k);
        const IndexEntry *e = index_lookup(&idx, mins.items[i].hash, id, &probes);
        CHECK(e != NULL && entry_has_pos(e, mins.items[i].pos));
        total_probes += probes;
    }

    for (b = 0; b < idx.nbuckets; b++) {
        size_t len = 0;
        const IndexEntry *e;
        for (e = idx.buckets[b]; e != NULL; e = e->next)
            len++;
        if (len > 0)
            used++;
        if (len > longest)
            longest = len;
    }

    printf("%s, k=%lu, w=%lu: %lu kept words, %lu distinct, %lu buckets (%lu used), "
           "longest chain %lu, average probes per lookup %.2f\n",
           path, (unsigned long)k, (unsigned long)w, (unsigned long)mins.count,
           (unsigned long)idx.nentries, (unsigned long)idx.nbuckets, (unsigned long)used,
           (unsigned long)longest, mins.count ? (double)total_probes / (double)mins.count : 0.0);

    index_free(&idx);
    minimizers_free(&mins);
    seq_free(&s);
}

/* Runs all checks, then the full pipeline on argv[1] (default sequence.fasta). */
int main(int argc, char **argv)
{
    test_word_id();
    test_repeats_and_collisions();
    test_bad_arguments();
    test_real_reference(argc > 1 ? argv[1] : "sequence.fasta");

    if (failures == 0)
        printf("all checks passed\n");
    else
        printf("%d check(s) failed\n", failures);
    return failures != 0;
}
