/*
 * End-to-end demo of all 9 stages on planted substitutions (option B).
 *
 * Usage: demo [reference.fasta] [mutations] [coverage] [error_rate] [seed]
 * Defaults: sequence.fasta 20 30 0.01 42
 *
 * Writes answer_key.txt (planted mutations) and calls.txt (what the pipeline
 * found, each marked TP or FP) in the current folder. Positions are 1-based.
 *
 * Build (from the repository root):
 *   gcc -std=c99 -O2 -Wall -Wextra -Isrc/stage1_load -Isrc/stage2_minimizers \
 *     -Isrc/stage3_index -Isrc/stage4_bloom -Isrc/stage5_vote -Isrc/stage6_top3 \
 *     -Isrc/stage7_verify -Isrc/stage8_pileup -Isrc/stage9_score -Isrc/simulate \
 *     src/stage1_load/seq.c src/stage2_minimizers/minimizer.c src/stage3_index/index.c \
 *     src/stage4_bloom/bloom.c src/stage5_vote/vote.c src/stage6_top3/top3.c \
 *     src/stage7_verify/verify.c src/stage8_pileup/pileup.c src/stage9_score/score.c \
 *     src/simulate/simulate.c src/demo/demo.c -lm -o demo
 */
#include "simulate.h"
#include "verify.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static const char LETTER[] = "ACGT";

/* Pipeline settings. k, w, band and threshold match the toy plan. */
#define K 15
#define W 10
#define READ_LEN 150
#define BLOOM_FP 0.01

/* Running totals over all reads, printed at the end. */
typedef struct {
    size_t words, bloom_skipped, lookups, hits;
    size_t placed, unplaced, misplaced;
} DemoTotals;

/* Seconds since start, for the timing lines. */
static double seconds_since(clock_t start)
{
    return (double)(clock() - start) / CLOCKS_PER_SEC;
}

/*
 * Stages 2 to 7 for one read: minimizers, Bloom + index vote, top 3, verify.
 * Returns the accepted start, or -1 if the read is not placed.
 */
static long place_read(const Seq *ref, const Seq *read, const Index *idx,
                       const BloomFilter *bloom, DemoTotals *t)
{
    MinimizerSet mins = {0};
    VoteTree tree = {0};
    VoteStats st = {0, 0, 0, 0};
    Top3Heap top;
    TopCandidate sorted[TOP3_CAPACITY];
    size_t n = 0;
    VerifyResult vr;
    long placed = -1;

    if (minimizers_generate(read, K, W, &mins) == 0 &&
        vote_read(read, &mins, idx, bloom, K, &tree, &st) == 0 &&
        top3_from_vote_tree(&tree, &top) == 0 &&
        top3_get_sorted(&top, sorted, &n) == 0 &&
        verify_candidates(ref, read, sorted, n, VERIFY_DEFAULT_BAND,
                          VERIFY_DEFAULT_MAX_EDIT_DISTANCE, &vr) == 0 &&
        vr.mapped)
        placed = vr.start;

    t->words += st.words;
    t->bloom_skipped += st.bloom_skipped;
    t->lookups += st.lookups;
    t->hits += st.hits;

    vote_tree_clear(&tree);
    minimizers_free(&mins);
    return placed;
}

/* Writes every call with 1-based position and a TP/FP mark from the key. Returns 0 or -1. */
static int save_calls(const char *path, const CallList *calls, const AnswerKey *key)
{
    FILE *f = fopen(path, "w");
    size_t i;

    if (f == NULL)
        return -1;
    fprintf(f, "# position(1-based) ref alt alt_reads depth result\n");
    for (i = 0; i < calls->count; i++) {
        const MutationCall *c = &calls->items[i];
        const KeyEntry *e = key_find(key, c->pos);
        fprintf(f, "%lu %c %c %lu %lu %s\n", (unsigned long)(c->pos + 1),
                LETTER[c->ref], LETTER[c->alt], (unsigned long)c->alt_count,
                (unsigned long)c->depth, (e != NULL && e->alt == c->alt) ? "TP" : "FP");
    }
    return fclose(f) == 0 ? 0 : -1;
}

/* Runs stages 1 to 9 once with the settings from the command line. */
int main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "sequence.fasta";
    size_t mutations = argc > 2 ? (size_t)strtoul(argv[2], NULL, 10) : 20;
    double coverage = argc > 3 ? atof(argv[3]) : 30.0;
    double error_rate = argc > 4 ? atof(argv[4]) : 0.01;
    unsigned long seed = argc > 5 ? strtoul(argv[5], NULL, 10) : 42;

    Seq ref = {0}, sample = {0};
    MinimizerSet ref_mins = {0};
    Index idx = {0};
    BloomFilter bloom = {0};
    AnswerKey key = {0};
    ReadSet reads = {0};
    Pileup pile = {0};
    CallList calls = {0};
    uint32_t *depth = NULL;
    ScoreResult s;
    DemoTotals t = {0, 0, 0, 0, 0, 0, 0};
    SimRng rng;
    clock_t t0, t_all = clock();
    size_t i;
    int rc = 1;

    /* Stage 1: load the reference. */
    if (seq_load_fasta(path, &ref) != 0) {
        fprintf(stderr, "could not load %s\n", path);
        return 1;
    }
    printf("Stage 1  load       %s: %lu letters, %lu N runs\n", path,
           (unsigned long)ref.len, (unsigned long)ref.nruns);

    /* Test data: plant substitutions, build the sample, cut reads. */
    sim_rng_seed(&rng, seed);
    if (sim_plant(&ref, mutations, &rng, &key) != 0 ||
        sim_make_sample(&ref, &key, &sample) != 0 ||
        sim_reads(&sample, READ_LEN, coverage, error_rate, &rng, &reads) != 0 ||
        key_save("answer_key.txt", &key) != 0) {
        fprintf(stderr, "could not generate test data\n");
        goto done;
    }
    printf("Generate            %lu planted substitutions, %lu reads of %d letters "
           "(%.0fx coverage, %.2f%% errors, seed %lu) -> answer_key.txt\n",
           (unsigned long)key.count, (unsigned long)reads.count, READ_LEN, coverage,
           100.0 * error_rate, seed);

    /* Stages 2 to 4 on the reference: minimizers, index, Bloom filter. */
    t0 = clock();
    if (minimizers_generate(&ref, K, W, &ref_mins) != 0 ||
        index_build(&ref, &ref_mins, K, 0, &idx) != 0 ||
        bloom_build_from_minimizers(&ref, &ref_mins, K, BLOOM_FP, &bloom) != 0) {
        fprintf(stderr, "could not build the reference index\n");
        goto done;
    }
    printf("Stage 2  minimizers %lu kept words (k=%d, w=%d)\n",
           (unsigned long)ref_mins.count, K, W);
    printf("Stage 3  index      %lu distinct words in %lu buckets\n",
           (unsigned long)idx.nentries, (unsigned long)idx.nbuckets);
    printf("Stage 4  bloom      built (target false-positive rate %.0f%%)   [%.3f s]\n",
           100.0 * BLOOM_FP, seconds_since(t0));

    /* Stages 5 to 8 for every read: place it, then add it to the pile-up. */
    t0 = clock();
    if (pileup_init(&pile, ref.len) != 0)
        goto done;
    for (i = 0; i < reads.count; i++) {
        long start = place_read(&ref, &reads.reads[i], &idx, &bloom, &t);
        if (start < 0) {
            t.unplaced++;
            continue;
        }
        t.placed++;
        if ((size_t)start != reads.true_start[i])
            t.misplaced++;
        if (pileup_add_read(&pile, &reads.reads[i], (size_t)start) != 0)
            goto done;
    }
    printf("Stage 5  vote       %lu read words, %lu skipped by Bloom, %lu lookups, %lu votes\n",
           (unsigned long)t.words, (unsigned long)t.bloom_skipped,
           (unsigned long)t.lookups, (unsigned long)t.hits);
    printf("Stage 6-7 top3+verify %lu placed, %lu not placed, %lu placed at the wrong start   [%.3f s]\n",
           (unsigned long)t.placed, (unsigned long)t.unplaced, (unsigned long)t.misplaced,
           seconds_since(t0));

    /* Stage 8: call mutations from the pile-up. */
    depth = (uint32_t *)malloc(ref.len * sizeof *depth);
    if (depth == NULL || pileup_depths(&pile, depth) != 0 ||
        pileup_call(&pile, &ref, PILEUP_DEFAULT_MIN_DEPTH,
                    PILEUP_DEFAULT_MIN_ALT_PERCENT, &calls) != 0)
        goto done;
    printf("Stage 8  pile-up    %lu calls (min depth %lu, alt >= %lu%%) -> calls.txt\n",
           (unsigned long)calls.count, (unsigned long)PILEUP_DEFAULT_MIN_DEPTH,
           (unsigned long)PILEUP_DEFAULT_MIN_ALT_PERCENT);

    /* Stage 9: score against the answer key. */
    if (score_calls(&calls, &key, depth, ref.len, PILEUP_DEFAULT_MIN_DEPTH, &s) != 0 ||
        save_calls("calls.txt", &calls, &key) != 0)
        goto done;
    printf("Stage 9  score      TP %lu  FP %lu  FN %lu  (wrong letter %lu)\n",
           (unsigned long)s.tp, (unsigned long)s.fp, (unsigned long)s.fn,
           (unsigned long)s.wrong_alt);
    printf("                    precision %.1f%%  recall %.1f%% of all planted, "
           "%.1f%% of %lu callable\n",
           s.precision_pct, s.recall_all_pct, s.recall_callable_pct,
           (unsigned long)s.callable);

    /* List any planted mutation that was missed, with its depth, to explain why. */
    for (i = 0; i < key.count; i++) {
        size_t c;
        int found = 0;
        for (c = 0; c < calls.count; c++)
            if (calls.items[c].pos == key.items[i].pos && calls.items[c].alt == key.items[i].alt)
                found = 1;
        if (!found)
            printf("  missed: position %lu %c -> %c, depth %lu\n",
                   (unsigned long)(key.items[i].pos + 1), LETTER[key.items[i].ref],
                   LETTER[key.items[i].alt], (unsigned long)depth[key.items[i].pos]);
    }
    printf("Total time %.3f s\n", seconds_since(t_all));
    rc = 0;

done:
    free(depth);
    calls_free(&calls);
    pileup_free(&pile);
    sim_reads_free(&reads);
    key_free(&key);
    bloom_free(&bloom);
    index_free(&idx);
    minimizers_free(&ref_mins);
    seq_free(&sample);
    seq_free(&ref);
    return rc;
}
