/*
 * Trace of one read through all 9 stages, printed as JSON for the web
 * visualizer (web/index.html, served by web/server.py).
 *
 * Usage: trace [reference.fasta] [read_index] [mutations] [coverage]
 *              [error_rate] [seed] [min_depth] [max_edit] [k] [w]
 * read_index -1 (the default) picks the first placed read that covers a
 * planted mutation. Defaults: 20 mutations, 30x, 1% errors, seed 42,
 * min depth 3, edit limit 3, k 15, w 10, so with no settings the trace
 * matches demo.exe.
 *
 * The pipeline runs over every read exactly as in demo.c. For the chosen read
 * the stage 4/5 loop and the stage 6 offers are repeated step by step so each
 * decision can be shown. Two things are recomputed for display only: the
 * stage 7 DP table and the stage 4 Bloom bit positions (the latter is checked
 * against bloom_maybe_contains for every word, see bloom_positions).
 *
 * Build (from the repository root): same command as demo.c, with
 * src/demo/trace.c in place of src/demo/demo.c and -o trace.
 */
#include "simulate.h"
#include "verify.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char LETTER[] = "ACGT";

#define READ_LEN 150
#define BLOOM_FP 0.01
#define DP_SHOW_INF ((size_t)-1)
#define MAX_BLOOM_HASHES 64

/* k-mer length and window size; set from the command line in main. */
static size_t K = 15, W = 10;

/* Settings that the visualizer can change; everything else is fixed above. */
typedef struct {
    size_t mutations;
    double coverage;
    double error_rate;
    unsigned long seed;
    uint32_t min_depth;
    size_t max_edit;
} TraceSettings;

/* Whole-run totals for stages 4 and 5, shown in the page's pipeline funnel. */
typedef struct {
    size_t words, bloom_skipped, lookups, hits;
} RunTotals;

/* Prints a 64-bit hash as 16 hex digits in two halves (old MinGW printf has no %llx). */
static void print_hash(uint64_t h)
{
    printf("\"%08lx%08lx\"", (unsigned long)(h >> 32), (unsigned long)(h & 0xFFFFFFFFu));
}

/* Prints letters [from, from+len) of s as a JSON string. */
static void print_letters(const Seq *s, size_t from, size_t len)
{
    size_t i;

    putchar('"');
    for (i = 0; i < len; i++)
        putchar(LETTER[seq_get(s, from + i)]);
    putchar('"');
}

/*
 * Stages 2 to 7 for one read, as in demo.c, with edit limit max_edit.
 * Adds to totals and writes how many candidates stage 6 short-listed.
 * Returns the accepted start, or -1.
 */
static long place_read(const Seq *ref, const Seq *read, const Index *idx,
                       const BloomFilter *bloom, size_t max_edit, RunTotals *totals,
                       size_t *candidates)
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
        verify_candidates(ref, read, sorted, n, VERIFY_DEFAULT_BAND, max_edit, &vr) == 0 &&
        vr.mapped)
        placed = vr.start;

    *candidates = n;
    totals->words += st.words;
    totals->bloom_skipped += st.bloom_skipped;
    totals->lookups += st.lookups;
    totals->hits += st.hits;
    vote_tree_clear(&tree);
    minimizers_free(&mins);
    return placed;
}

/* Prints a vote tree node and its subtrees as nested JSON (null for no node). */
static void print_tree(const VoteNode *node)
{
    if (node == NULL) {
        printf("null");
        return;
    }
    printf("{\"start\":%ld,\"count\":%lu,\"left\":", node->start, (unsigned long)node->count);
    print_tree(node->left);
    printf(",\"right\":");
    print_tree(node->right);
    putchar('}');
}

/* State shared with the stage 6 walk callback. */
typedef struct {
    Top3Heap *top;
    int first;
} OfferTrace;

/* Prints the heap array as [[start, votes], ...] in array order (index 0 = root). */
static void print_heap(const Top3Heap *top)
{
    size_t i;

    putchar('[');
    for (i = 0; i < top->count; i++)
        printf("%s[%ld,%lu]", i ? "," : "", top->heap[i].start, (unsigned long)top->heap[i].votes);
    putchar(']');
}

/* Walk callback: offers one node to the heap and prints what happened plus the heap after. */
static void trace_offer(long start, size_t count, void *ctx)
{
    OfferTrace *t = (OfferTrace *)ctx;
    size_t before = t->top->count;
    TopCandidate root_before = t->top->heap[0];
    const char *action;

    top3_offer(t->top, start, count);
    if (before < TOP3_CAPACITY)
        action = "insert";
    else if (t->top->heap[0].start != root_before.start || t->top->heap[0].votes != root_before.votes)
        action = "replace";
    else
        action = "discard";

    printf("%s{\"start\":%ld,\"votes\":%lu,\"action\":\"%s\",\"root_before\":",
           t->first ? "" : ",", start, (unsigned long)count, action);
    if (before > 0)
        printf("[%ld,%lu]", root_before.start, (unsigned long)root_before.votes);
    else
        printf("null");
    printf(",\"heap\":");
    print_heap(t->top);
    putchar('}');
    t->first = 0;
}

/* Returns the smaller of three DP values. */
static size_t min3(size_t a, size_t b, size_t c)
{
    size_t m = a < b ? a : b;
    return m < c ? m : c;
}

/*
 * Recomputes stage 7's banded table for display and prints only the band:
 * one row per read prefix i, 2*band+1 cells for j = i-band .. i+band
 * (null where j is outside the window or the cell is unreachable).
 */
static void print_dp_band(const Seq *read, const Seq *ref, size_t start,
                          size_t window, size_t band)
{
    size_t rows = read->len + 1, width = 2 * band + 1;
    size_t *dp = (size_t *)malloc(rows * width * sizeof *dp);
    size_t i, d;

    if (dp == NULL) {
        printf("null");
        return;
    }
    for (i = 0; i < rows * width; i++)
        dp[i] = DP_SHOW_INF;

    /* Cell (i, j) lives at dp[i * width + (j - i + band)]. */
    for (i = 0; i < rows; i++) {
        for (d = 0; d < width; d++) {
            long jl = (long)i - (long)band + (long)d;
            size_t j, best = DP_SHOW_INF;
            if (jl < 0 || (size_t)jl > window)
                continue;
            j = (size_t)jl;
            if (i == 0) {
                best = j;
            } else if (j == 0) {
                best = i;
            } else {
                size_t diag = dp[(i - 1) * width + d];                         /* (i-1, j-1) */
                size_t up = d + 1 < width ? dp[(i - 1) * width + d + 1] : DP_SHOW_INF;   /* (i-1, j) */
                size_t left = d > 0 ? dp[i * width + d - 1] : DP_SHOW_INF;     /* (i, j-1) */
                size_t cost = seq_get(read, i - 1) != seq_get(ref, start + j - 1);
                best = min3(diag == DP_SHOW_INF ? DP_SHOW_INF : diag + cost,
                            up == DP_SHOW_INF ? DP_SHOW_INF : up + 1,
                            left == DP_SHOW_INF ? DP_SHOW_INF : left + 1);
            }
            dp[i * width + d] = best;
        }
    }

    putchar('[');
    for (i = 0; i < rows; i++) {
        printf("%s[", i ? "," : "");
        for (d = 0; d < width; d++) {
            size_t v = dp[i * width + d];
            if (v == DP_SHOW_INF)
                printf("%snull", d ? "," : "");
            else
                printf("%s%lu", d ? "," : "", (unsigned long)v);
        }
        putchar(']');
    }
    putchar(']');
    free(dp);
}

/* Name of a stage 7 candidate status, for display. */
static const char *status_name(VerifyCandidateStatus s)
{
    switch (s) {
    case VERIFY_CANDIDATE_EVALUATED: return "evaluated";
    case VERIFY_CANDIDATE_NEGATIVE_START: return "negative start";
    case VERIFY_CANDIDATE_START_OUT_OF_RANGE: return "start out of range";
    case VERIFY_CANDIDATE_UNKNOWN_BASE: return "touches N";
    case VERIFY_CANDIDATE_OUTSIDE_BAND: return "outside band";
    }
    return "unknown";
}

/* Returns 1 if a read with true start s covers any planted position, else 0. */
static int covers_mutation(const AnswerKey *key, size_t s)
{
    size_t i;

    for (i = 0; i < key->count; i++)
        if (key->items[i].pos >= s && key->items[i].pos < s + READ_LEN)
            return 1;
    return 0;
}

/* Display copy of bloom.c's SplitMix64 mixer (that one is static). */
static uint64_t mix64(uint64_t x)
{
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

/* Display copy of bloom.c's rotate-left. */
static uint64_t rotl64(uint64_t x, unsigned int r)
{
    return (x << r) | (x >> (64U - r));
}

/*
 * Writes the bit positions bloom.c tests for (hash, word_id), using the same
 * double hashing h1 + i*h2 mod m. Returns 1 if they agree with
 * bloom_maybe_contains (all bits set <=> "maybe"), 0 if this display copy has
 * drifted from bloom.c.
 */
static int bloom_positions(const BloomFilter *bloom, uint64_t hash, uint64_t word_id,
                           size_t *out)
{
    uint64_t h1 = mix64(hash ^ rotl64(word_id, 23) ^ UINT64_C(0x243f6a8885a308d3));
    uint64_t h2 = mix64(word_id ^ rotl64(hash, 41) ^ UINT64_C(0x13198a2e03707344)) | UINT64_C(1);
    size_t i;
    int all_set = 1;

    for (i = 0; i < bloom->hash_count && i < MAX_BLOOM_HASHES; i++) {
        out[i] = (size_t)((h1 + (uint64_t)i * h2) % (uint64_t)bloom->bit_count);
        if (!(bloom->bits[out[i] / 8] & (1u << (out[i] % 8))))
            all_set = 0;
    }
    return all_set == (bloom_maybe_contains(bloom, hash, word_id) == 1);
}

/* Runs the whole pipeline, then prints the JSON trace for one read. */
int main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "sequence.fasta";
    long want = argc > 2 ? strtol(argv[2], NULL, 10) : -1;
    TraceSettings set;
    RunTotals totals = {0, 0, 0, 0};
    Seq ref = {0}, sample = {0};
    MinimizerSet ref_mins = {0};
    Index idx = {0};
    BloomFilter bloom = {0};
    AnswerKey key = {0};
    ReadSet reads = {0};
    Pileup pile = {0};
    CallList calls = {0};
    uint32_t *depth = NULL;
    long *placed_at = NULL;
    size_t *ncand = NULL;
    ScoreResult sc;
    SimRng rng;
    size_t i, placed = 0, misplaced = 0, chosen;
    int rc = 1;

    set.mutations = argc > 3 ? (size_t)strtoul(argv[3], NULL, 10) : 20;
    set.coverage = argc > 4 ? atof(argv[4]) : 30.0;
    set.error_rate = argc > 5 ? atof(argv[5]) : 0.01;
    set.seed = argc > 6 ? strtoul(argv[6], NULL, 10) : 42;
    set.min_depth = argc > 7 ? (uint32_t)strtoul(argv[7], NULL, 10) : PILEUP_DEFAULT_MIN_DEPTH;
    set.max_edit = argc > 8 ? (size_t)strtoul(argv[8], NULL, 10) : VERIFY_DEFAULT_MAX_EDIT_DISTANCE;
    K = argc > 9 ? (size_t)strtoul(argv[9], NULL, 10) : 15;
    W = argc > 10 ? (size_t)strtoul(argv[10], NULL, 10) : 10;
    if (K < 4 || K > 32 || W < 1) {
        printf("{\"error\":\"k must be 4..32 and w at least 1\"}\n");
        return 1;
    }

    if (seq_load_fasta(path, &ref) != 0) {
        printf("{\"error\":\"could not load %s\"}\n", path);
        return 1;
    }
    sim_rng_seed(&rng, set.seed);
    if (sim_plant(&ref, set.mutations, &rng, &key) != 0 ||
        sim_make_sample(&ref, &key, &sample) != 0 ||
        sim_reads(&sample, READ_LEN, set.coverage, set.error_rate, &rng, &reads) != 0 ||
        minimizers_generate(&ref, K, W, &ref_mins) != 0 ||
        index_build(&ref, &ref_mins, K, 0, &idx) != 0 ||
        bloom_build_from_minimizers(&ref, &ref_mins, K, BLOOM_FP, &bloom) != 0 ||
        pileup_init(&pile, ref.len) != 0) {
        printf("{\"error\":\"pipeline setup failed\"}\n");
        goto done;
    }

    /* Stages 2 to 8 for every read, remembering where each one was placed. */
    placed_at = (long *)malloc(reads.count * sizeof *placed_at);
    ncand = (size_t *)malloc(reads.count * sizeof *ncand);
    depth = (uint32_t *)malloc(ref.len * sizeof *depth);
    if (placed_at == NULL || ncand == NULL || depth == NULL)
        goto done;
    for (i = 0; i < reads.count; i++) {
        placed_at[i] = place_read(&ref, &reads.reads[i], &idx, &bloom, set.max_edit,
                                  &totals, &ncand[i]);
        if (placed_at[i] < 0)
            continue;
        placed++;
        if ((size_t)placed_at[i] != reads.true_start[i])
            misplaced++;
        pileup_add_read(&pile, &reads.reads[i], (size_t)placed_at[i]);
    }
    if (pileup_depths(&pile, depth) != 0 ||
        pileup_call(&pile, &ref, set.min_depth, PILEUP_DEFAULT_MIN_ALT_PERCENT, &calls) != 0 ||
        score_calls(&calls, &key, depth, ref.len, set.min_depth, &sc) != 0)
        goto done;

    /* Pick the read: as asked, or the first placed read that covers a mutation. */
    chosen = 0;
    if (want >= 0 && (size_t)want < reads.count) {
        chosen = (size_t)want;
    } else {
        for (i = 0; i < reads.count; i++)
            if (placed_at[i] >= 0 && covers_mutation(&key, reads.true_start[i])) {
                chosen = i;
                break;
            }
    }

    printf("{\"settings\":{\"k\":%lu,\"w\":%lu,\"read_len\":%d,\"mutations\":%lu,\"coverage\":%g,"
           "\"error_rate\":%g,\"seed\":%lu,\"band\":%lu,\"max_edit\":%lu,\"min_depth\":%lu,"
           "\"min_alt_pct\":%lu,\"bloom_fp\":%.2f},\n",
           (unsigned long)K, (unsigned long)W, READ_LEN, (unsigned long)set.mutations,
           set.coverage, set.error_rate, set.seed,
           (unsigned long)VERIFY_DEFAULT_BAND, (unsigned long)set.max_edit,
           (unsigned long)set.min_depth, (unsigned long)PILEUP_DEFAULT_MIN_ALT_PERCENT,
           BLOOM_FP);
    printf("\"summary\":{\"ref_len\":%lu,\"ref_nruns\":%lu,\"reads\":%lu,\"placed\":%lu,"
           "\"unplaced\":%lu,\"misplaced\":%lu,\"ref_minimizers\":%lu,\"index_entries\":%lu,"
           "\"nbuckets\":%lu,\"bloom_bits\":%lu,\"bloom_hashes\":%lu,\"calls\":%lu,\"tp\":%lu,"
           "\"fp\":%lu,\"fn\":%lu,\"callable\":%lu,\"precision\":%.1f,\"recall_all\":%.1f,"
           "\"recall_callable\":%.1f,\"words\":%lu,\"bloom_skipped\":%lu,\"lookups\":%lu,\"votes\":%lu},\n",
           (unsigned long)ref.len, (unsigned long)ref.nruns, (unsigned long)reads.count,
           (unsigned long)placed, (unsigned long)(reads.count - placed), (unsigned long)misplaced,
           (unsigned long)ref_mins.count, (unsigned long)idx.nentries, (unsigned long)idx.nbuckets,
           (unsigned long)bloom.bit_count, (unsigned long)bloom.hash_count,
           (unsigned long)calls.count, (unsigned long)sc.tp, (unsigned long)sc.fp,
           (unsigned long)sc.fn, (unsigned long)sc.callable, sc.precision_pct,
           sc.recall_all_pct, sc.recall_callable_pct, (unsigned long)totals.words,
           (unsigned long)totals.bloom_skipped, (unsigned long)totals.lookups,
           (unsigned long)totals.hits);

    /* Reads worth picking: over a mutation, with 2+ candidate starts, and not placed. */
    {
        size_t shown = 0, multi = 0;
        printf("\"interesting\":{\"covers_mutation\":[");
        for (i = 0; i < reads.count && shown < 40; i++)
            if (placed_at[i] >= 0 && covers_mutation(&key, reads.true_start[i]))
                printf("%s%lu", shown++ ? "," : "", (unsigned long)i);
        printf("],\"multi_candidate\":[");
        for (i = 0; i < reads.count && multi < 20; i++)
            if (ncand[i] >= 2)
                printf("%s%lu", multi++ ? "," : "", (unsigned long)i);
        for (i = 0, multi = 0; i < reads.count; i++)
            multi += ncand[i] >= 2;
        printf("],\"multi_candidate_total\":%lu,\"not_placed\":[", (unsigned long)multi);
        for (i = 0, shown = 0; i < reads.count && shown < 15; i++)
            if (placed_at[i] < 0)
                printf("%s%lu", shown++ ? "," : "", (unsigned long)i);
        printf("]},\n");
    }

    /* Stage 1: the read itself, with its planted mutations and sequencing errors marked. */
    {
        const Seq *rd = &reads.reads[chosen];
        size_t s = reads.true_start[chosen], first = 1;

        printf("\"read\":{\"index\":%lu,\"true_start\":%lu,\"placed_start\":%ld,\"letters\":",
               (unsigned long)chosen, (unsigned long)s, placed_at[chosen]);
        print_letters(rd, 0, rd->len);
        printf(",\"ref_letters\":");
        print_letters(&ref, s, rd->len);
        printf(",\"mutations\":[");
        for (i = 0; i < key.count; i++)
            if (key.items[i].pos >= s && key.items[i].pos < s + rd->len) {
                printf("%s%lu", first ? "" : ",", (unsigned long)(key.items[i].pos - s));
                first = 0;
            }
        printf("],\"errors\":[");
        for (i = 0, first = 1; i < rd->len; i++)
            if (seq_get(rd, i) != seq_get(&sample, s + i)) {
                printf("%s%lu", first ? "" : ",", (unsigned long)i);
                first = 0;
            }
        printf("],\"packed_bytes\":[");
        for (i = 0; i < (rd->len + 3) / 4; i++)
            printf("%s%u", i ? "," : "", (unsigned)rd->cells[i]);
        printf("]},\n");
    }

    /* Stages 2 to 7 for the chosen read, step by step. */
    {
        const Seq *rd = &reads.reads[chosen];
        MinimizerSet all = {0}, mins = {0};
        VoteTree tree = {0};
        Top3Heap top;
        OfferTrace ot;
        TopCandidate sorted[TOP3_CAPACITY];
        size_t n = 0, c, bits[MAX_BLOOM_HASHES];
        VerifyResult vr;
        int bloom_ok = 1;

        /* Window of 1 k-mer keeps every k-mer, giving every hash with stage 2's own function. */
        if (minimizers_generate(rd, K, 1, &all) != 0 || minimizers_generate(rd, K, W, &mins) != 0)
            goto done;
        printf("\"stage2\":{\"kmers\":[");
        for (i = 0; i < all.count; i++) {
            printf("%s{\"pos\":%lu,\"hash\":", i ? "," : "", (unsigned long)all.items[i].pos);
            print_hash(all.items[i].hash);
            putchar('}');
        }
        printf("],\"minimizers\":[");
        for (i = 0; i < mins.count; i++)
            printf("%s%lu", i ? "," : "", (unsigned long)mins.items[i].pos);
        printf("]},\n");

        /* Stages 3 to 5: Bloom check, index lookup, votes, word by word (same order as vote_read). */
        printf("\"words\":[");
        for (i = 0; i < mins.count; i++) {
            size_t pos = mins.items[i].pos, probes = 0;
            uint64_t hash = mins.items[i].hash;
            uint64_t wid = word_id_at(rd, pos, K);
            int maybe = bloom_maybe_contains(&bloom, hash, wid);
            const IndexEntry *e = NULL;
            const PosNode *p;

            bloom_ok &= bloom_positions(&bloom, hash, wid, bits);
            printf("%s{\"pos\":%lu,\"word\":", i ? "," : "", (unsigned long)pos);
            print_letters(rd, pos, K);
            printf(",\"hash\":");
            print_hash(hash);
            printf(",\"bloom\":%d,\"bits\":[", maybe);
            for (c = 0; c < bloom.hash_count && c < MAX_BLOOM_HASHES; c++)
                printf("%s%lu", c ? "," : "", (unsigned long)bits[c]);
            printf("],\"bucket\":%lu", (unsigned long)(hash % (uint64_t)idx.nbuckets));
            if (maybe == 1)
                e = index_lookup(&idx, hash, wid, &probes);
            printf(",\"looked_up\":%d,\"probes\":%lu,\"ref_positions\":[", maybe == 1,
                   (unsigned long)probes);
            for (p = e ? e->pos_head : NULL, c = 0; p != NULL; p = p->next, c++) {
                printf("%s%lu", c ? "," : "", (unsigned long)p->pos);
                vote_tree_add(&tree, (long)p->pos - (long)pos);
            }
            printf("],\"votes\":[");
            for (p = e ? e->pos_head : NULL, c = 0; p != NULL; p = p->next, c++)
                printf("%s%ld", c ? "," : "", (long)p->pos - (long)pos);
            printf("]}");
        }
        printf("],\"bloom_display_matches_c\":%d,\n", bloom_ok);

        printf("\"stage5\":{\"nodes\":%lu,\"height\":%lu,\"tree\":",
               (unsigned long)tree.nodes, (unsigned long)vote_tree_height(&tree));
        print_tree(tree.root);
        printf("},\n");

        /* Stage 6: the same in-order walk as top3_from_vote_tree, one offer at a time. */
        top3_init(&top);
        ot.top = &top;
        ot.first = 1;
        printf("\"stage6\":{\"offers\":[");
        vote_tree_walk(&tree, trace_offer, &ot);
        top3_get_sorted(&top, sorted, &n);
        printf("],\"sorted\":[");
        for (c = 0; c < n; c++)
            printf("%s[%ld,%lu]", c ? "," : "", sorted[c].start, (unsigned long)sorted[c].votes);
        printf("]},\n");

        /* Stage 7: verify.c's decision, plus the band of the table for each candidate. */
        verify_candidates(&ref, rd, sorted, n, VERIFY_DEFAULT_BAND, set.max_edit, &vr);
        printf("\"stage7\":{\"mapped\":%d,\"start\":%ld,\"candidates\":[", vr.mapped,
               vr.mapped ? vr.start : -1L);
        for (c = 0; c < vr.candidate_count; c++) {
            const VerifyCandidateResult *r = &vr.candidates[c];
            printf("%s{\"start\":%ld,\"votes\":%lu,\"status\":\"%s\",\"accepted\":%d,\"distance\":",
                   c ? "," : "", r->start, (unsigned long)r->votes, status_name(r->status),
                   r->accepted);
            if (r->edit_distance == VERIFY_DISTANCE_INF)
                printf("null");
            else
                printf("%lu", (unsigned long)r->edit_distance);
            if (r->start >= 0 && (size_t)r->start < ref.len) {
                size_t avail = ref.len - (size_t)r->start;
                size_t win = rd->len < avail ? rd->len : avail;
                printf(",\"window\":");
                print_letters(&ref, (size_t)r->start, win);
                printf(",\"dp\":");
                if (r->status == VERIFY_CANDIDATE_EVALUATED)
                    print_dp_band(rd, &ref, (size_t)r->start, win, VERIFY_DEFAULT_BAND);
                else
                    printf("null");
            }
            putchar('}');
        }
        printf("]},\n");

        vote_tree_clear(&tree);
        minimizers_free(&all);
        minimizers_free(&mins);
    }

    /* Stage 3's whole table: every entry in bucket order, then chain order. */
    printf("\"index\":[");
    for (i = 0, rc = 0; i < idx.nbuckets; i++) {      /* rc doubles as "printed one already" */
        const IndexEntry *e;
        for (e = idx.buckets[i]; e != NULL; e = e->next) {
            const PosNode *p;
            size_t c = 0;
            printf("%s[%lu,", rc ? "," : "", (unsigned long)i);
            print_letters(&ref, e->pos_head->pos, K);
            printf(",[");
            for (p = e->pos_head; p != NULL; p = p->next)
                printf("%s%lu", c++ ? "," : "", (unsigned long)p->pos);
            printf("]]");
            rc = 1;
        }
    }
    rc = 1;                                            /* back to "failed" until the end */
    printf("],\n");

    /* Stage 4's whole bit array as hex, byte by byte (bit 0 of a byte = lowest bit). */
    printf("\"bloom\":{\"m\":%lu,\"h\":%lu,\"hex\":\"", (unsigned long)bloom.bit_count,
           (unsigned long)bloom.hash_count);
    for (i = 0; i < bloom.byte_count; i++)
        printf("%02x", (unsigned)bloom.bits[i]);
    printf("\"},\n");

    /* Stage 8: the pile-up under this read (all placed reads counted). */
    {
        size_t s = placed_at[chosen] >= 0 ? (size_t)placed_at[chosen] : reads.true_start[chosen];
        size_t len = READ_LEN < ref.len - s ? READ_LEN : ref.len - s;
        size_t first = 1, shown;

        printf("\"stage8\":{\"span_start\":%lu,\"depth_before\":%lu,\"ref_letters\":",
               (unsigned long)s, (unsigned long)(s > 0 ? depth[s - 1] : 0));
        print_letters(&ref, s, len);
        printf(",\"depth\":[");
        for (i = 0; i < len; i++)
            printf("%s%lu", i ? "," : "", (unsigned long)depth[s + i]);
        printf("],\"diff\":[");
        for (i = 0; i < len; i++)
            printf("%s%ld", i ? "," : "", (long)pile.diff[s + i]);
        printf("],\"counts\":[");
        for (i = 0; i < len; i++) {
            const uint32_t *ct = pile.cells[s + i].count;
            printf("%s[%lu,%lu,%lu,%lu]", i ? "," : "", (unsigned long)ct[0],
                   (unsigned long)ct[1], (unsigned long)ct[2], (unsigned long)ct[3]);
        }
        printf("],\"calls\":[");
        for (i = 0; i < calls.count; i++)
            if (calls.items[i].pos >= s && calls.items[i].pos < s + len) {
                printf("%s%lu", first ? "" : ",", (unsigned long)(calls.items[i].pos - s));
                first = 0;
            }
        /* The placed reads that overlap this span (up to 60), so the page can stack them. */
        printf("],\"reads\":[");
        for (i = 0, shown = 0; i < reads.count && shown < 60; i++) {
            long ps = placed_at[i];
            if (ps < 0 || (size_t)ps >= s + len || (size_t)ps + READ_LEN <= s)
                continue;
            printf("%s{\"index\":%lu,\"start\":%ld,\"letters\":", shown ? "," : "",
                   (unsigned long)i, ps);
            print_letters(&reads.reads[i], 0, reads.reads[i].len);
            putchar('}');
            shown++;
        }
        printf("]},\n");
    }

    /* Stage 9: the whole key and every call, for the genome-wide picture. */
    printf("\"stage9\":{\"key\":[");
    for (i = 0; i < key.count; i++) {
        /* Context letters, pile-up counts and one placed read over it, for the page's zoom and jump. */
        size_t p = key.items[i].pos, from = p >= 12 ? p - 12 : 0, r;
        size_t to = p + 13 < ref.len ? p + 13 : ref.len;
        const uint32_t *ct = pile.cells[p].count;
        long over = -1;
        for (r = 0; r < reads.count && over < 0; r++)
            if (placed_at[r] >= 0 && (size_t)placed_at[r] <= p && (size_t)placed_at[r] + READ_LEN > p)
                over = (long)r;
        printf("%s{\"pos\":%lu,\"ref\":\"%c\",\"alt\":\"%c\",\"depth\":%lu,\"read\":%ld,"
               "\"context_start\":%lu,\"context\":", i ? "," : "", (unsigned long)p,
               LETTER[key.items[i].ref], LETTER[key.items[i].alt], (unsigned long)depth[p], over,
               (unsigned long)from);
        print_letters(&ref, from, to - from);
        printf(",\"counts\":[%lu,%lu,%lu,%lu]}", (unsigned long)ct[0], (unsigned long)ct[1],
               (unsigned long)ct[2], (unsigned long)ct[3]);
    }
    printf("],\"calls\":[");
    for (i = 0; i < calls.count; i++) {
        const MutationCall *m = &calls.items[i];
        const KeyEntry *e = key_find(&key, m->pos);
        printf("%s{\"pos\":%lu,\"ref\":\"%c\",\"alt\":\"%c\",\"alt_count\":%lu,\"depth\":%lu,\"tp\":%d}",
               i ? "," : "", (unsigned long)m->pos, LETTER[m->ref], LETTER[m->alt],
               (unsigned long)m->alt_count, (unsigned long)m->depth,
               e != NULL && e->alt == m->alt);
    }
    printf("],\"depth_sample\":[");
    for (i = 0; i < ref.len; i += 100)
        printf("%s%lu", i ? "," : "", (unsigned long)depth[i]);
    /* The whole reference (about 30 KB), so the page's helix shows real letters at any zoom. */
    printf("],\"genome\":");
    print_letters(&ref, 0, ref.len);
    printf("}}\n");
    rc = 0;

done:
    free(placed_at);
    free(ncand);
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
