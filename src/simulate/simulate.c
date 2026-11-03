#include "simulate.h"

#include <stdlib.h>
#include <string.h>

/* Mixes the seed so nearby seeds give unrelated streams, and avoids the all-zero state. See simulate.h. */
void sim_rng_seed(SimRng *rng, uint64_t seed)
{
    uint64_t z = seed + 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    rng->state = z ? z : 1;
}

/* xorshift64*: three shifts and a multiply. See simulate.h. */
uint64_t sim_rng_next(SimRng *rng)
{
    uint64_t x = rng->state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    rng->state = x;
    return x * 0x2545F4914F6CDD1Dull;
}

/* Random number in [0, n). The tiny modulo bias is irrelevant for test data. See simulate.h. */
size_t sim_rng_below(SimRng *rng, size_t n)
{
    return (size_t)(sim_rng_next(rng) % (uint64_t)n);
}

/* Top 53 bits as a fraction in [0, 1). See simulate.h. */
double sim_rng_unit(SimRng *rng)
{
    return (double)(sim_rng_next(rng) >> 11) / 9007199254740992.0;
}

/* Allocates an all-A Seq of len letters with no N runs. Returns 0 or -1. */
static int seq_alloc(Seq *s, size_t len)
{
    memset(s, 0, sizeof *s);
    s->cells_cap = len / 4 + 1;
    s->cells = (uint8_t *)calloc(s->cells_cap, 1);
    if (s->cells == NULL)
        return -1;
    s->len = len;
    return 0;
}

/* Overwrites the letter at pos with code, using stage 1's packing (top bits first). */
static void seq_put(Seq *s, size_t pos, uint8_t code)
{
    unsigned int shift = 6u - 2u * (unsigned int)(pos % 4);
    s->cells[pos / 4] = (uint8_t)((s->cells[pos / 4] & ~(0x03u << shift)) | (code << shift));
}

/* Returns 1 if [start, start+len) overlaps an N run of s, else 0. Runs are sorted. */
static int overlaps_n(const Seq *s, size_t start, size_t len)
{
    size_t i;

    for (i = 0; i < s->nruns; i++) {
        if (s->runs[i].start >= start + len)
            break;
        if (s->runs[i].start + s->runs[i].length > start)
            return 1;
    }
    return 0;
}

/* Marks taken positions in a byte array and re-draws on a clash. See simulate.h. */
int sim_plant(const Seq *reference, size_t count, SimRng *rng, AnswerKey *key)
{
    uint8_t *taken;
    size_t usable = 0, i, placed = 0;

    if (reference == NULL || rng == NULL || key == NULL || reference->len == 0)
        return -1;
    key->count = 0;

    for (i = 0; i < reference->nruns; i++)
        usable += reference->runs[i].length;
    usable = reference->len - usable;
    if (count > usable)
        return -1;

    taken = (uint8_t *)calloc(reference->len, 1);
    if (taken == NULL)
        return -1;
    for (i = 0; i < reference->nruns; i++)
        memset(taken + reference->runs[i].start, 1, reference->runs[i].length);

    while (placed < count) {
        size_t pos = sim_rng_below(rng, reference->len);
        uint8_t ref, alt;

        if (taken[pos])
            continue;
        taken[pos] = 1;
        ref = seq_get(reference, pos);
        alt = (uint8_t)((ref + 1 + sim_rng_below(rng, 3)) % 4);   /* one of the 3 other letters */
        if (key_add(key, pos, ref, alt) != 0) {
            free(taken);
            return -1;
        }
        placed++;
    }
    free(taken);
    return key_sort(key);
}

/* Copies cells and N runs, then writes each planted alt letter. See simulate.h. */
int sim_make_sample(const Seq *reference, const AnswerKey *key, Seq *sample)
{
    size_t i;

    if (reference == NULL || key == NULL || sample == NULL)
        return -1;
    if (seq_alloc(sample, reference->len) != 0)
        return -1;
    memcpy(sample->cells, reference->cells, (reference->len + 3) / 4);

    if (reference->nruns > 0) {
        sample->runs = (NRun *)malloc(reference->nruns * sizeof *sample->runs);
        if (sample->runs == NULL) {
            seq_free(sample);
            return -1;
        }
        memcpy(sample->runs, reference->runs, reference->nruns * sizeof *sample->runs);
        sample->nruns = sample->runs_cap = reference->nruns;
    }

    for (i = 0; i < key->count; i++) {
        if (key->items[i].pos >= sample->len) {
            seq_free(sample);
            return -1;
        }
        seq_put(sample, key->items[i].pos, key->items[i].alt);
    }
    return 0;
}

/* Draws a start, copies read_len letters, and flips each with probability error_rate. See simulate.h. */
int sim_reads(const Seq *sample, size_t read_len, double coverage,
              double error_rate, SimRng *rng, ReadSet *out)
{
    size_t n, r;

    if (out == NULL)
        return -1;
    memset(out, 0, sizeof *out);
    if (sample == NULL || rng == NULL || read_len == 0 || read_len > sample->len ||
        coverage <= 0.0 || error_rate < 0.0 || error_rate >= 1.0)
        return -1;

    n = (size_t)(coverage * (double)sample->len / (double)read_len + 0.5);
    if (n == 0)
        n = 1;
    out->reads = (Seq *)calloc(n, sizeof *out->reads);
    out->true_start = (size_t *)calloc(n, sizeof *out->true_start);
    if (out->reads == NULL || out->true_start == NULL) {
        sim_reads_free(out);
        return -1;
    }

    for (r = 0; r < n; r++) {
        size_t start, i, tries = 0;
        Seq *read = &out->reads[r];

        do {                                       /* stage 7 never places reads with N */
            start = sim_rng_below(rng, sample->len - read_len + 1);
        } while (overlaps_n(sample, start, read_len) && ++tries < 1000);
        if (tries >= 1000 || seq_alloc(read, read_len) != 0) {
            sim_reads_free(out);
            return -1;
        }
        out->count = r + 1;
        out->true_start[r] = start;

        for (i = 0; i < read_len; i++) {
            uint8_t code = seq_get(sample, start + i);
            if (sim_rng_unit(rng) < error_rate)
                code = (uint8_t)((code + 1 + sim_rng_below(rng, 3)) % 4);
            seq_put(read, i, code);
        }
    }
    return 0;
}

/* Frees every read and the arrays, and empties the set. See simulate.h. */
void sim_reads_free(ReadSet *set)
{
    size_t i;

    if (set == NULL)
        return;
    for (i = 0; i < set->count; i++)
        seq_free(&set->reads[i]);
    free(set->reads);
    free(set->true_start);
    memset(set, 0, sizeof *set);
}
