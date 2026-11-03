#include "pileup.h"

#include <stdlib.h>
#include <string.h>

/* Allocates zeroed cells and a zeroed difference array. See pileup.h. */
int pileup_init(Pileup *pile, size_t ref_len)
{
    if (pile == NULL)
        return -1;
    memset(pile, 0, sizeof *pile);
    if (ref_len == 0)
        return -1;

    pile->cells = (PileCell *)calloc(ref_len, sizeof *pile->cells);
    pile->diff = (int64_t *)calloc(ref_len + 1, sizeof *pile->diff);
    if (pile->cells == NULL || pile->diff == NULL) {
        pileup_free(pile);
        return -1;
    }
    pile->len = ref_len;
    return 0;
}

/* One counter up per read letter, then +1 / -1 at the covered range's ends. See pileup.h. */
int pileup_add_read(Pileup *pile, const Seq *read, size_t start)
{
    size_t covered;
    size_t i;

    if (pile == NULL || pile->cells == NULL || read == NULL || read->len == 0)
        return -1;
    if (start >= pile->len || read->nruns > 0)
        return -1;

    /* Clip at the reference end: only letters with start + i < len count. */
    covered = read->len < pile->len - start ? read->len : pile->len - start;

    for (i = 0; i < covered; i++)
        pile->cells[start + i].count[seq_get(read, i)]++;

    pile->diff[start] += 1;
    pile->diff[start + covered] -= 1;
    pile->reads++;
    return 0;
}

/* Running sum of the difference array gives the depth at each position. See pileup.h. */
int pileup_depths(const Pileup *pile, uint32_t *depth_out)
{
    int64_t running = 0;
    size_t p;

    if (pile == NULL || pile->diff == NULL || depth_out == NULL)
        return -1;
    for (p = 0; p < pile->len; p++) {
        running += pile->diff[p];
        depth_out[p] = (uint32_t)running;
    }
    return 0;
}

/* Appends one call, doubling the list's capacity when it is full. Returns 0 or -1. */
static int push_call(CallList *calls, const MutationCall *call)
{
    if (calls->count == calls->capacity) {
        size_t new_cap = calls->capacity ? calls->capacity * 2 : 16;
        MutationCall *grown = (MutationCall *)realloc(calls->items,
                                                      new_cap * sizeof *grown);
        if (grown == NULL)
            return -1;
        calls->items = grown;
        calls->capacity = new_cap;
    }
    calls->items[calls->count++] = *call;
    return 0;
}

/* Scans positions in order, keeping a running depth and a pointer into the reference N runs. See pileup.h. */
int pileup_call(const Pileup *pile, const Seq *reference,
                uint32_t min_depth, uint32_t min_alt_percent, CallList *out)
{
    int64_t running = 0;
    size_t run = 0;                       /* next reference N run not yet passed */
    size_t p;

    if (out == NULL)
        return -1;
    out->count = 0;

    if (pile == NULL || pile->cells == NULL || reference == NULL ||
        reference->len != pile->len || min_alt_percent > 100)
        return -1;

    for (p = 0; p < pile->len; p++) {
        const PileCell *cell = &pile->cells[p];
        uint8_t ref_code;
        uint8_t alt = 0;
        uint32_t alt_count = 0;
        uint32_t depth;
        uint8_t c;

        running += pile->diff[p];
        depth = (uint32_t)running;

        /* Skip positions inside a reference N run: there is no reference letter to differ from. */
        while (run < reference->nruns &&
               reference->runs[run].start + reference->runs[run].length <= p)
            run++;
        if (run < reference->nruns && reference->runs[run].start <= p)
            continue;

        if (depth == 0 || depth < min_depth)
            continue;

        /* Most common non-reference letter; strict > keeps the smaller code on a tie. */
        ref_code = seq_get(reference, p);
        for (c = 0; c < 4; c++) {
            if (c != ref_code && cell->count[c] > alt_count) {
                alt = c;
                alt_count = cell->count[c];
            }
        }
        if (alt_count == 0)
            continue;

        if ((uint64_t)alt_count * 100u >= (uint64_t)min_alt_percent * depth) {
            MutationCall call;
            call.pos = p;
            call.ref = ref_code;
            call.alt = alt;
            call.alt_count = alt_count;
            call.depth = depth;
            if (push_call(out, &call) != 0)
                return -1;
        }
    }
    return 0;
}

/* Frees the call list and empties it. See pileup.h. */
void calls_free(CallList *calls)
{
    if (calls == NULL)
        return;
    free(calls->items);
    memset(calls, 0, sizeof *calls);
}

/* Frees the pile-up and empties it. See pileup.h. */
void pileup_free(Pileup *pile)
{
    if (pile == NULL)
        return;
    free(pile->cells);
    free(pile->diff);
    memset(pile, 0, sizeof *pile);
}
