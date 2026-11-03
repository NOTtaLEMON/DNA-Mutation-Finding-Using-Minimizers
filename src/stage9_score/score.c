#include "score.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char LETTER[] = "ACGT";

/* Appends one entry, doubling capacity when full. See score.h. */
int key_add(AnswerKey *key, size_t pos, uint8_t ref, uint8_t alt)
{
    if (key == NULL || ref > 3 || alt > 3 || ref == alt)
        return -1;
    if (key->count == key->capacity) {
        size_t new_cap = key->capacity ? key->capacity * 2 : 16;
        KeyEntry *grown = (KeyEntry *)realloc(key->items, new_cap * sizeof *grown);
        if (grown == NULL)
            return -1;
        key->items = grown;
        key->capacity = new_cap;
    }
    key->items[key->count].pos = pos;
    key->items[key->count].ref = ref;
    key->items[key->count].alt = alt;
    key->count++;
    return 0;
}

/* qsort comparison: orders entries by position. */
static int compare_entries(const void *a, const void *b)
{
    size_t pa = ((const KeyEntry *)a)->pos;
    size_t pb = ((const KeyEntry *)b)->pos;
    return (pa > pb) - (pa < pb);
}

/* Sorts by position, then checks neighbours for duplicates. See score.h. */
int key_sort(AnswerKey *key)
{
    size_t i;

    if (key == NULL)
        return -1;
    if (key->count > 1)
        qsort(key->items, key->count, sizeof *key->items, compare_entries);
    for (i = 1; i < key->count; i++)
        if (key->items[i].pos == key->items[i - 1].pos)
            return -1;
    return 0;
}

/* Binary search on [lo, hi): halve the range until pos is found or the range is empty. See score.h. */
const KeyEntry *key_find(const AnswerKey *key, size_t pos)
{
    size_t lo = 0, hi;

    if (key == NULL)
        return NULL;
    hi = key->count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (key->items[mid].pos == pos)
            return &key->items[mid];
        if (key->items[mid].pos < pos)
            lo = mid + 1;
        else
            hi = mid;
    }
    return NULL;
}

/* Writes the header line and one 1-based "position ref alt" line per entry. See score.h. */
int key_save(const char *path, const AnswerKey *key)
{
    FILE *f;
    size_t i;

    if (path == NULL || key == NULL)
        return -1;
    f = fopen(path, "w");
    if (f == NULL)
        return -1;
    fprintf(f, "# position(1-based) ref alt\n");
    for (i = 0; i < key->count; i++)
        fprintf(f, "%lu %c %c\n", (unsigned long)(key->items[i].pos + 1),
                LETTER[key->items[i].ref], LETTER[key->items[i].alt]);
    return fclose(f) == 0 ? 0 : -1;
}

/* Parses "position ref alt" lines, skipping '#' lines and blank lines. See score.h. */
int key_load(const char *path, AnswerKey *key)
{
    FILE *f;
    char line[256];
    int ok = 1;

    if (path == NULL || key == NULL)
        return -1;
    memset(key, 0, sizeof *key);
    f = fopen(path, "r");
    if (f == NULL)
        return -1;

    while (ok && fgets(line, sizeof line, f) != NULL) {
        unsigned long pos;
        char ref, alt;
        int r, a;

        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
            continue;
        if (sscanf(line, "%lu %c %c", &pos, &ref, &alt) != 3 || pos == 0) {
            ok = 0;
            break;
        }
        r = base_to_code(ref);
        a = base_to_code(alt);
        if (r < 0 || a < 0 || key_add(key, (size_t)(pos - 1), (uint8_t)r, (uint8_t)a) != 0)
            ok = 0;
    }
    fclose(f);

    if (!ok || key_sort(key) != 0) {
        key_free(key);
        return -1;
    }
    return 0;
}

/* Frees the key and empties it. See score.h. */
void key_free(AnswerKey *key)
{
    if (key == NULL)
        return;
    free(key->items);
    memset(key, 0, sizeof *key);
}

/* Returns 100 * part / whole, or 0 when whole is 0. */
static double percent(size_t part, size_t whole)
{
    return whole == 0 ? 0.0 : 100.0 * (double)part / (double)whole;
}

/* Binary-searches each call in the key, then counts callable planted positions. See score.h. */
int score_calls(const CallList *calls, const AnswerKey *key,
                const uint32_t *depth, size_t ref_len, uint32_t min_depth,
                ScoreResult *out)
{
    size_t i;

    if (out == NULL)
        return -1;
    memset(out, 0, sizeof *out);
    if (calls == NULL || key == NULL || (calls->count > 0 && calls->items == NULL))
        return -1;

    out->calls = calls->count;
    out->planted = key->count;

    for (i = 0; i < calls->count; i++) {
        const KeyEntry *e = key_find(key, calls->items[i].pos);
        if (e == NULL)
            continue;
        if (e->alt == calls->items[i].alt)
            out->tp++;
        else
            out->wrong_alt++;
    }
    out->fp = out->calls - out->tp;
    out->fn = out->planted - out->tp;

    if (depth != NULL) {
        for (i = 0; i < key->count; i++)
            if (key->items[i].pos < ref_len && depth[key->items[i].pos] >= min_depth)
                out->callable++;
    }

    out->precision_pct = percent(out->tp, out->calls);
    out->recall_all_pct = percent(out->tp, out->planted);
    out->recall_callable_pct = percent(out->tp, out->callable);
    return 0;
}
