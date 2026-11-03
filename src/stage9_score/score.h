#ifndef SCORE_H
#define SCORE_H

#include <stddef.h>
#include <stdint.h>

#include "pileup.h"

/*
 * Stage 9: grade mutation calls against the planted answer key (test data only).
 *
 * The key is a sorted array; each call is looked up by binary search on
 * position. A call is a true positive when the key has that position AND the
 * same alt letter. Positions are 0-based in memory and 1-based in key files
 * and printed output (like VCF/SAM).
 */

/* One planted mutation. ref and alt are 2-bit codes (A=0, C=1, G=2, T=3). */
typedef struct {
    size_t pos;
    uint8_t ref;
    uint8_t alt;
} KeyEntry;

/* Growable array of planted mutations; sorted by pos after key_sort(). */
typedef struct {
    KeyEntry *items;
    size_t count;
    size_t capacity;
} AnswerKey;

/*
 * tp:          calls with the right position and the right alt letter.
 * fp:          all other calls (calls - tp).
 * fn:          planted mutations not called correctly (key size - tp).
 * wrong_alt:   calls at a planted position but with a different alt (counted in fp).
 * callable:    planted mutations whose depth is >= min_depth, i.e. the pile-up
 *              had enough reads to call them at all.
 * The *_pct fields are percentages; 0 when the denominator is 0.
 */
typedef struct {
    size_t calls;
    size_t planted;
    size_t tp;
    size_t fp;
    size_t fn;
    size_t wrong_alt;
    size_t callable;
    double precision_pct;
    double recall_all_pct;
    double recall_callable_pct;
} ScoreResult;

/* Appends one entry (0-based pos). Returns 0, or -1 on bad codes or no memory. */
int key_add(AnswerKey *key, size_t pos, uint8_t ref, uint8_t alt);

/* Sorts the key by position. Returns -1 if two entries share a position. */
int key_sort(AnswerKey *key);

/* Binary search for pos in a sorted key. Returns the entry, or NULL. */
const KeyEntry *key_find(const AnswerKey *key, size_t pos);

/*
 * Writes the key as text: a '#' header line, then one "position ref alt"
 * line per entry, position 1-based, letters as A/C/G/T. Returns 0 or -1.
 */
int key_save(const char *path, const AnswerKey *key);

/*
 * Reads a file written by key_save (lines starting with '#' are skipped),
 * converts positions to 0-based and sorts the key. Returns 0, or -1 on a
 * malformed line, a duplicate position, or a file error.
 */
int key_load(const char *path, AnswerKey *key);

/* Frees the key and empties it. */
void key_free(AnswerKey *key);

/*
 * Scores calls against a sorted key. depth (one entry per reference
 * position, from pileup_depths) is used for the callable count; pass NULL
 * and ref_len 0 to skip it. Returns 0, or -1 on bad arguments.
 */
int score_calls(const CallList *calls, const AnswerKey *key,
                const uint32_t *depth, size_t ref_len, uint32_t min_depth,
                ScoreResult *out);

#endif /* SCORE_H */
