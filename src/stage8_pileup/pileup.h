#ifndef PILEUP_H
#define PILEUP_H

#include <stddef.h>
#include <stdint.h>

#include "seq.h"

/*
 * Stage 8: pile up all placed reads and call mutations.
 *
 * Every read that Stage 7 accepts is added with its start. For each read
 * letter i, the counter for that letter at reference position start + i goes
 * up by one (substitution-only rule; see verify_indel.h for the inactive
 * indel path). Depth (reads covering a position) comes from a difference
 * array, so it stays correct once deletions are added later.
 *
 * Call rule (toy defaults = JS reference): at each position take the most
 * common letter that is not the reference letter. Call it when
 *   depth >= min_depth  and  alt_count * 100 >= min_alt_percent * depth.
 */
#define PILEUP_DEFAULT_MIN_DEPTH ((uint32_t)3)
#define PILEUP_DEFAULT_MIN_ALT_PERCENT ((uint32_t)50)

/* The four letter counters for one reference position (index = 2-bit code). */
typedef struct {
    uint32_t count[4];
} PileCell;

/*
 * cells:  one PileCell per reference position.
 * diff:   difference array of length len + 1; the running sum of
 *         diff[0..p] is the depth at position p.
 * len:    reference length.
 * reads:  number of reads added.
 */
typedef struct {
    PileCell *cells;
    int64_t *diff;
    size_t len;
    size_t reads;
} Pileup;

/* One called mutation. ref and alt are 2-bit codes (A=0, C=1, G=2, T=3). */
typedef struct {
    size_t pos;
    uint8_t ref;
    uint8_t alt;
    uint32_t alt_count;
    uint32_t depth;
} MutationCall;

/* Growable list of calls, in increasing position order. */
typedef struct {
    MutationCall *items;
    size_t count;
    size_t capacity;
} CallList;

/* Allocates an empty pile-up for a reference of ref_len letters. Returns 0, or -1 on failure. */
int pileup_init(Pileup *pile, size_t ref_len);

/*
 * Adds one accepted read placed at start. Letters past the reference end are
 * clipped. Returns -1 if start is outside the reference or the read contains
 * an N (Stage 7 never places such reads).
 */
int pileup_add_read(Pileup *pile, const Seq *read, size_t start);

/* Writes the depth of every position into depth_out (pile->len entries). Returns 0 or -1. */
int pileup_depths(const Pileup *pile, uint32_t *depth_out);

/*
 * Scans every position and appends calls to out (which is reset first).
 * Positions where the reference letter is N are skipped. When two non-reference
 * letters tie, the smaller code wins (A, C, G, T order) so results are
 * deterministic. Returns 0, or -1 on bad arguments or allocation failure.
 */
int pileup_call(const Pileup *pile, const Seq *reference,
                uint32_t min_depth, uint32_t min_alt_percent, CallList *out);

/* Frees the call list and empties it. */
void calls_free(CallList *calls);

/* Frees the pile-up and empties it. */
void pileup_free(Pileup *pile);

#endif /* PILEUP_H */
