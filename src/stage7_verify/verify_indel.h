/*
 * INACTIVE: indel-aware version of Stage 7 (kept for the paper track).
 *
 * Everything in this file is switched off unless the compiler is given
 * -DVERIFY_INDEL_ENABLED. The normal build and the Nov 4 submission use
 * verify.h / verify.c only. #ifdef is used instead of a comment block
 * because the code has its own comments, which a block comment cannot hold.
 *
 * What is different from verify.c:
 *   1. The reference window is read_length + band letters long, and the
 *      alignment may end anywhere in the last band of columns. A read with an
 *      extra or missing letter no longer pays a second edit at the window end,
 *      so one insertion or deletion costs 1 (like a substitution).
 *   2. A traceback through the DP table records the alignment path as a list
 *      of operations, so Stage 8 can place every read letter at the right
 *      reference position even after an indel:
 *        'M' read letter sits on a reference letter (match or substitution)
 *        'I' read letter has no reference letter (insertion in the sample)
 *        'D' reference letter has no read letter (deletion in the sample)
 */
#ifndef VERIFY_INDEL_H
#define VERIFY_INDEL_H

#ifdef VERIFY_INDEL_ENABLED

#include "verify.h"

/* Alignment operation codes stored in IndelCandidateResult.ops. */
#define ALIGN_MATCH  'M'
#define ALIGN_INSERT 'I'
#define ALIGN_DELETE 'D'

/*
 * Verification result for one Stage 6 candidate.
 *
 * start, votes, edit_distance, accepted, status: same meaning as in
 *                verify.h.
 * ref_span:      number of reference letters the alignment covers, starting
 *                at start ('M' plus 'D' operations).
 * ops / nops:    the alignment path, first read letter first. Allocated only
 *                for accepted candidates; NULL otherwise.
 */
typedef struct {
    long start;
    size_t votes;
    size_t edit_distance;
    size_t ref_span;
    char *ops;
    size_t nops;
    int accepted;
    VerifyCandidateStatus status;
} IndelCandidateResult;

/*
 * Final result for one read.
 *
 * mapped == 1 means candidates[best] passed the threshold and is the chosen
 * placement (same selection rule as verify.c: smallest distance, then more
 * votes, then smaller start). Free with verify_indel_result_free().
 */
typedef struct {
    int mapped;
    size_t best;
    size_t candidate_count;
    IndelCandidateResult candidates[TOP3_CAPACITY];
} IndelVerifyResult;

/*
 * Indel-aware replacement for verify_candidates(). Same inputs, same N and
 * start checks. Returns 0 on success (including when nothing is accepted),
 * -1 for invalid arguments or allocation failure.
 */
int verify_candidates_indel(const Seq *reference,
                            const Seq *read,
                            const TopCandidate *candidates,
                            size_t candidate_count,
                            size_t band,
                            size_t max_edit_distance,
                            IndelVerifyResult *out);

/* Frees every ops array inside out and resets it. Safe on a zeroed result. */
void verify_indel_result_free(IndelVerifyResult *out);

#endif /* VERIFY_INDEL_ENABLED */

#endif /* VERIFY_INDEL_H */
