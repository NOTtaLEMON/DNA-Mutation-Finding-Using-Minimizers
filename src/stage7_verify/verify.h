#ifndef VERIFY_H
#define VERIFY_H

#include <stddef.h>
#include <stdint.h>

#include "seq.h"
#include "top3.h"

/*
 * Stage 7 defaults from the project design.
 * The toy configuration uses a DP band of 2 and accepts edit distance <= 3.
 * For real experiments, make these configurable and benchmark them; they are
 * not universal biological constants.
 */
#define VERIFY_DEFAULT_BAND ((size_t)2)
#define VERIFY_DEFAULT_MAX_EDIT_DISTANCE ((size_t)3)

/*
 * Public representation of an unreachable / unusable candidate distance.
 * Do not compare this sentinel to a threshold without checking it first.
 */
#define VERIFY_DISTANCE_INF SIZE_MAX

/* Why a candidate did or did not receive a finite DP comparison. */
typedef enum {
    VERIFY_CANDIDATE_EVALUATED = 0,
    VERIFY_CANDIDATE_NEGATIVE_START,
    VERIFY_CANDIDATE_START_OUT_OF_RANGE,
    VERIFY_CANDIDATE_UNKNOWN_BASE,
    VERIFY_CANDIDATE_OUTSIDE_BAND
} VerifyCandidateStatus;

/*
 * Verification result for one Stage 6 candidate.
 *
 * start:         proposed zero-based reference start from Stage 6.
 * votes:         Stage 5 vote count carried through Stage 6.
 * edit_distance: finite distance, or VERIFY_DISTANCE_INF when no valid
 *                banded comparison was possible.
 * accepted:      1 iff a finite distance is <= max_edit_distance.
 * status:        records whether DP was performed or why it was skipped.
 */
typedef struct {
    long start;
    size_t votes;
    size_t edit_distance;
    int accepted;
    VerifyCandidateStatus status;
} VerifyCandidateResult;

/*
 * Final result for one read.
 *
 * mapped == 1 means a candidate passed the edit-distance threshold.
 * If mapped == 0, no candidate was accepted and the read must not be placed.
 * candidates[] preserves the input candidate order (normally best-first from
 * top3_get_sorted), so the caller can display each candidate's verification.
 */
typedef struct {
    int mapped;
    long start;
    size_t votes;
    size_t edit_distance;
    size_t candidate_count;
    VerifyCandidateResult candidates[TOP3_CAPACITY];
} VerifyResult;

/*
 * Verify up to TOP3_CAPACITY candidate starts for one read.
 *
 * Inputs:
 *   reference       Packed reference sequence from Stage 1.
 *   read            Packed sample read from Stage 1.
 *   candidates      Stage 6 candidates, normally returned by
 *                   top3_get_sorted() in best-first order.
 *   candidate_count Number of entries in candidates (0..TOP3_CAPACITY).
 *   band            Calculate DP cells only where |i-j| <= band.
 *   max_edit_distance Candidate accepted when finite distance <= this value.
 *   out             Receives per-candidate details and the selected placement.
 *
 * Reference window policy (matches the project brief): for a candidate start,
 * compare the read with reference[start .. start+read_length). If fewer than
 * read_length bases remain before the reference ends, use the shorter suffix.
 * If the length difference exceeds band, distance is infinite.
 *
 * DP policy: Levenshtein distance with substitution, insertion and deletion
 * cost 1; exact match cost 0. Only cells inside the diagonal band are
 * calculated. The implementation uses a flattened 2D matrix for clarity. It
 * stores the full matrix but leaves out-of-band cells at infinity; this saves
 * DP computation time, not matrix allocation. See README for complexity.
 *
 * Unknown-base policy / explicit assumption: Stage 1 records every non-ACGT
 * character in Seq.runs and stores a placeholder A in the packed cells. To
 * avoid pretending that placeholder A is a real base, this function marks a
 * candidate unusable if the read or its reference comparison window overlaps
 * any NRun. Ambiguous-base alignment is outside the current project scope.
 *
 * Selection rule: among accepted candidates choose (1) smallest edit distance,
 * (2) more votes, then (3) smaller reference start. Rule 3 is a deterministic
 * project choice because the brief does not specify a final tie-break after
 * equal distance and equal votes.
 *
 * Returns 0 on success (including when no candidate is accepted), -1 for
 * invalid arguments, inconsistent Seq metadata, or allocation failure.
 * `out` is reset before results are populated.
 */
int verify_candidates(const Seq *reference,
                      const Seq *read,
                      const TopCandidate *candidates,
                      size_t candidate_count,
                      size_t band,
                      size_t max_edit_distance,
                      VerifyResult *out);

#endif /* VERIFY_H */
