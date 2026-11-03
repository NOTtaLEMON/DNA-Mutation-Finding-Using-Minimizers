#include "verify.h"

#include <stdlib.h>
#include <string.h>

/* A private infinity value that leaves room for safe +1 operations. */
#define DP_INF (SIZE_MAX / 4u)

/*
 * Check that the packed sequence and its N-run metadata are internally safe.
 * Stage 1's Seq contract says runs are sorted, non-overlapping and non-touching.
 */
static int sequence_is_valid(const Seq *seq)
{
    size_t required_bytes;
    size_t i;
    size_t previous_end = 0;

    if (seq == NULL || seq->len == 0 || seq->cells == NULL)
        return 0;

    /* Compute ceil(len / 4) without adding 3 (which could overflow). */
    required_bytes = seq->len / 4 + (seq->len % 4 != 0);
    if (seq->cells_cap < required_bytes)
        return 0;

    if (seq->nruns > 0 && seq->runs == NULL)
        return 0;

    for (i = 0; i < seq->nruns; i++) {
        size_t start = seq->runs[i].start;
        size_t length = seq->runs[i].length;
        size_t end;

        if (length == 0 || start >= seq->len || length > seq->len - start)
            return 0;

        end = start + length; /* Safe because length <= seq->len - start. */
        if (i > 0 && start <= previous_end)
            return 0; /* Stage 1 promises sorted runs that neither overlap nor touch. */
        previous_end = end;
    }

    return 1;
}

/*
 * Return 1 if [start, start+length) overlaps an NRun, 0 if it does not,
 * and -1 if the run metadata is invalid. The range must already be in bounds.
 */
static int range_contains_unknown(const Seq *seq, size_t start, size_t length)
{
    size_t end;
    size_t i;

    if (seq == NULL || start > seq->len || length > seq->len - start)
        return -1;

    end = start + length; /* Safe after the bounds check above. */

    for (i = 0; i < seq->nruns; i++) {
        size_t run_start = seq->runs[i].start;
        size_t run_end = run_start + seq->runs[i].length;

        if (run_start >= end)
            break; /* Runs are sorted, so no later run can overlap. */
        if (run_end > start)
            return 1;
    }

    return 0;
}

/* Add one to a DP value without allowing the infinity sentinel to overflow. */
static size_t dp_plus_one(size_t value)
{
    return value >= DP_INF ? DP_INF : value + 1;
}

/* Return the smaller of three candidate DP costs. */
static size_t minimum3(size_t a, size_t b, size_t c)
{
    size_t result = a < b ? a : b;
    return result < c ? result : c;
}

/*
 * Calculate Levenshtein edit distance between:
 *   - the complete read [0, read->len), and
 *   - reference[reference_start .. reference_start+reference_length).
 *
 * This function assumes both ranges are valid and contain only known A/C/G/T
 * bases. It computes only |i-j| <= band, leaving all other cells at DP_INF.
 *
 * Returns 0 on success and writes either a finite distance or
 * VERIFY_DISTANCE_INF. Returns -1 if matrix-size arithmetic overflows or
 * allocation fails.
 */
static int banded_distance(const Seq *read,
                           const Seq *reference,
                           size_t reference_start,
                           size_t reference_length,
                           size_t band,
                           size_t *distance_out)
{
    const size_t read_length = read->len;
    const size_t rows = read_length + 1;
    const size_t columns = reference_length + 1;
    size_t cell_count;
    size_t *dp;
    size_t i, j;
    size_t initial_limit;
    size_t final_distance;

    if (distance_out == NULL)
        return -1;
    *distance_out = VERIFY_DISTANCE_INF;

    /* A path cannot reach the bottom-right cell if the length difference
       already exceeds the allowed diagonal band. */
    if (read_length > reference_length) {
        if (read_length - reference_length > band)
            return 0;
    } else if (reference_length - read_length > band) {
        return 0;
    }

    /* Guard every dimension calculation before allocating the matrix. */
    if (rows == 0 || columns == 0 || rows > SIZE_MAX / columns)
        return -1;
    cell_count = rows * columns;
    if (cell_count > SIZE_MAX / sizeof *dp)
        return -1;

    dp = (size_t *)malloc(cell_count * sizeof *dp);
    if (dp == NULL)
        return -1;

    /* Initially mark every cell unreachable. Only the band is filled below. */
    for (i = 0; i < cell_count; i++)
        dp[i] = DP_INF;

    /* Base case: empty read prefix against reference prefixes. */
    initial_limit = reference_length < band ? reference_length : band;
    for (j = 0; j <= initial_limit; j++)
        dp[j] = j;

    /* Base case: read prefixes against an empty reference prefix. */
    initial_limit = read_length < band ? read_length : band;
    for (i = 0; i <= initial_limit; i++)
        dp[i * columns] = i;

    /*
     * Fill only the diagonal band. For row i:
     *   first column = max(1, i-band)
     *   last column  = min(reference_length, i+band)
     * The calculations below avoid unsigned underflow/overflow.
     */
    for (i = 1; i <= read_length; i++) {
        size_t first_j = (i > band) ? (i - band) : 1;
        size_t last_j;

        if (i >= reference_length || band >= reference_length - i)
            last_j = reference_length;
        else
            last_j = i + band;

        if (first_j > last_j)
            continue;

        for (j = first_j; j <= last_j; j++) {
            size_t above = dp[(i - 1) * columns + j];
            size_t left = dp[i * columns + (j - 1)];
            size_t diagonal = dp[(i - 1) * columns + (j - 1)];
            uint8_t read_code = seq_get(read, i - 1);
            uint8_t reference_code = seq_get(reference, reference_start + j - 1);
            size_t substitution_cost = (read_code == reference_code) ? 0u : 1u;

            dp[i * columns + j] = minimum3(
                dp_plus_one(above),
                dp_plus_one(left),
                (diagonal >= DP_INF) ? DP_INF : diagonal + substitution_cost
            );
        }
    }

    final_distance = dp[read_length * columns + reference_length];
    free(dp);

    if (final_distance >= DP_INF)
        *distance_out = VERIFY_DISTANCE_INF;
    else
        *distance_out = final_distance;

    return 0;
}

/*
 * Is candidate a better accepted placement than the current best?
 * Primary rule comes from the project brief: smaller edit distance wins;
 * ties go to more votes. Smaller start is a deterministic final tie-break.
 */
static int candidate_is_better(const VerifyCandidateResult *candidate,
                               const VerifyResult *current_best)
{
    if (!current_best->mapped)
        return 1;
    if (candidate->edit_distance != current_best->edit_distance)
        return candidate->edit_distance < current_best->edit_distance;
    if (candidate->votes != current_best->votes)
        return candidate->votes > current_best->votes;
    return candidate->start < current_best->start;
}

int verify_candidates(const Seq *reference,
                      const Seq *read,
                      const TopCandidate *candidates,
                      size_t candidate_count,
                      size_t band,
                      size_t max_edit_distance,
                      VerifyResult *out)
{
    size_t i;
    int read_has_unknown;

    if (out == NULL)
        return -1;

    memset(out, 0, sizeof *out);
    out->edit_distance = VERIFY_DISTANCE_INF;

    if (!sequence_is_valid(reference) || !sequence_is_valid(read))
        return -1;
    if (candidate_count > TOP3_CAPACITY)
        return -1;
    if (candidate_count > 0 && candidates == NULL)
        return -1;

    read_has_unknown = range_contains_unknown(read, 0, read->len);
    if (read_has_unknown < 0)
        return -1;

    out->candidate_count = candidate_count;

    for (i = 0; i < candidate_count; i++) {
        VerifyCandidateResult *result = &out->candidates[i];
        size_t start;
        size_t available;
        size_t reference_window_length;
        int reference_has_unknown;

        result->start = candidates[i].start;
        result->votes = candidates[i].votes;
        result->edit_distance = VERIFY_DISTANCE_INF;
        result->accepted = 0;
        result->status = VERIFY_CANDIDATE_EVALUATED;

        if (candidates[i].votes == 0)
            return -1; /* Stage 6 candidates must represent at least one vote. */

        /* Negative candidate starts are legal in Stage 5/6 but cannot map. */
        if (candidates[i].start < 0) {
            result->status = VERIFY_CANDIDATE_NEGATIVE_START;
            continue;
        }

        /* Convert safely after checking the signed value is non-negative. */
        if ((uintmax_t)candidates[i].start > (uintmax_t)SIZE_MAX) {
            result->status = VERIFY_CANDIDATE_START_OUT_OF_RANGE;
            continue;
        }
        start = (size_t)candidates[i].start;

        if (start >= reference->len) {
            result->status = VERIFY_CANDIDATE_START_OUT_OF_RANGE;
            continue;
        }

        /* Stage 7 compares a read-length reference window, clipped at the end. */
        available = reference->len - start;
        reference_window_length = read->len < available ? read->len : available;

        /* Unknown bases are placeholder As in the packed array; never treat
           those placeholders as genuine DNA during verification. */
        if (read_has_unknown) {
            result->status = VERIFY_CANDIDATE_UNKNOWN_BASE;
            continue;
        }
        reference_has_unknown = range_contains_unknown(
            reference, start, reference_window_length
        );
        if (reference_has_unknown < 0)
            return -1;
        if (reference_has_unknown) {
            result->status = VERIFY_CANDIDATE_UNKNOWN_BASE;
            continue;
        }

        if (banded_distance(read, reference, start,
                            reference_window_length, band,
                            &result->edit_distance) != 0)
            return -1;

        if (result->edit_distance == VERIFY_DISTANCE_INF) {
            result->status = VERIFY_CANDIDATE_OUTSIDE_BAND;
            continue;
        }

        result->status = VERIFY_CANDIDATE_EVALUATED;
        if (result->edit_distance <= max_edit_distance) {
            result->accepted = 1;
            if (candidate_is_better(result, out)) {
                out->mapped = 1;
                out->start = result->start;
                out->votes = result->votes;
                out->edit_distance = result->edit_distance;
            }
        }
    }

    return 0;
}
