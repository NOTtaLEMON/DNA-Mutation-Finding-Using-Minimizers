/*
 * INACTIVE: indel-aware Stage 7. Compiled only with -DVERIFY_INDEL_ENABLED.
 * See verify_indel.h for what differs from verify.c.
 */
#ifdef VERIFY_INDEL_ENABLED

#include "verify_indel.h"

#include <stdlib.h>
#include <string.h>

/* A private infinity value that leaves room for safe +1 operations. */
#define DP_INF (SIZE_MAX / 4u)

/*
 * Check that the packed sequence and its N-run metadata are internally safe.
 * Same rule as verify.c (that copy is static, so it is repeated here).
 */
static int sequence_is_valid(const Seq *seq)
{
    size_t required_bytes;
    size_t i;
    size_t previous_end = 0;

    if (seq == NULL || seq->len == 0 || seq->cells == NULL)
        return 0;

    required_bytes = seq->len / 4 + (seq->len % 4 != 0);
    if (seq->cells_cap < required_bytes)
        return 0;

    if (seq->nruns > 0 && seq->runs == NULL)
        return 0;

    for (i = 0; i < seq->nruns; i++) {
        size_t start = seq->runs[i].start;
        size_t length = seq->runs[i].length;

        if (length == 0 || start >= seq->len || length > seq->len - start)
            return 0;
        if (i > 0 && start <= previous_end)
            return 0;
        previous_end = start + length;
    }

    return 1;
}

/*
 * Return 1 if [start, start+length) overlaps an NRun, 0 if it does not,
 * and -1 if the range is out of bounds.
 */
static int range_contains_unknown(const Seq *seq, size_t start, size_t length)
{
    size_t end;
    size_t i;

    if (seq == NULL || start > seq->len || length > seq->len - start)
        return -1;

    end = start + length;

    for (i = 0; i < seq->nruns; i++) {
        size_t run_start = seq->runs[i].start;
        size_t run_end = run_start + seq->runs[i].length;

        if (run_start >= end)
            break;
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

/* Distance between two column numbers, used to prefer ends near the diagonal. */
static size_t column_gap(size_t a, size_t b)
{
    return a > b ? a - b : b - a;
}

/*
 * Walk back from dp[read_length][end_column] to dp[0][0] and write the
 * alignment path into result->ops (first read letter first).
 * At each cell it prefers 'M', then 'D', then 'I' when costs are equal.
 * Returns 0 on success, -1 on allocation failure or a broken table.
 */
static int traceback(const size_t *dp,
                     size_t columns,
                     const Seq *read,
                     const Seq *reference,
                     size_t reference_start,
                     size_t end_column,
                     IndelCandidateResult *result)
{
    size_t i = read->len;
    size_t j = end_column;
    size_t capacity = read->len + end_column;
    size_t n = 0;
    size_t k;
    char *reversed = (char *)malloc(capacity);

    if (reversed == NULL)
        return -1;

    while (i > 0 || j > 0) {
        size_t current = dp[i * columns + j];

        if (i > 0 && j > 0) {
            size_t diagonal = dp[(i - 1) * columns + (j - 1)];
            size_t cost = seq_get(read, i - 1) !=
                          seq_get(reference, reference_start + j - 1);
            if (diagonal < DP_INF && diagonal + cost == current) {
                reversed[n++] = ALIGN_MATCH;
                i--;
                j--;
                continue;
            }
        }
        if (j > 0 && dp_plus_one(dp[i * columns + (j - 1)]) == current) {
            reversed[n++] = ALIGN_DELETE;
            j--;
            continue;
        }
        if (i > 0 && dp_plus_one(dp[(i - 1) * columns + j]) == current) {
            reversed[n++] = ALIGN_INSERT;
            i--;
            continue;
        }

        free(reversed);
        return -1; /* No neighbour explains this cell: the table is broken. */
    }

    result->ops = (char *)malloc(n > 0 ? n : 1);
    if (result->ops == NULL) {
        free(reversed);
        return -1;
    }
    for (k = 0; k < n; k++)
        result->ops[k] = reversed[n - 1 - k];
    result->nops = n;
    result->ref_span = end_column;

    free(reversed);
    return 0;
}

/*
 * Banded DP of the whole read against reference[start .. start+window_length).
 * Unlike verify.c, the alignment may end at any column in the last band
 * (|read_length - j| <= band), so a trailing reference letter left over by an
 * indel is not charged. The best end is the smallest cost; on a tie, the end
 * closest to the diagonal, then the smaller column.
 *
 * Writes the distance (or VERIFY_DISTANCE_INF). When want_path is non-zero
 * and the distance is finite, also fills result->ops via traceback().
 * Returns 0 on success, -1 on overflow or allocation failure.
 */
static int banded_alignment(const Seq *read,
                            const Seq *reference,
                            size_t reference_start,
                            size_t window_length,
                            size_t band,
                            size_t max_edit_distance,
                            IndelCandidateResult *result)
{
    const size_t read_length = read->len;
    const size_t rows = read_length + 1;
    const size_t columns = window_length + 1;
    size_t cell_count;
    size_t *dp;
    size_t i, j;
    size_t initial_limit;
    size_t first_end;
    size_t best_column = 0;
    size_t best_distance = DP_INF;
    int status = 0;

    result->edit_distance = VERIFY_DISTANCE_INF;

    /* Even the longest allowed end column cannot reach the read's last row. */
    if (window_length + band < read_length)
        return 0;

    if (rows > SIZE_MAX / columns)
        return -1;
    cell_count = rows * columns;
    if (cell_count > SIZE_MAX / sizeof *dp)
        return -1;

    dp = (size_t *)malloc(cell_count * sizeof *dp);
    if (dp == NULL)
        return -1;

    for (i = 0; i < cell_count; i++)
        dp[i] = DP_INF;

    /* Base cases inside the band: leading deletions and leading insertions. */
    initial_limit = window_length < band ? window_length : band;
    for (j = 0; j <= initial_limit; j++)
        dp[j] = j;
    initial_limit = read_length < band ? read_length : band;
    for (i = 0; i <= initial_limit; i++)
        dp[i * columns] = i;

    /* Fill only the diagonal band, as in verify.c. */
    for (i = 1; i <= read_length; i++) {
        size_t first_j = (i > band) ? (i - band) : 1;
        size_t last_j;

        if (i >= window_length || band >= window_length - i)
            last_j = window_length;
        else
            last_j = i + band;

        for (j = first_j; j <= last_j; j++) {
            size_t above = dp[(i - 1) * columns + j];
            size_t left = dp[i * columns + (j - 1)];
            size_t diagonal = dp[(i - 1) * columns + (j - 1)];
            size_t substitution_cost =
                seq_get(read, i - 1) != seq_get(reference, reference_start + j - 1);

            dp[i * columns + j] = minimum3(
                dp_plus_one(above),
                dp_plus_one(left),
                (diagonal >= DP_INF) ? DP_INF : diagonal + substitution_cost
            );
        }
    }

    /* Free reference end: pick the best cell of the last row inside the band. */
    first_end = read_length > band ? read_length - band : 0;
    for (j = first_end; j <= window_length; j++) {
        size_t value = dp[read_length * columns + j];

        if (value < best_distance ||
            (value == best_distance && value < DP_INF &&
             column_gap(j, read_length) < column_gap(best_column, read_length))) {
            best_distance = value;
            best_column = j;
        }
    }

    if (best_distance < DP_INF) {
        result->edit_distance = best_distance;
        /* The path is only needed when Stage 8 may use this candidate. */
        if (best_distance <= max_edit_distance)
            status = traceback(dp, columns, read, reference, reference_start,
                               best_column, result);
    }

    free(dp);
    return status;
}

/*
 * Is candidate a better accepted placement than the current best?
 * Same rule as verify.c: smaller distance, then more votes, then smaller start.
 */
static int candidate_is_better(const IndelCandidateResult *candidate,
                               const IndelVerifyResult *out)
{
    const IndelCandidateResult *best;

    if (!out->mapped)
        return 1;
    best = &out->candidates[out->best];
    if (candidate->edit_distance != best->edit_distance)
        return candidate->edit_distance < best->edit_distance;
    if (candidate->votes != best->votes)
        return candidate->votes > best->votes;
    return candidate->start < best->start;
}

/* Frees every ops array inside out and resets it. See verify_indel.h. */
void verify_indel_result_free(IndelVerifyResult *out)
{
    size_t i;

    if (out == NULL)
        return;
    for (i = 0; i < TOP3_CAPACITY; i++)
        free(out->candidates[i].ops);
    memset(out, 0, sizeof *out);
}

/* Checks each candidate with the indel-aware alignment. See verify_indel.h. */
int verify_candidates_indel(const Seq *reference,
                            const Seq *read,
                            const TopCandidate *candidates,
                            size_t candidate_count,
                            size_t band,
                            size_t max_edit_distance,
                            IndelVerifyResult *out)
{
    size_t i;
    int read_has_unknown;

    if (out == NULL)
        return -1;
    memset(out, 0, sizeof *out);

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
        IndelCandidateResult *result = &out->candidates[i];
        size_t start;
        size_t available;
        size_t wanted;
        size_t window_length;
        int reference_has_unknown;

        result->start = candidates[i].start;
        result->votes = candidates[i].votes;
        result->edit_distance = VERIFY_DISTANCE_INF;
        result->status = VERIFY_CANDIDATE_EVALUATED;

        if (candidates[i].votes == 0) {
            verify_indel_result_free(out);
            return -1;
        }
        if (candidates[i].start < 0) {
            result->status = VERIFY_CANDIDATE_NEGATIVE_START;
            continue;
        }
        if ((uintmax_t)candidates[i].start >= (uintmax_t)reference->len) {
            result->status = VERIFY_CANDIDATE_START_OUT_OF_RANGE;
            continue;
        }
        start = (size_t)candidates[i].start;

        /* Window = read length + band, clipped at the reference end, so a
           read with up to `band` missing letters still fits. */
        available = reference->len - start;
        wanted = (band > SIZE_MAX - read->len) ? SIZE_MAX : read->len + band;
        window_length = wanted < available ? wanted : available;

        if (read_has_unknown) {
            result->status = VERIFY_CANDIDATE_UNKNOWN_BASE;
            continue;
        }
        reference_has_unknown = range_contains_unknown(reference, start,
                                                       window_length);
        if (reference_has_unknown < 0) {
            verify_indel_result_free(out);
            return -1;
        }
        if (reference_has_unknown) {
            result->status = VERIFY_CANDIDATE_UNKNOWN_BASE;
            continue;
        }

        if (banded_alignment(read, reference, start, window_length, band,
                             max_edit_distance, result) != 0) {
            verify_indel_result_free(out);
            return -1;
        }

        if (result->edit_distance == VERIFY_DISTANCE_INF) {
            result->status = VERIFY_CANDIDATE_OUTSIDE_BAND;
            continue;
        }

        if (result->edit_distance <= max_edit_distance) {
            result->accepted = 1;
            if (candidate_is_better(result, out)) {
                out->mapped = 1;
                out->best = i;
            }
        }
    }

    return 0;
}

#endif /* VERIFY_INDEL_ENABLED */
