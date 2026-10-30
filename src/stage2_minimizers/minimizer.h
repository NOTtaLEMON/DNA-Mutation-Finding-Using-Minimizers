
#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <stddef.h>
#include <stdint.h>

#include "seq.h"

/*
 * One selected minimizer occurrence.
 *
 * hash: stable 64-bit hash of the selected k-mer. This is an ordering key,
 *       NOT proof that two k-mers are identical (hash collisions are possible).
 * pos:  zero-based start position of this k-mer in the source Seq.
 */
typedef struct {
    uint64_t hash;
    size_t pos;
} Minimizer;

/* Dynamic array of selected minimizer occurrences. */
typedef struct {
    Minimizer *items;
    size_t count;
    size_t capacity;
} MinimizerSet;

/*
 * Generate minimizers from a packed sequence.
 *
 * PROJECT CONVENTION / ASSUMPTION:
 * This module implements the common read-mapping scheme in which:
 *   1. Every consecutive k-base substring is a k-mer.
 *   2. Each k-mer is ranked by a deterministic hash.
 *   3. Each window contains `window_kmers` consecutive k-mers.
 *   4. The lowest-hash k-mer in each full window is selected.
 *
 * This is NOT the variable-length finimizer scheme from the 2025 paper, and
 * it is not the paper's alternate formalism where an m-mer is selected inside
 * a k-mer. We choose this convention because it matches this project's
 * candidate-position calculation (reference k-mer start - read k-mer start).
 * Revisit this choice if the project specification defines k/m/windows
 * differently.
 *
 * N/ambiguous-base policy:
 * Stage 1 stores every non-ACGT character as a placeholder A in `cells` and
 * records its original location in `Seq.runs`. This function checks those
 * runs and NEVER forms a k-mer that overlaps an unknown position. Thus the
 * placeholder A cannot create a false minimizer.
 *
 * Window boundary policy:
 * Windows are counted in k-mers, not DNA bases. A full window has exactly
 * `window_kmers` valid consecutive k-mers. Unknown bases break the valid
 * k-mer stream; a window never crosses such a break. If a valid segment has
 * fewer than `window_kmers` k-mers, this implementation emits its single
 * lowest-hash k-mer. This short-segment rule is an explicit project choice
 * (some minimizer implementations instead emit no minimizer for that case).
 *
 * Repeated-selection policy:
 * If the same k-mer occurrence is the minimum in adjacent windows, it is
 * emitted only once. Equal hash values are resolved in favor of the leftmost
 * occurrence in the window.
 *
 * Strand policy / limitation:
 * Hashes are computed in the input's forward orientation. Reverse-complement
 * canonicalization and strand tracking are NOT implemented in Step 2. If the
 * input reads can be in either orientation, a later mapping stage must add
 * reverse-strand handling or this project must extend this module.
 *
 * Parameters:
 *   seq           Loaded sequence from stage1_load/seq.h.
 *   k             k-mer length in bases; must be > 0.
 *   window_kmers  Number of consecutive k-mers per minimizer window; > 0.
 *   out           Must point to a zero-initialized or previously-freed set.
 *
 * Returns 0 on success, -1 for invalid arguments or allocation failure.
 * If k > seq->len, success is returned with an empty output set.
 */
int minimizers_generate(const Seq *seq,
                        size_t k,
                        size_t window_kmers,
                        MinimizerSet *out);

/* Frees a set and resets it to empty. Safe for a zero-initialized set. */
void minimizers_free(MinimizerSet *set);

#endif /* MINIMIZER_H */
