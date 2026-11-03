# Stage 7 — Read verification using banded edit distance (C)

## 1. Folder placement

Copy `src/stage7_verify/` into the repository:

```text
src/
  stage1_load/
    seq.c
    seq.h
  stage2_minimizers/
    minimizer.h
  stage3_index/
    index.h
  stage4_bloom/
    bloom.h
  stage5_vote/
    vote.h
  stage6_top3/
    top3.h
  stage7_verify/
    verify.c
    verify.h
    test_verify.c
```

No previous-stage source file needs to be modified for this module.

## 2. Purpose

Stage 5 produces a vote count for each distinct candidate read start. Stage 6 retains up to three candidates in a min-heap. Stage 7 checks the complete read against the reference sequence at each candidate start. This prevents a read from being placed using a few coincidental minimizer matches alone.

Stage 7 calculates Levenshtein edit distance, where:

- matching bases cost 0;
- substitution costs 1;
- insertion costs 1;
- deletion costs 1.

The toy defaults from the project specification are:

```c
#define VERIFY_DEFAULT_BAND ((size_t)2)
#define VERIFY_DEFAULT_MAX_EDIT_DISTANCE ((size_t)3)
```

These are experiment settings, not universal constants. Keep them configurable for later benchmarking.

## 3. How the banded DP works

Let `read` be the complete sample read, of length `m`. For a candidate start `s`, the reference window is:

```text
reference[s .. s + m)
```

If fewer than `m` bases remain before the end of the reference, the window is clipped to the available suffix. Its length is then `n <= m`.

The DP grid has `(m + 1) * (n + 1)` cells. The recurrence for cell `(i,j)` is:

```text
D[i][j] = min(
    D[i-1][j]   + 1,                         // deletion
    D[i][j-1]   + 1,                         // insertion
    D[i-1][j-1] + (read[i-1] == ref[j-1] ? 0 : 1)
)
```

The initial conditions are:

```text
D[0][0] = 0
D[0][j] = j
D[i][0] = i
```

Only cells satisfying `abs(i - j) <= band` are calculated. Other cells remain infinity. If the length difference alone exceeds the band, the final distance is infinity without allocating the DP matrix.

**Implementation detail:** the current version uses a flattened allocation for the complete 2D matrix because it is straightforward to explain and inspect. It calculates only in-band cells, so its DP work is approximately `O(m * band)`, but its matrix memory is `O(m * n)`. If memory becomes a measured bottleneck, a later optimization can use rolling rows or compact band storage. Do not describe the current version as `O(m * band)` memory.

## 4. Candidate selection

For each candidate supplied by Stage 6:

1. Reject a negative start.
2. Reject a start at or beyond the reference length.
3. Extract up to `read length` bases from the reference, clipping at the right edge.
4. Calculate banded edit distance.
5. Accept that candidate only if the distance is finite and `distance <= max_edit_distance`.

After examining every candidate, choose the accepted candidate using:

1. Smallest edit distance.
2. If distances tie, highest Stage 5 vote count.
3. If both tie, smallest reference start (a deterministic project choice; the brief does not specify this final tie-break).

If none passes, `mapped` remains 0. The read must not be sent to Stage 8 as a placed read.

## 5. Unknown-base policy

Stage 1 packs every non-ACGT symbol as placeholder A and records its original position in `Seq.runs`. The placeholder must not be mistaken for a biological A. Therefore, this implementation marks a candidate unusable when the read or its compared reference window overlaps an NRun.

This is a conservative policy. Ambiguous-base alignment is outside the current project scope; revisit it explicitly if that scope changes.

## 6. API

```c
int verify_candidates(const Seq *reference,
                      const Seq *read,
                      const TopCandidate candidates[TOP3_CAPACITY],
                      size_t candidate_count,
                      size_t band,
                      size_t max_edit_distance,
                      VerifyResult *out);
```

The candidate array should normally be obtained from Stage 6:

```c
TopCandidate shortlist[TOP3_CAPACITY];
size_t shortlist_count = 0;
VerifyResult verification;

if (top3_get_sorted(&top3, shortlist, &shortlist_count) != 0) {
    /* Handle Stage 6 output error. */
}

if (verify_candidates(&reference, &read,
                      shortlist, shortlist_count,
                      VERIFY_DEFAULT_BAND,
                      VERIFY_DEFAULT_MAX_EDIT_DISTANCE,
                      &verification) != 0) {
    /* Handle verification error. */
}

if (verification.mapped) {
    /* Read is accepted at verification.start. */
} else {
    /* No candidate passed; leave this read unmapped. */
}
```

`VerifyResult.candidates[i]` retains the distance and status for every candidate checked, which is useful for debugging and the eventual frontend.

### Candidate status values

- `VERIFY_CANDIDATE_EVALUATED`: DP produced a finite distance. `accepted` separately says whether it passed the threshold.
- `VERIFY_CANDIDATE_NEGATIVE_START`: candidate start was negative.
- `VERIFY_CANDIDATE_START_OUT_OF_RANGE`: candidate starts at or after the end of the reference.
- `VERIFY_CANDIDATE_UNKNOWN_BASE`: read or candidate reference window contains a recorded unknown base.
- `VERIFY_CANDIDATE_OUTSIDE_BAND`: no DP path can reach the final cell within the configured band; distance is `VERIFY_DISTANCE_INF`.

## 7. Tests

`test_verify.c` covers:

- exact matching;
- choosing the lowest edit distance even if another candidate has more votes;
- one substitution;
- one extra base (insertion/deletion cost);
- rejection above the configured threshold;
- sequence-length difference larger than the band;
- clipping a reference window at the right edge;
- refusing to treat an N placeholder as a real A;
- negative and out-of-range candidate starts;
- no candidates;
- vote-count and final coordinate tie-breaking;
- invalid arguments.

The tests use the actual Stage 1 FASTA loader. Therefore, `base_to_code()` and `seq_get()` in `stage1_load/seq.c` must already be implemented (the TODO versions that always return -1 and 0 are not sufficient).

## 8. Compile and run

From the repository root, using GCC or Clang:

```sh
cc -std=c11 -Wall -Wextra -Wpedantic \
  -Isrc/stage1_load -Isrc/stage2_minimizers \
  -Isrc/stage3_index -Isrc/stage4_bloom \
  -Isrc/stage5_vote -Isrc/stage6_top3 -Isrc/stage7_verify \
  src/stage1_load/seq.c \
  src/stage7_verify/verify.c \
  src/stage7_verify/test_verify.c \
  -o test_verify

./test_verify
```

The test uses the `TopCandidate` type declared by `top3.h`, but does not need to link `top3.c`: Stage 6 is responsible for obtaining the sorted shortlist before it calls Stage 7.

## 9. Important scope notes

- The stage verifies read placement; it does not call mutations. Stage 8 does that by combining bases from all accepted reads.
- Edit distance permits substitutions, insertions and deletions during comparison. However, the planned Stage 8 pile-up maps read offset `i` to reference position `start + i`; that pile-up does not correctly interpret indel-shifted alignments. The whole pipeline therefore remains substitution-focused unless Stage 8 is extended.
- The current read/reference orientation is forward-only, matching the project scope. Reverse-complement support is not implemented.
- `band` limits alignment shifts. A correct alignment requiring a path outside this band will be reported as unreachable even if its unrestricted edit distance is small.
- The default threshold of 3 and band of 2 are toy settings. They should be benchmarked against read length, sequencing error rate and expected true differences.
