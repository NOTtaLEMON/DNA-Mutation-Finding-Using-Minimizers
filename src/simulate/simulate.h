#ifndef SIMULATE_H
#define SIMULATE_H

#include <stddef.h>
#include <stdint.h>

#include "seq.h"
#include "score.h"

/*
 * Test-data generator for option B (planted substitutions).
 *
 *   reference --sim_plant--> answer key --sim_make_sample--> sample genome
 *   sample --sim_reads--> reads (with sequencing errors) for stages 2 to 8
 *
 * A small xorshift random generator with a seed is used instead of rand(),
 * so the same seed gives the same data on every compiler and OS.
 */

/* Random number generator state. Seed it with sim_rng_seed before use. */
typedef struct {
    uint64_t state;
} SimRng;

/* Reads cut from the sample. true_start[i] is where read i really came from. */
typedef struct {
    Seq *reads;
    size_t *true_start;
    size_t count;
} ReadSet;

/* Sets the generator state from a seed (any value, including 0). */
void sim_rng_seed(SimRng *rng, uint64_t seed);

/* Returns the next 64-bit random number. */
uint64_t sim_rng_next(SimRng *rng);

/* Returns a random number in [0, n). n must be > 0. */
size_t sim_rng_below(SimRng *rng, size_t n);

/* Returns a random number in [0, 1). */
double sim_rng_unit(SimRng *rng);

/*
 * Picks count distinct positions in the reference (never inside an N run)
 * and a random different letter for each. The key is reset, filled and
 * sorted. Returns 0, or -1 if there are not enough usable positions.
 */
int sim_plant(const Seq *reference, size_t count, SimRng *rng, AnswerKey *key);

/* Copies the reference into sample and applies every substitution in key. Returns 0 or -1. */
int sim_make_sample(const Seq *reference, const AnswerKey *key, Seq *sample);

/*
 * Cuts coverage * sample_len / read_len reads of read_len letters at random
 * starts (reads touching an N in the sample are re-drawn). Each letter is
 * replaced by a different random letter with probability error_rate
 * (a sequencing error). Returns 0, or -1 on bad arguments or no memory.
 */
int sim_reads(const Seq *sample, size_t read_len, double coverage,
              double error_rate, SimRng *rng, ReadSet *out);

/* Frees every read and the arrays, and empties the set. */
void sim_reads_free(ReadSet *set);

#endif /* SIMULATE_H */
