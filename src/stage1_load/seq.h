#ifndef SEQ_H
#define SEQ_H

#include <stddef.h>
#include <stdint.h>

/* A run of consecutive non-ACGT letters (normally N): positions start .. start+length-1. */
typedef struct {
    size_t start;
    size_t length;
} NRun;

/*
 * A DNA sequence stored 2 bits per letter, 4 letters per byte.
 * Position p lives in byte p / 4. The first letter of each byte is in the top two bits.
 * Every position has exactly one cell, so the index of a cell is the position of its letter.
 * Under an N the cell holds the placeholder 00; the N list is the only record that it was an N.
 */
typedef struct {
    uint8_t *cells;      /* packed letters */
    size_t cells_cap;    /* capacity of cells, in bytes */
    size_t len;          /* number of letters (positions) stored */
    NRun *runs;          /* sorted by start, never overlapping or touching */
    size_t nruns;
    size_t runs_cap;
} Seq;

/* Reads the first record of a FASTA file. Returns 0 on success, -1 on failure. */
int seq_load_fasta(const char *path, Seq *s);

void seq_free(Seq *s);

/* YOUR PART (defined in seq.c) */

/* 'A','C','G','T' (either case) -> 0,1,2,3. Anything else -> -1. */
int base_to_code(char c);

/* The 2-bit code (0..3) stored at position pos. The caller guarantees pos < s->len. */
uint8_t seq_get(const Seq *s, size_t pos);

#endif
