#include "seq.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 * YOUR PART. Two small functions. The contracts are in seq.h.
 * ------------------------------------------------------------------------- */

/* Converts one letter to its 2-bit code: A=0, C=1, G=2, T=3 (either case); -1 for anything else. */
int base_to_code(char c)
{
    (void)c;
    /* TODO */
    return -1;
}

/* Returns the 2-bit code stored at position pos (byte pos / 4, top two bits first). */
uint8_t seq_get(const Seq *s, size_t pos)
{
    (void)s;
    (void)pos;
    /* TODO */
    return 0;
}

/* ---------------------------------------------------------------------------
 * Loader
 * ------------------------------------------------------------------------- */

/* Appends one 2-bit code at position s->len. Returns 0 on success, -1 if memory ran out. */
static int push_code(Seq *s, uint8_t code)
{
    size_t byte = s->len / 4;

    if (byte >= s->cells_cap) {
        size_t new_cap = s->cells_cap ? s->cells_cap * 2 : 1024;
        uint8_t *p = realloc(s->cells, new_cap);
        if (p == NULL)
            return -1;
        memset(p + s->cells_cap, 0, new_cap - s->cells_cap);   /* new bytes start as 00 */
        s->cells = p;
        s->cells_cap = new_cap;
    }

    /* cells are zero to begin with, so OR-ing the code into its slot is enough */
    s->cells[byte] |= (uint8_t)(code << (6 - 2 * (s->len % 4)));
    s->len++;
    return 0;
}

/* Records that position pos is an N. Extends the last run if pos touches it. */
static int note_n(Seq *s, size_t pos)
{
    if (s->nruns > 0) {
        NRun *last = &s->runs[s->nruns - 1];
        if (last->start + last->length == pos) {
            last->length++;
            return 0;
        }
    }

    if (s->nruns == s->runs_cap) {
        size_t new_cap = s->runs_cap ? s->runs_cap * 2 : 16;
        NRun *p = realloc(s->runs, new_cap * sizeof *p);
        if (p == NULL)
            return -1;
        s->runs = p;
        s->runs_cap = new_cap;
    }

    s->runs[s->nruns].start = pos;
    s->runs[s->nruns].length = 1;
    s->nruns++;
    return 0;
}

/*
 * Reads the first record of the FASTA file at path into s in one pass: the header line and line
 * breaks are skipped, A/C/G/T go straight into the 2-bit array, and every other letter gets the
 * placeholder 00 plus an entry in the N list. Returns 0 on success, -1 if the file can't be
 * opened or memory runs out (s is left empty and freed in that case).
 */
int seq_load_fasta(const char *path, Seq *s)
{
    FILE *f;
    int ch;
    int in_header = 0;
    int seen_record = 0;

    memset(s, 0, sizeof *s);

    f = fopen(path, "rb");
    if (f == NULL)
        return -1;

    while ((ch = fgetc(f)) != EOF) {
        if (ch == '>') {
            if (seen_record)          /* a second record: stop, only the first is loaded */
                break;
            seen_record = 1;
            in_header = 1;
            continue;
        }
        if (in_header) {
            if (ch == '\n')
                in_header = 0;
            continue;
        }
        if (ch == '\r' || ch == '\n' || ch == ' ' || ch == '\t')
            continue;

        {
            int code = base_to_code((char)ch);
            size_t pos = s->len;

            if (code >= 0) {
                if (push_code(s, (uint8_t)code) != 0)
                    goto fail;
            } else {
                /* N or any other non-ACGT letter: placeholder cell plus an N-list entry */
                if (push_code(s, 0) != 0 || note_n(s, pos) != 0)
                    goto fail;
            }
        }
    }

    fclose(f);
    return 0;

fail:
    fclose(f);
    seq_free(s);
    return -1;
}

/* Frees the packed array and the N list and resets s to empty. Safe on an already-empty Seq. */
void seq_free(Seq *s)
{
    free(s->cells);
    free(s->runs);
    memset(s, 0, sizeof *s);
}
