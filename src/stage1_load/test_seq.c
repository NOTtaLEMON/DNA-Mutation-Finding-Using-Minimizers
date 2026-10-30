#include "seq.h"

#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("FAIL line %d: %s\n", __LINE__, #cond);                   \
            failures++;                                                      \
        }                                                                    \
    } while (0)

static const char LETTER[] = "ACGT";

/* Writes the first n bytes of text to a new file at path (used to build small test inputs). */
static void write_file(const char *path, const char *text, size_t n)
{
    FILE *f = fopen(path, "wb");
    if (f) {
        fwrite(text, 1, n, f);
        fclose(f);
    }
}

/* CRLF line endings, a lowercase line, and an N run that crosses a line break */
static void test_demo(void)
{
    static const char text[] =
        ">demo toy reference\r\nGACATCCANN\r\nNNACGAAGCT\r\ntcgcggtg\r\n";
    /* what the cells should decode to: the four Ns (positions 8..11) show as the placeholder A */
    static const char expect[] = "GACATCCAAAAAACGAAGCTTCGCGGTG";
    Seq s;
    size_t i;

    write_file("test_demo.fasta", text, sizeof text - 1);
    CHECK(seq_load_fasta("test_demo.fasta", &s) == 0);

    CHECK(s.len == 28);
    CHECK(s.nruns == 1);
    if (s.nruns == 1) {
        CHECK(s.runs[0].start == 8);
        CHECK(s.runs[0].length == 4);
    }
    for (i = 0; i < s.len && i < sizeof expect - 1; i++)
        CHECK(LETTER[seq_get(&s, i)] == expect[i]);

    seq_free(&s);
    remove("test_demo.fasta");
}

/* only the first record is loaded; separate N runs stay separate */
static void test_two_records(void)
{
    static const char text[] = ">one\nACNGT\nNA\n>two\nTTTTTTTT\n";
    Seq s;

    write_file("test_two.fasta", text, sizeof text - 1);
    CHECK(seq_load_fasta("test_two.fasta", &s) == 0);
    CHECK(s.len == 7);
    CHECK(s.nruns == 2);
    if (s.nruns == 2) {
        CHECK(s.runs[0].start == 2 && s.runs[0].length == 1);
        CHECK(s.runs[1].start == 5 && s.runs[1].length == 1);
    }
    seq_free(&s);
    remove("test_two.fasta");
}

/* Loading a file that doesn't exist must return -1 instead of crashing. */
static void test_missing_file(void)
{
    Seq s;
    CHECK(seq_load_fasta("no_such_file.fasta", &s) == -1);
}

/* Loads a real FASTA file and prints its length, its first 10 N runs and its first 30 letters. */
static void report_real_file(const char *path)
{
    Seq s;
    size_t i;

    if (seq_load_fasta(path, &s) != 0) {
        printf("(could not open %s, skipping)\n", path);
        return;
    }
    printf("%s: %lu letters, %lu N run(s)\n", path, (unsigned long)s.len, (unsigned long)s.nruns);
    for (i = 0; i < s.nruns && i < 10; i++)
        printf("  N run: start %lu, length %lu\n",
               (unsigned long)s.runs[i].start, (unsigned long)s.runs[i].length);
    printf("  first 30 letters: ");
    for (i = 0; i < 30 && i < s.len; i++)
        putchar(LETTER[seq_get(&s, i)]);
    printf("\n");
    seq_free(&s);
}

/* Runs all checks, then reports on the FASTA file given as argv[1] (default sequence.fasta). */
int main(int argc, char **argv)
{
    test_demo();
    test_two_records();
    test_missing_file();

    if (failures == 0)
        printf("all checks passed\n");
    else
        printf("%d check(s) failed\n", failures);

    report_real_file(argc > 1 ? argv[1] : "sequence.fasta");
    return failures != 0;
}
