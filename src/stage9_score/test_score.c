#include "score.h"

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

/* Appends one call to a CallList backed by a fixed test array. */
static void add_call(CallList *calls, size_t pos, uint8_t ref, uint8_t alt, uint32_t depth)
{
    MutationCall *c = &calls->items[calls->count++];
    c->pos = pos;
    c->ref = ref;
    c->alt = alt;
    c->alt_count = depth;
    c->depth = depth;
}

/* Unsorted adds come back sorted; binary search finds every entry and misses others. */
static void test_sort_and_find(void)
{
    AnswerKey key = {0};
    static const size_t pos[] = {500, 20, 9000, 37, 1};
    size_t i;

    for (i = 0; i < 5; i++)
        CHECK(key_add(&key, pos[i], 0, 1) == 0);
    CHECK(key_sort(&key) == 0);
    for (i = 1; i < key.count; i++)
        CHECK(key.items[i - 1].pos < key.items[i].pos);
    for (i = 0; i < 5; i++)
        CHECK(key_find(&key, pos[i]) != NULL && key_find(&key, pos[i])->pos == pos[i]);
    CHECK(key_find(&key, 0) == NULL);
    CHECK(key_find(&key, 21) == NULL);
    CHECK(key_find(&key, 9001) == NULL);

    CHECK(key_add(&key, 37, 2, 3) == 0);         /* duplicate position */
    CHECK(key_sort(&key) == -1);
    CHECK(key_add(&key, 5, 2, 2) == -1);         /* alt equal to ref is not a mutation */
    key_free(&key);
}

/* Toy scoring: 2 right, 1 wrong alt, 1 at an unplanted position, 1 planted position missed. */
static void test_score_counts(void)
{
    AnswerKey key = {0};
    MutationCall items[8];
    CallList calls = {items, 0, 8};
    uint32_t depth[100];
    ScoreResult s;

    memset(depth, 0, sizeof depth);
    depth[20] = depth[37] = depth[50] = depth[60] = 10;   /* position 80 has no reads */

    CHECK(key_add(&key, 20, 0, 3) == 0);
    CHECK(key_add(&key, 37, 1, 2) == 0);
    CHECK(key_add(&key, 50, 2, 0) == 0);
    CHECK(key_add(&key, 80, 3, 1) == 0);
    CHECK(key_sort(&key) == 0);

    add_call(&calls, 20, 0, 3, 10);              /* TP */
    add_call(&calls, 37, 1, 2, 10);              /* TP */
    add_call(&calls, 50, 2, 1, 10);              /* right place, wrong letter: FP */
    add_call(&calls, 60, 0, 1, 10);              /* not planted: FP */

    CHECK(score_calls(&calls, &key, depth, 100, 3, &s) == 0);
    CHECK(s.calls == 4 && s.planted == 4);
    CHECK(s.tp == 2 && s.fp == 2 && s.fn == 2 && s.wrong_alt == 1);
    CHECK(s.callable == 3);                      /* 80 was never covered */
    CHECK(s.precision_pct == 50.0 && s.recall_all_pct == 50.0);
    CHECK(s.recall_callable_pct > 66.6 && s.recall_callable_pct < 66.7);

    CHECK(score_calls(&calls, &key, NULL, 0, 3, &s) == 0);
    CHECK(s.callable == 0 && s.recall_callable_pct == 0.0);
    key_free(&key);
}

/* A saved key loads back identical (1-based on disk, 0-based in memory); bad files are refused. */
static void test_save_load(void)
{
    AnswerKey key = {0}, back = {0};
    FILE *f;
    size_t i;

    CHECK(key_add(&key, 1069, 0, 1) == 0);
    CHECK(key_add(&key, 0, 2, 3) == 0);
    CHECK(key_sort(&key) == 0);
    CHECK(key_save("test_score_key.txt", &key) == 0);
    CHECK(key_load("test_score_key.txt", &back) == 0);
    CHECK(back.count == 2);
    for (i = 0; i < back.count && i < key.count; i++)
        CHECK(back.items[i].pos == key.items[i].pos && back.items[i].ref == key.items[i].ref &&
              back.items[i].alt == key.items[i].alt);

    f = fopen("test_score_key.txt", "r");
    if (f != NULL) {
        char line[64];
        CHECK(fgets(line, sizeof line, f) != NULL && line[0] == '#');
        CHECK(fgets(line, sizeof line, f) != NULL && strcmp(line, "1 G T\n") == 0);
        fclose(f);
    }

    f = fopen("test_score_key.txt", "w");
    if (f != NULL) {
        fputs("# bad\n0 A C\n", f);              /* position 0 does not exist in 1-based files */
        fclose(f);
    }
    CHECK(key_load("test_score_key.txt", &back) == -1);
    remove("test_score_key.txt");

    key_free(&key);
    key_free(&back);
}

/* Runs all checks. */
int main(void)
{
    test_sort_and_find();
    test_score_counts();
    test_save_load();

    if (failures == 0)
        printf("all checks passed\n");
    else
        printf("%d check(s) failed\n", failures);
    return failures != 0;
}
