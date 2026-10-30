#include "index.h"

#include <stdlib.h>
#include <string.h>

/* Packs k letters starting at pos into one integer, 2 bits per letter, first letter highest. */
uint64_t word_id_at(const Seq *seq, size_t pos, size_t k)
{
    uint64_t id = 0;
    size_t i;

    for (i = 0; i < k; i++)
        id = (id << 2) | (uint64_t)seq_get(seq, pos + i);

    return id;
}

/* Appends pos at the end of the entry's position list. Returns 0 on success, -1 if memory ran out. */
static int add_position(IndexEntry *e, size_t pos)
{
    PosNode *node = (PosNode *)malloc(sizeof *node);
    if (node == NULL)
        return -1;

    node->pos = pos;
    node->next = NULL;

    if (e->pos_tail == NULL)
        e->pos_head = node;
    else
        e->pos_tail->next = node;
    e->pos_tail = node;
    e->npos++;
    return 0;
}

/*
 * Inserts one kept (hash, pos) into the table.
 * Walks the bucket's chain: if an entry with the same hash AND the same word ID exists, the position
 * is appended to it; otherwise a new entry is added at the end of the chain.
 * Returns 0 on success, -1 if memory ran out.
 */
static int insert_kept_word(Index *idx, const Seq *seq, uint64_t hash, size_t pos)
{
    uint64_t id = word_id_at(seq, pos, idx->k);
    size_t b = (size_t)(hash % idx->nbuckets);
    IndexEntry *e = idx->buckets[b];
    IndexEntry *last = NULL;
    IndexEntry *fresh;

    for (; e != NULL; e = e->next) {
        if (e->hash == hash && e->word_id == id)   /* hash first (cheap filter), then the real check */
            return add_position(e, pos);
        last = e;
    }

    fresh = (IndexEntry *)calloc(1, sizeof *fresh);
    if (fresh == NULL)
        return -1;
    fresh->word_id = id;
    fresh->hash = hash;
    if (add_position(fresh, pos) != 0) {
        free(fresh);
        return -1;
    }

    if (last == NULL)
        idx->buckets[b] = fresh;
    else
        last->next = fresh;
    idx->nentries++;
    return 0;
}

/* Builds the index from the reference's minimizers. See index.h for the full contract. */
int index_build(const Seq *seq, const MinimizerSet *mins, size_t k, size_t nbuckets, Index *out)
{
    Index idx;
    size_t i;

    memset(out, 0, sizeof *out);

    if (seq == NULL || mins == NULL || k == 0 || k > 32)
        return -1;
    if (mins->count > 0 && mins->items == NULL)
        return -1;

    /* every kept word must lie fully inside the sequence, or word_id_at would read past the end */
    for (i = 0; i < mins->count; i++)
        if (k > seq->len || mins->items[i].pos > seq->len - k)
            return -1;

    /* size once, before building: about count / 0.75 buckets, at least 1 */
    if (nbuckets == 0)
        nbuckets = mins->count + mins->count / 3 + 1;

    memset(&idx, 0, sizeof idx);
    idx.k = k;
    idx.nbuckets = nbuckets;
    idx.buckets = (IndexEntry **)calloc(nbuckets, sizeof *idx.buckets);
    if (idx.buckets == NULL)
        return -1;

    for (i = 0; i < mins->count; i++) {
        if (insert_kept_word(&idx, seq, mins->items[i].hash, mins->items[i].pos) != 0) {
            index_free(&idx);
            return -1;
        }
    }

    *out = idx;
    return 0;
}

/* Finds the entry with this hash and word ID, or NULL. Counts the entries checked into *probes. */
const IndexEntry *index_lookup(const Index *idx, uint64_t hash, uint64_t word_id, size_t *probes)
{
    const IndexEntry *e;
    size_t checked = 0;

    if (idx == NULL || idx->buckets == NULL || idx->nbuckets == 0) {
        if (probes != NULL)
            *probes = 0;
        return NULL;
    }

    for (e = idx->buckets[hash % idx->nbuckets]; e != NULL; e = e->next) {
        checked++;
        if (e->hash == hash && e->word_id == word_id)
            break;
    }

    if (probes != NULL)
        *probes = checked;
    return e;
}

/* Frees all entries, position nodes and the bucket array; safe on an empty Index. */
void index_free(Index *idx)
{
    size_t b;

    if (idx == NULL)
        return;

    for (b = 0; idx->buckets != NULL && b < idx->nbuckets; b++) {
        IndexEntry *e = idx->buckets[b];
        while (e != NULL) {
            IndexEntry *next_entry = e->next;
            PosNode *p = e->pos_head;
            while (p != NULL) {
                PosNode *next_pos = p->next;
                free(p);
                p = next_pos;
            }
            free(e);
            e = next_entry;
        }
    }

    free(idx->buckets);
    memset(idx, 0, sizeof *idx);
}
