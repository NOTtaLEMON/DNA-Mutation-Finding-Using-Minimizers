#include "vote.h"

#include <stdlib.h>
#include <string.h>

/* Allocates a new node for start with count 1 and no children. Returns NULL if memory ran out. */
static VoteNode *new_node(long start)
{
    VoteNode *node = (VoteNode *)malloc(sizeof *node);
    if (node == NULL)
        return NULL;
    node->start = start;
    node->count = 1;
    node->left = NULL;
    node->right = NULL;
    return node;
}

/* Adds one vote for start using a loop down one path of the tree. See vote.h. */
int vote_tree_add(VoteTree *tree, long start)
{
    VoteNode *cur;

    if (tree->root == NULL) {
        tree->root = new_node(start);
        if (tree->root == NULL)
            return -1;
        tree->nodes = 1;
        return 0;
    }

    cur = tree->root;
    for (;;) {
        if (start == cur->start) {          /* seen before: one more vote, no new node */
            cur->count++;
            return 0;
        }
        if (start < cur->start) {           /* smaller: go left, or attach here if left is empty */
            if (cur->left == NULL) {
                cur->left = new_node(start);
                if (cur->left == NULL)
                    return -1;
                tree->nodes++;
                return 0;
            }
            cur = cur->left;
        } else {                            /* larger: go right, or attach here if right is empty */
            if (cur->right == NULL) {
                cur->right = new_node(start);
                if (cur->right == NULL)
                    return -1;
                tree->nodes++;
                return 0;
            }
            cur = cur->right;
        }
    }
}

/* Recursive in-order visit: left subtree, this node, right subtree. */
static void walk_node(const VoteNode *node, VoteVisitFn visit, void *ctx)
{
    if (node == NULL)
        return;
    walk_node(node->left, visit, ctx);
    visit(node->start, node->count, ctx);
    walk_node(node->right, visit, ctx);
}

/* Calls visit for every node in increasing order of start. See vote.h. */
void vote_tree_walk(const VoteTree *tree, VoteVisitFn visit, void *ctx)
{
    if (tree == NULL || visit == NULL)
        return;
    walk_node(tree->root, visit, ctx);
}

/* Recursive post-order free: both children first, then the node, so no pointer is used after free. */
static void free_node(VoteNode *node)
{
    if (node == NULL)
        return;
    free_node(node->left);
    free_node(node->right);
    free(node);
}

/* Frees all nodes and empties the tree. See vote.h. */
void vote_tree_clear(VoteTree *tree)
{
    if (tree == NULL)
        return;
    free_node(tree->root);
    tree->root = NULL;
    tree->nodes = 0;
}

/* Height of the subtree under node: 0 for none, else 1 + the taller child's height. */
static size_t node_height(const VoteNode *node)
{
    size_t hl, hr;

    if (node == NULL)
        return 0;
    hl = node_height(node->left);
    hr = node_height(node->right);
    return 1 + (hl > hr ? hl : hr);
}

/* Height of the whole tree. See vote.h. */
size_t vote_tree_height(const VoteTree *tree)
{
    return tree == NULL ? 0 : node_height(tree->root);
}

/* Stage 5 for one read: Bloom check, index lookup, one vote per reference position. See vote.h. */
int vote_read(const Seq *read, const MinimizerSet *read_mins, const Index *idx,
              const BloomFilter *bloom, size_t k, VoteTree *tree, VoteStats *stats)
{
    VoteStats st;
    size_t i;

    memset(&st, 0, sizeof st);

    if (read == NULL || read_mins == NULL || idx == NULL || bloom == NULL || tree == NULL)
        return -1;
    if (k == 0 || k > 32 || idx->k != k)
        return -1;
    if (read_mins->count > 0 && read_mins->items == NULL)
        return -1;

    for (i = 0; i < read_mins->count; i++) {
        size_t read_pos = read_mins->items[i].pos;
        uint64_t hash = read_mins->items[i].hash;
        uint64_t word_id;
        int maybe;
        const IndexEntry *e;
        const PosNode *p;

        if (k > read->len || read_pos > read->len - k)   /* word must lie inside the read */
            return -1;

        st.words++;
        word_id = word_id_at(read, read_pos, k);

        maybe = bloom_maybe_contains(bloom, hash, word_id);
        if (maybe < 0)
            return -1;
        if (maybe == 0) {                     /* definitely not in the reference: skip the lookup */
            st.bloom_skipped++;
            continue;
        }

        st.lookups++;
        e = index_lookup(idx, hash, word_id, NULL);
        if (e == NULL)                        /* Bloom false positive: nothing to vote for */
            continue;

        for (p = e->pos_head; p != NULL; p = p->next) {
            long start = (long)p->pos - (long)read_pos;
            if (vote_tree_add(tree, start) != 0)
                return -1;
            st.hits++;
        }
    }

    if (stats != NULL)
        *stats = st;
    return 0;
}
