/*******************************************************************************
uspr_neighbors_numbered.h

SPR neighbor generation returning tree numbers instead of strings/utrees.
Eliminates string serialization and deep tree copies from neighbor generation.

Based on uspr_neighbors.h by Chris Whidden.
*******************************************************************************/

#ifndef INCLUDE_USPR_NEIGHBORS_NUMBERED
#define INCLUDE_USPR_NEIGHBORS_NUMBERED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>

#include "utree.h"
#include "tree_numbering.h"

using namespace std;
using TreeTools::tree_num_t;
using TreeTools::tree_num_hash;

// Function prototypes
list<tree_num_t> get_neighbors_numbered(utree *T,
    unordered_set<tree_num_t, tree_num_hash> *known_trees);
void get_neighbors_numbered(utree *T, unode *prev, unode *current,
    list<tree_num_t> &neighbors,
    unordered_set<tree_num_t, tree_num_hash> *known_trees);
void get_neighbors_numbered(utree *T, unode *x, unode *y,
    unode *prev, unode *current,
    list<tree_num_t> &neighbors,
    unordered_set<tree_num_t, tree_num_hash> *known_trees);
void add_neighbor_numbered(utree *T, unode *x, unode *y,
    unode *w, unode *z,
    list<tree_num_t> &neighbors,
    unordered_set<tree_num_t, tree_num_hash> *known_trees);


list<tree_num_t> get_neighbors_numbered(utree *T,
    unordered_set<tree_num_t, tree_num_hash> *known_trees) {
	list<tree_num_t> neighbors;
	unode *root = T->get_node(T->get_smallest_leaf());
	get_neighbors_numbered(T, NULL, root, neighbors, known_trees);
	return neighbors;
}

// enumerate the source edges
void get_neighbors_numbered(utree *T, unode *prev, unode *current,
    list<tree_num_t> &neighbors,
    unordered_set<tree_num_t, tree_num_hash> *known_trees) {
	list<unode *> c_neighbors = current->get_neighbors();
	for (unode *next : c_neighbors) {
		if (next != prev) {
			get_neighbors_numbered(T, current, next, neighbors, known_trees);
		}
	}
	if (prev != NULL) {
		get_neighbors_numbered(T, prev, current, prev, current, neighbors, known_trees);
		get_neighbors_numbered(T, current, prev, current, prev, neighbors, known_trees);
	}
}

// enumerate the target edges
void get_neighbors_numbered(utree *T, unode *x, unode *y,
    unode *prev, unode *current,
    list<tree_num_t> &neighbors,
    unordered_set<tree_num_t, tree_num_hash> *known_trees) {
	list<unode *> c_neighbors = current->get_neighbors();
	for (unode *next : c_neighbors) {
		if (next != prev) {
			get_neighbors_numbered(T, x, y, current, next, neighbors, known_trees);
		}
	}
	if (prev != NULL) {
		add_neighbor_numbered(T, x, y, prev, current, neighbors, known_trees);
	}
}

void add_neighbor_numbered(utree *T, unode *x, unode *y,
    unode *w, unode *z,
    list<tree_num_t> &neighbors,
    unordered_set<tree_num_t, tree_num_hash> *known_trees) {
	// check for duplicate SPR moves
	if (x == y || y == w || y == z) {
		return;
	}
	if (w == y->get_parent()->get_parent() &&
			z == y->get_parent()) {
		return;
	}
	if (z == y->get_parent()->get_parent() &&
			w == y->get_parent()) {
		return;
	}
	if (z->get_parent() == y &&
			w->get_parent() == z) {
		return;
	}
	if (w->get_parent() == y &&
			z->get_parent() == w) {
		return;
	}

	// node info so the uspr can be reversed
	unode *yprime = NULL;
	unode *y1 = NULL;
	unode *y2 = NULL;

	// apply the spr
	T->uspr(x, y, w, z, &yprime, &y1, &y2);
	// normalize the tree
	distances_from_leaf_decorator(*T, T->get_smallest_leaf());
	T->normalize_order();

	// compute tree number instead of string
	tree_num_t num = utree_to_tree_number(*T);
	if (known_trees->find(num) == known_trees->end()) {
		known_trees->insert(num);
		neighbors.push_back(num);
	}

	// revert the SPR
	T->uspr(x, yprime, y1, y2);
	distances_from_leaf_decorator(*T, T->get_smallest_leaf());
	T->normalize_order();
}

#endif
