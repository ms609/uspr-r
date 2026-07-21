/*******************************************************************************
uspr.h

Unrooted SPR distance computation and data structures

Copyright 2018 Chris Whidden
cwhidden@fredhutch.org
https://github.com/cwhidden/uspr
May 1, 2018
Version 1.0.1

This file is part of uspr.

uspr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

uspr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uspr.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef INCLUDE_USPR
#define INCLUDE_USPR


#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>
#include <list>
#include <memory>
#include <ctime>
#include <cstdlib>
#include "utree.h"
#include "unode.h"
#include "uforest.h"
#include "tbr.h"
#include "uspr_neighbors.h"
#include "tree_numbering.h"
#include "uspr_neighbors_numbered.h"
#include "../spr/utree_splits.h"

//#define DEBUG_USPR 1
#ifdef DEBUG_USPR
	#define debug_uspr(x) x
#else
	#define debug_uspr(x)
#endif

// options
bool USE_TBR_APPROX_ESTIMATE = true;
bool USE_TBR_ESTIMATE = true;
bool USE_REPLUG_ESTIMATE = true;

// classes

typedef enum {REPLUG, TBR, TBR_APPROX, BFS} estimator_t;
string estimator_t_name[] = {"REPLUG", "TBR", "TBR_APPROX", "BFS"};

// distance estimate (string-based, for fallback)
class tree_distance {
	public:
	int cost;
	int estimate;
	int distance;
	string tree;
	estimator_t estimator;

	tree_distance(int c, int d, string t, estimator_t e) {
		cost = c;
		estimate = d;
		distance = c + d;
		tree = t;
		estimator = e;
	}

	friend bool operator < (tree_distance a, tree_distance b);
	friend bool operator <= (tree_distance a, tree_distance b);

};

bool operator < (tree_distance a, tree_distance b) {
	if (a.distance == b.distance) {
		if (a.estimate == b.estimate) {
			return a.estimator < b.estimator;
		}
		else {
			return a.estimate < b.estimate;
		}
	}
	return a.distance < b.distance; }
bool operator <= (tree_distance a, tree_distance b) {
	if (a.distance == b.distance) {
		if (a.estimate == b.estimate) {
			return a.estimator <= b.estimator;
		}
		else {
			return a.estimate <= b.estimate;
		}
	}
	return a.distance <= b.distance;
}

// distance estimate (tree-number-based, for n <= 51)
class tree_distance_num {
	public:
	int cost;
	int estimate;
	int distance;
	tree_num_t tree_number;
	estimator_t estimator;

	tree_distance_num(int c, int d, tree_num_t t, estimator_t e) {
		cost = c;
		estimate = d;
		distance = c + d;
		tree_number = t;
		estimator = e;
	}

	friend bool operator < (const tree_distance_num& a, const tree_distance_num& b);
	friend bool operator > (const tree_distance_num& a, const tree_distance_num& b);
};

bool operator < (const tree_distance_num& a, const tree_distance_num& b) {
	if (a.distance != b.distance) return a.distance < b.distance;
	if (a.estimate != b.estimate) return a.estimate < b.estimate;
	return a.estimator < b.estimator;
}
bool operator > (const tree_distance_num& a, const tree_distance_num& b) {
	return b < a;
}


// function prototypes
int uspr_distance(uforest &T1, uforest &T2);
int uspr_distance_string_based(uforest &T1, uforest &T2);
int uspr_distance_numbered(uforest &T1, uforest &T2, int n_tip);

// functions

// Main entry point: dispatches to numbered or string-based A* search.
int uspr_distance(uforest &T1_original, uforest &T2_original) {

	uforest T1 = uforest(T1_original);
	uforest T2 = uforest(T2_original);

	// normalize tree order
	T1.normalize_order();
	T2.normalize_order();

	// check if trees are equal
	if (utree(T1).str() == utree(T2).str()) {
		return 0;
	}

	// leaf reduction
	map<string, int> label_map = map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();
	leaf_reduction(&T1, &T2, &label_map, &reverse_label_map);
	T1.normalize_order();
	T2.normalize_order();

	debug_uspr(
		Rcout << "T1R: " << T1 << endl;
		Rcout << "T2R: " << T2 << endl;
	)

	// exact lookup for small reduced trees (4-9 leaves)
	{
		int lookup_result = spr_lookup::lookup_utrees(T1, T2);
		if (lookup_result >= 0) {
			return lookup_result;
		}
	}

	int n_tip = spr_lookup::count_leaves(T1);

	// Use tree-number A* for n <= 51, fall back to strings for larger trees
	if (n_tip <= TreeTools::TREE_NUM_MAX_TIP) {
		return uspr_distance_numbered(T1, T2, n_tip);
	} else {
		return uspr_distance_string_based(T1, T2);
	}
}

// Tree-number-based A* search (n <= 51 leaves after reduction).
// Uses uint256 tree numbers instead of Newick strings for:
// - O(1) equality comparison and hashing
// - zero heap allocation for tree identity
// - no redundant serialization in neighbor generation
int uspr_distance_numbered(uforest &T1, uforest &T2, int n_tip) {

	// visited set: O(1) amortized lookup via hash
	unordered_set<tree_num_t, tree_num_hash> visited_trees;

	// target tree number
	tree_num_t target_number = utree_to_tree_number(T2);

	// priority queue: binary min-heap (better cache behaviour than multiset)
	priority_queue<tree_distance_num,
	               vector<tree_distance_num>,
	               greater<tree_distance_num>> distance_priority_queue;

	// start with T1
	tree_num_t start_number = utree_to_tree_number(T1);
	visited_trees.insert(start_number);
	distance_priority_queue.push(tree_distance_num(0, 1, start_number, BFS));

	// final estimator
	estimator_t final_estimator = BFS;
	if (USE_TBR_APPROX_ESTIMATE) {
		final_estimator = TBR_APPROX;
	}
	if (USE_TBR_ESTIMATE) {
		final_estimator = TBR;
	}
	if (USE_REPLUG_ESTIMATE) {
		final_estimator = REPLUG;
	}

	unsigned long interrupt_counter = 0;
	while (!distance_priority_queue.empty()) {
		if ((++interrupt_counter & 1023UL) == 0UL) {
			Rcpp::checkUserInterrupt();
		}

		const tree_distance_num top = distance_priority_queue.top();
		distance_priority_queue.pop();

		int cost = top.cost;
		tree_num_t tn = top.tree_number;
		estimator_t prev_estimator = top.estimator;

		// decode tree number → uforest directly (no Newick roundtrip)
		// distances_from_leaf_decorator is NOT called here: normalize_order()
		// and utree_to_tree_number() don't use get_distance(), and every
		// estimator function (tbr_high_lower_bound, tbr_distance,
		// replug_distance) makes its own copy and calls root() +
		// distances_from_leaf_decorator internally.
		uforest T(tree_number_to_utree(tn, n_tip));
		T.normalize_order();

		if (prev_estimator != final_estimator) {
			// Phase 3: try exact lookup on reduced pair at first pop (BFS)
			if (prev_estimator == BFS) {
				uforest T_copy(T);
				uforest T2_copy(T2);

				vector<int> red_leaves = T_copy.find_leaves();
				nodemapping twins(red_leaves);
				map<int, int> sibling_pairs = T_copy.find_sibling_pairs();
				T_copy.root(T_copy.get_smallest_leaf());
				T2_copy.root(T2_copy.get_smallest_leaf());
				distances_from_leaf_decorator(T_copy, T_copy.get_smallest_leaf());
				distances_from_leaf_decorator(T2_copy, T2_copy.get_smallest_leaf());
				for (unode* u : T_copy.get_leaves())
					if (u != nullptr) u->set_terminal(true);
				for (unode* u : T2_copy.get_leaves())
					if (u != nullptr) u->set_terminal(true);

				leaf_reduction_hlpr(T_copy, T2_copy, twins, sibling_pairs);

				unode* anchor1 = T_copy.get_node(
					T_copy.get_smallest_leaf())->find_uncontracted_node();
				unode* root1 = nullptr;
				if (!anchor1->get_terminal()) {
					for (unode* nbr : anchor1->get_neighbors()) {
						if (nbr->get_terminal()) { root1 = nbr; break; }
					}
				}
				if (root1 != nullptr) {
					vector<int> t1_labels;
					t1_labels.push_back(root1->get_label());
					spr_lookup::collect_terminal_labels(anchor1, root1, t1_labels);

					int n_red = static_cast<int>(t1_labels.size());
					if (n_red >= 4 && n_red <= 9) {
						map<int, int> remap1, remap2;
						bool remap_ok = true;
						for (int ri = 0; ri < n_red; ri++) {
							remap1[t1_labels[ri]] = ri;
							int t2_lab = twins.get_forward(t1_labels[ri]);
							if (t2_lab == -1) { remap_ok = false; break; }
							remap2[t2_lab] = ri;
						}
						if (remap_ok) {
							int root2_lab = twins.get_forward(root1->get_label());
							unode* root2 = T2_copy.get_node(root2_lab);
							int exact = spr_lookup::lookup_reduced_utrees(
								T_copy, T2_copy, root1, root2,
								remap1, remap2, n_red);
							if (exact >= 0) {
								distance_priority_queue.push(
									tree_distance_num(cost, exact, tn, final_estimator));
								continue;
							}
						}
					}
				}
			}
			// Compute the next estimate, re-queue with SAME tree number
			int distance = 1;
			if (prev_estimator > TBR_APPROX &&
					USE_TBR_APPROX_ESTIMATE) {
				distance = tbr_high_lower_bound(T, T2);
				distance_priority_queue.push(tree_distance_num(cost, distance, tn, TBR_APPROX));
			}
			else if (prev_estimator > TBR &&
					USE_TBR_ESTIMATE) {
				distance = tbr_distance(T, T2);
				distance_priority_queue.push(tree_distance_num(cost, distance, tn, TBR));
			}
			else if (prev_estimator > REPLUG &&
					USE_REPLUG_ESTIMATE) {
				distance = replug_distance(T, T2);
				distance_priority_queue.push(tree_distance_num(cost, distance, tn, REPLUG));
			}
			continue;
		}

		// Final estimator: expand neighbors
		list<tree_num_t> neighbors = get_neighbors_numbered(&T, &visited_trees);
		debug_uspr(
			Rcout << "examining " << neighbors.size() << " neighbors" << endl;
		)
		for (const tree_num_t& nbr_num : neighbors) {
			if (nbr_num == target_number) {
				debug_uspr(
					Rcout << "examined " << visited_trees.size() << " trees" << endl;
				)
				return cost + 1;
			}
			distance_priority_queue.push(tree_distance_num(cost + 1, 1, nbr_num, BFS));
		}
	}

	return -1;
}

// String-based A* search (fallback for n > 51).
int uspr_distance_string_based(uforest &T1, uforest &T2) {

	set<string> visited_trees = set<string>();
	string target = utree(T2).str();
	multiset<tree_distance> distance_priority_queue = multiset<tree_distance>();

	visited_trees.insert(T1.str());
	distance_priority_queue.insert(tree_distance(0, 1, utree(T1).str(), BFS));

	estimator_t final_estimator = BFS;
	if (USE_TBR_APPROX_ESTIMATE) {
		final_estimator = TBR_APPROX;
	}
	if (USE_TBR_ESTIMATE) {
		final_estimator = TBR;
	}
	if (USE_REPLUG_ESTIMATE) {
		final_estimator = REPLUG;
	}

	unsigned long interrupt_counter = 0;
	while (!distance_priority_queue.empty()) {
		if ((++interrupt_counter & 1023UL) == 0UL) {
			Rcpp::checkUserInterrupt();
		}

		multiset<tree_distance>::iterator it = distance_priority_queue.begin();

		int cost = it->cost;
		string tree = it->tree;
		estimator_t prev_estimator = it->estimator;
		distance_priority_queue.erase(it);

		uforest T = uforest(tree);
		distances_from_leaf_decorator(T, T.get_smallest_leaf());
		T.normalize_order();

		if (prev_estimator != final_estimator) {
			int distance = 1;
			if (prev_estimator > TBR_APPROX &&
					USE_TBR_APPROX_ESTIMATE) {
				distance = tbr_high_lower_bound(T, T2);
				distance_priority_queue.insert(tree_distance(cost, distance, tree, TBR_APPROX));
			}
			else if (prev_estimator > TBR &&
					USE_TBR_ESTIMATE) {
				distance = tbr_distance(T, T2);
				distance_priority_queue.insert(tree_distance(cost, distance, tree, TBR));
			}
			else if (prev_estimator > REPLUG &&
					USE_REPLUG_ESTIMATE) {
				distance = replug_distance(T, T2);
				distance_priority_queue.insert(tree_distance(cost, distance, tree, REPLUG));
			}
			continue;
		}

		list<utree> neighbors = get_neighbors(&T, &visited_trees);
		for (utree tree : neighbors) {
			string tree_string = tree.str();
			if (tree_string == target) {
				return cost + 1;
			}
			distance_priority_queue.insert(tree_distance(cost + 1, 1, tree_string, BFS));
		}
	}

	return -1;
}

#endif
