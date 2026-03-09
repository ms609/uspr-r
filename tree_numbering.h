#ifndef INCLUDE_TREE_NUMBERING
#define INCLUDE_TREE_NUMBERING

// Bridge between utree and TreeTools tree numbers.
// Converts utree ↔ tree_num_t for use in the A* search.

#include <TreeTools/tree_number.h>
#include "utree.h"
#include "unode.h"

using TreeTools::tree_num_t;

// ------------------------------------------------------------
// utree → tree_num_t
// ------------------------------------------------------------

// Postorder DFS helper: emits parent/child edge pairs in the format
// expected by edges_to_tree_number.
// node_map: maps utree node labels → 1-indexed R-convention labels.
// internal_counter: next R-style internal node label to assign (n+1..2n-1).
// parent/child/idx: output arrays and current write position.
inline void postorder_dfs(
    unode* node, unode* from,
    int n_tip,
    int* node_map,
    int& internal_counter,
    intx* parent, intx* child, int& idx) {

  // Leaf: nothing to recurse into
  if (node->get_label() >= 0) return;

  // Assign R-style label to this internal node.
  // We assign in postorder: children are numbered before parent.
  // The actual assignment happens after visiting children.
  int children_visited = 0;
  int child_r_labels[2]; // R-convention labels of the two children

  for (unode* nbr : node->const_neighbors()) {
    if (nbr == from) continue;
    // Recurse into child subtrees first (postorder)
    postorder_dfs(nbr, node, n_tip, node_map, internal_counter,
                  parent, child, idx);
    // Record this child's R-convention label
    child_r_labels[children_visited] = node_map[nbr->get_label() + n_tip];
    ++children_visited;
  }

  // Assign this internal node its R label (postorder: after children)
  int r_label = internal_counter++;
  // Store in node_map: utree internal labels are negative (-2, -3, ...),
  // mapped to node_map index = label + n_tip.
  node_map[node->get_label() + n_tip] = r_label;

  // Emit the two edges (parent, child1) and (parent, child2)
  parent[idx]     = r_label;
  child[idx]      = child_r_labels[0];
  parent[idx + 1] = r_label;
  child[idx + 1]  = child_r_labels[1];
  idx += 2;
}

// Convert a utree to a tree_num_t.
// The utree must be an unrooted bifurcating tree with leaves labeled 0..n-1.
// Returns tree_num_t() (all zeros) for trees with < 4 leaves.
inline tree_num_t utree_to_tree_number(utree& T) {
  int n_tip = 0;
  for (unode* u : T.get_leaves()) {
    if (u != nullptr) ++n_tip;
  }
  if (n_tip < 4) return tree_num_t();

  // node_map: maps utree node labels to 1-indexed R-convention labels.
  // Leaves 0..n-1 → tips 1..n.
  // Internal nodes -2,-3,... → mapped via offset: index = label + n_tip.
  // Allocated to cover labels from -(n-1) to (n-1), i.e. indices 1 to 2n-1.
  int node_map[TreeTools::TREE_NUM_MAX_TIP * 2 + 1];

  // Map leaves: utree label k → R tip (k+1)
  for (int k = 0; k < n_tip; ++k) {
    node_map[k + n_tip] = k + 1;
  }

  // Root at leaf 0 (smallest leaf) for consistent orientation.
  // Leaf 0 → R tip 1 (the root tip in the R convention).
  unode* root_leaf = T.get_leaf(0);
  if (root_leaf == nullptr) {
    // Find the actual smallest leaf
    for (int k = 0; k < n_tip; ++k) {
      if (T.get_leaf(k) != nullptr) {
        root_leaf = T.get_leaf(k);
        break;
      }
    }
  }
  if (root_leaf == nullptr) return tree_num_t();

  unode* root_internal = root_leaf->const_neighbors().front();

  // Build postorder edge arrays (n-2 internal nodes × 2 edges each)
  intx par[TreeTools::TREE_NUM_MAX_TIP * 2];
  intx chi[TreeTools::TREE_NUM_MAX_TIP * 2];
  int idx = 0;
  int internal_counter = n_tip + 1; // first R internal node label

  postorder_dfs(root_internal, root_leaf, n_tip, node_map,
                internal_counter, par, chi, idx);

  // The root internal node should have gotten label 2*n_tip - 1
  // (last internal node in postorder = root).
  // Its edges to root_leaf are NOT emitted by postorder_dfs because
  // root_leaf was the "from" parameter. We don't need them — the
  // edges_to_tree_number algorithm processes n_edge - 2 edges
  // (all except the root's pair).

  return TreeTools::edges_to_tree_number(par, chi, static_cast<intx>(n_tip));
}

// ------------------------------------------------------------
// tree_num_t → utree
// ------------------------------------------------------------

// Build a utree directly from a tree number, bypassing Newick serialization.
// Parent vector (1-indexed R convention):
//   par[i] = parent of node (i+1), for i = 0..2n-3
//   Nodes 1..n_tip are tips; n_tip+1..2n_tip-2 are internal; 2n_tip-1 is root.
// utree convention:
//   Leaves labeled 0..n_tip-1; internal nodes labeled -2, -3, ...
//   Edges are bidirectional (unrooted), all internal nodes degree 3.
//
// The R root node (2n-1) has degree 2 (rooted representation). To produce
// a proper unrooted tree, we skip the root node entirely and connect its
// two children directly — equivalent to the Newick parser's dummy/contract.
inline utree tree_number_to_utree(const tree_num_t& num, int n_tip) {
  intx par[TreeTools::TREE_NUM_MAX_TIP * 2];
  TreeTools::tree_number_to_parent(num, static_cast<intx>(n_tip), par);

  utree T;

  for (int k = 0; k < n_tip; ++k) {
    T.add_leaf(static_cast<unsigned int>(k));
  }

  // n_tip - 2 internal nodes for an unrooted binary tree.
  // R internal nodes n_tip+1 .. 2*n_tip-2 each get a utree node.
  // The R root (2*n_tip-1) is NOT created — it is dissolved below.
  for (int j = 1; j <= n_tip - 2; ++j) {
    T.add_internal_node(); // labels -2, -3, ..., -(n_tip-1)
  }

  // R node → utree label.  Root is excluded (caller must not pass it).
  // R tips 1..n_tip → utree 0..n_tip-1
  // R internal n_tip+j (j=1..n_tip-2) → utree -(j+1)
  auto r_to_utree = [n_tip](int r_node) -> int {
    if (r_node <= n_tip) return r_node - 1;
    return -(r_node - n_tip) - 1;
  };

  const int root_r = 2 * n_tip - 1;
  int root_children[2];
  int rc = 0;

  for (int i = 0; i < 2 * n_tip - 2; ++i) {
    int r_child  = i + 1;
    int r_parent = par[i];

    if (r_parent == root_r) {
      root_children[rc++] = r_child;
      continue;
    }

    unode* c = T.get_node(r_to_utree(r_child));
    unode* p = T.get_node(r_to_utree(r_parent));
    p->add_neighbor(c);
    c->add_neighbor(p);
  }

  // Connect the root's two children directly (unrooting)
  unode* a = T.get_node(r_to_utree(root_children[0]));
  unode* b = T.get_node(r_to_utree(root_children[1]));
  a->add_neighbor(b);
  b->add_neighbor(a);

  T.set_smallest_leaf(0);
  return T;
}

#endif
