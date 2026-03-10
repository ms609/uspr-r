/*******************************************************************************
unode.h

Unrooted tree node data structure and functions

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

#ifndef INCLUDE_UNODE
#define INCLUDE_UNODE

//#define DEBUG 1
#ifdef DEBUG
	#define debug(x) x
#else
	#define debug(x)
#endif

#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <cstdio>
#include <climits>

using namespace std;

class unode;

// Fixed-capacity inline list for neighbor pointers.
// Unrooted binary tree nodes have degree <= 3, so a 3-element array
// replaces std::list<unode*> with zero heap allocation.
class neighbor_list {
	// Unrooted binary tree nodes have degree <= 3. Capacity 4 accommodates
	// the transient dummy parent added during Newick parsing of the
	// trifurcating root (3 children + 1 dummy, removed immediately).
	unode* data_[4];
	int size_ = 0;
public:
	using iterator = unode**;
	using const_iterator = unode* const*;
	using reverse_iterator = std::reverse_iterator<iterator>;
	using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	neighbor_list() = default;
	neighbor_list(const neighbor_list&) = default;
	neighbor_list& operator=(const neighbor_list&) = default;

	iterator begin() { return data_; }
	iterator end()   { return data_ + size_; }
	const_iterator begin() const { return data_; }
	const_iterator end()   const { return data_ + size_; }

	reverse_iterator rbegin() { return reverse_iterator(end()); }
	reverse_iterator rend()   { return reverse_iterator(begin()); }
	const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
	const_reverse_iterator rend()   const { return const_reverse_iterator(begin()); }

	unode* front() const { return data_[0]; }
	unode* back()  const { return data_[size_ - 1]; }

	bool empty() const { return size_ == 0; }
	int size()   const { return size_; }

	void push_front(unode* n) {
		for (int i = size_; i > 0; --i) data_[i] = data_[i - 1];
		data_[0] = n;
		++size_;
	}

	void push_back(unode* n) {
		data_[size_++] = n;
	}

	// Remove first occurrence by value, preserving order.
	void remove(unode* n) {
		for (int i = 0; i < size_; ++i) {
			if (data_[i] == n) {
				for (int j = i; j < size_ - 1; ++j)
					data_[j] = data_[j + 1];
				--size_;
				return;
			}
		}
	}

	void clear() { size_ = 0; }
};


class unode {
	private:
	int label;
	neighbor_list neighbors;
	vector<unode *> contracted_neighbors;
	int component;
	bool terminal;
	int distance;
	bool b_protected;
	bool phi;

	public:
	unode() {
		label = -1;
		component = -1;
		terminal = false;
		distance = -1;
		b_protected = false;
		phi = false;
	}
	unode(int l) {
		label = l;
		component = -1;
		terminal = false;
		distance = -1;
		b_protected = false;
		phi = false;
	}
	unode(const unode &n, bool include_neighbors = true) {
		label = n.label;
		// don't include neighbors when copying as they will be updated later
		if (include_neighbors) {
			neighbors = n.neighbors;
			contracted_neighbors = n.contracted_neighbors;
		}
		component = n.component;
		terminal = n.terminal;
		distance = n.distance;
		b_protected = n.b_protected;
		phi = n.phi;
	}
	~unode() {
	}

	void add_neighbor(unode *n) {
		if (neighbors.size() > 0 && neighbors.front()->get_distance() > n->get_distance()) {
			neighbors.push_front(n);
		}
		else {
			neighbors.push_back(n);
		}
	}

	void add_contracted_neighbor(unode *n) {
		contracted_neighbors.push_back(n);
	}

	void add_parent(unode *n) {
		neighbors.push_front(n);
	}

	bool remove_neighbor(unode *n) {
		for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
			if (*it == n) {
				neighbors.remove(n);
				return true;
			}
		}
		return false;
	}

	bool remove_contracted_neighbor(unode *n) {
		auto it = std::find(contracted_neighbors.begin(),
		                    contracted_neighbors.end(), n);
		if (it != contracted_neighbors.end()) {
			contracted_neighbors.erase(it);
			return true;
		}
		return false;
	}

	bool contract_neighbor(unode *n) {
		bool ret = remove_neighbor(n);
		if (ret) {
			contracted_neighbors.push_back(n);
		}
		return ret;
	}

	string str(map<int, string> *reverse_label_map = NULL) const {
		stringstream ss;
		if (phi) {
			ss << "*";
		}
		else {
			if (reverse_label_map != NULL &&
					reverse_label_map->find(label) != reverse_label_map->end()) {
				ss << (*reverse_label_map)[label];
			}
			else {
				ss << label;
			}
		}
		return ss.str();
	}

	bool operator ==(const unode n) const {
		return this->label == n.label;
	}
	bool operator !=(const unode n) const {
		return this->label != n.label;
	}


	int get_label() const {
		return label;
	}

	const neighbor_list &const_neighbors() const {
		return neighbors;
	}

	const vector<unode *> &const_contracted_neighbors() const {
		return contracted_neighbors;
	}

	neighbor_list &get_neighbors() {
		return neighbors;
	}

	vector<unode *> &get_contracted_neighbors() {
		return contracted_neighbors;
	}

	int get_num_neighbors() {
		return neighbors.size();
	}

	int get_num_all_neighbors() {
		return neighbors.size() + static_cast<int>(contracted_neighbors.size());
	}

	bool is_leaf() {
		return (neighbors.size() == 1);
	}

	void set_component(int c) {
		component = c;
	}

	int get_component() {
		return component;
	}

	void set_phi(bool b) {
		phi = b;
	}

	bool is_phi() {
		return phi;
	}

	void root(int l) {
		unode *p = NULL;
		for (unode *n : neighbors) {
			if (n->get_label() == l) {
				p = n;
			}
			else {
				n->root(get_label());
			}
		}
		if (p != NULL) {
			neighbors.remove(p);
			neighbors.push_front(p);
		}
	}

	void rotate(int l) {
		unode *p = NULL;
		for (unode *n : neighbors) {
			if (n->get_label() == l) {
				p = n;
			}
		}
		if (p != NULL) {
			neighbors.remove(p);
			neighbors.push_front(p);
		}
	}

	unode *get_parent() {
		if (neighbors.empty()) {
			return NULL;
		}
		return neighbors.front();
	}

	bool is_adjacent(unode *a) {
//		Rcout << label << "->is_adjacent(" << a->get_label() << ")" << endl;
		for (unode *n : neighbors) {
			if (n == a) {
//				Rcout << "true" << endl;
				return true;
			}
		}
		for (unode *n : contracted_neighbors) {
			if (n == a) {
//				Rcout << "true" << endl;
				return true;
			}
		}
//		Rcout << "false" << endl;
		return false;
	}

	unode *get_neighbor_not(unode *a) {
		return get_neighbor_not(a, a);
	}

	unode *get_neighbor_not(unode *a, unode *b) {
		for (auto x = neighbors.rbegin(); x != neighbors.rend(); ++x) {
			if (*x != a && *x != b) {
				return *x;
			}
		}
		return NULL;
	}

	void set_terminal(bool t) {
		terminal = t;
	}

	bool get_terminal() {
		return terminal;
	}

	unode *get_sibling() {
		unode *parent = get_parent();
		int count = 0;
		for (unode *x : parent->get_neighbors()) {
			if (count > 0 && x != this) {
				return x;
			}
			count++;
		}
		return parent; // Never reached, but appeases compiler
	}

	void clear_neighbors() {
		neighbors.clear();
	}


	void clear_contracted_neighbors() {
		contracted_neighbors.clear();
	}

	void uncontract_neighbors() {
		for (unode *x: contracted_neighbors) {
			add_neighbor(x);
		}
		clear_contracted_neighbors();
	}

	void uncontract_subtree(unode *last = NULL) {
		for (unode *n : neighbors) {
			if (last == NULL || n != last) {
				n->uncontract_subtree(this);
			}
		}
		for (unode *n : contracted_neighbors) {
			if (last == NULL || n != last) {
				n->uncontract_subtree(this);
			}
		}
		uncontract_neighbors();
	}

	unode *contract_degree_two_subtree(unode *last = NULL) {
		debug(
			Rcout << label << ".contract_degree_two_subtree()" << endl;
		)
		neighbor_list neighbor_copy(neighbors);
		for (unode *n : neighbor_copy) {
			if (last == NULL || n != last) {
				n->contract_degree_two_subtree(this);
			}
		}
		return contract();
	}

	unode *contract() {
		debug(
			Rcout << "n: " << neighbors.size() << endl;
			Rcout << "c_n: " << contracted_neighbors.size() << endl;
		)
		if (neighbors.size() == 1 && contracted_neighbors.empty()) {
			unode *p = neighbors.front();
			if (p->is_leaf() && this->get_label() < -1) {
				p->remove_neighbor(this);
				this->remove_neighbor(p);
				if (component > -1) {
					p->set_component(component);
				}
				if (is_protected()) {
					p->set_protected(true);
				}
				return p;
			}
		}
//		/*
		else if (neighbors.size() == 0 && contracted_neighbors.size() == 2) {
//			uncontract_neighbors();
			unode *p = contracted_neighbors.front();
			unode *c = contracted_neighbors[1];
			debug(
				Rcout << "contracting:" << endl;
				Rcout << p << "\t" << p->get_num_all_neighbors() << endl;
				Rcout << c << "\t" << c->get_num_all_neighbors() << endl;
			)
			if (p->get_num_all_neighbors() < c->get_num_all_neighbors()) {
				unode *temp = p;
				p = c;
				c = temp;
			}
			if (p->get_num_all_neighbors() > 1) {
				clear_contracted_neighbors();
				p->remove_neighbor(this);
				p->remove_contracted_neighbor(this);
				c->remove_neighbor(this);
				c->remove_contracted_neighbor(this);
				c->add_parent(p);
				p->add_contracted_neighbor(c);
				if (p->get_distance() > distance &&
						c->get_distance() > distance) {
					p->set_distance(distance-1);
					c->set_distance(distance);
				}
				else {
					c->set_distance(p->get_distance()+1);
				}
				if (!get_terminal()) {
					p->set_terminal(false);
				}
				else {
					p->set_terminal(true);
				}
				if (component > -1) {
					p->set_component(component);
				}
				if (is_protected()) {
					c->set_protected(true);
				}
				return p;
			}
		}
//		*/
		else if (neighbors.size() == 2 && contracted_neighbors.empty()) {
			unode *p = neighbors.front();
			unode *c = *(neighbors.begin() + 1);
			debug(
				Rcout << "contracting:" << endl;
				Rcout << p << endl;
				Rcout << c << endl;
			)
			if (!p->is_leaf() ||
						!(p->get_contracted_neighbors().empty()) ||
						!c->is_leaf()) {
				clear_neighbors();
				p->remove_neighbor(this);
				c->remove_neighbor(this);
				c->add_parent(p);
				p->add_neighbor(c);
				if (p->get_distance() > distance &&
						c->get_distance() > distance) {
					p->set_distance(distance-1);
					c->set_distance(distance);
				}
				else {
					c->set_distance(p->get_distance()+1);
				}
				if (!get_terminal()) {
					p->set_terminal(false);
				}
				if (component > -1) {
					p->set_component(component);
				}
				if (is_protected()) {
					c->set_protected(true);
				}
				return p;
			}
		}
		return this;
	}

	void set_distance(int d) {
		distance = d;
	}

	int get_distance() {
		return distance;
	}
	bool is_singleton() {
		if (neighbors.size() == 0) {
			return true;
		}
		return false;
	}

	bool is_protected() {
		return b_protected;
	}

	void set_protected(bool b) {
		b_protected = b;
	}

	void get_connected_nodes(list<unode *> &connected_nodes, unode *last = NULL) {
		for (unode *n : neighbors) {
			if (last == NULL || n != last) {
				n->get_connected_nodes(connected_nodes,this);
			}
		}
		for (unode *n : contracted_neighbors) {
			if (last == NULL || n != last) {
				n->get_connected_nodes(connected_nodes,this);
			}
		}
		connected_nodes.push_back(this);
	}

	int normalize_order_hlpr(unode *prev = NULL) {
		// return leaf label
		if (label >= 0 && prev != NULL) {
			return label;
		}
		unode *parent = NULL;
		int min_descendant = INT_MAX;

		// Active neighbors: at most 3 children.
		// Inline array + insertion sort (already inline — no extra allocation).
		unode* ch[3];
		int ch_min[3];
		int nch = 0;

		for (unode *n : neighbors) {
			if (n != prev) {
				int m = n->normalize_order_hlpr(this);
				int i = nch++;
				while (i > 0 && ch_min[i - 1] > m) {
					ch[i]     = ch[i - 1];
					ch_min[i] = ch_min[i - 1];
					--i;
				}
				ch[i]     = n;
				ch_min[i] = m;
				if (m < min_descendant) min_descendant = m;
			}
			else {
				parent = n;
			}
		}

		// re-add in sorted order
		clear_neighbors();
		if (parent != NULL) add_neighbor(parent);
		for (int i = 0; i < nch; ++i) add_neighbor(ch[i]);

		// Contracted neighbors
		if (!contracted_neighbors.empty()) {
			vector<pair<int, unode*>> ordered;
			ordered.reserve(contracted_neighbors.size());
			for (unode *n : contracted_neighbors) {
				int m = n->normalize_order_hlpr(this);
				ordered.push_back({m, n});
				if (m < min_descendant) min_descendant = m;
			}
			sort(ordered.begin(), ordered.end());
			clear_contracted_neighbors();
			for (auto& p : ordered) add_contracted_neighbor(p.second);
		}

		return min_descendant;
	}



	// normalize branching order by smallest subtree leaf
	// guaranteed unique if started at the smallest leaf
	void normalize_order() {
		// normalize order
		normalize_order_hlpr();
	}

	unode *find_uncontracted_node() {
		unode *ret = this;
		unode *prev = this;
		while (ret->is_leaf()) {
			unode *next = ret->get_parent();
			if (prev == next) {
				return ret;
			}
			prev = ret;
			ret = next;
		}
		return ret;
	}

};

#endif
