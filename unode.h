#ifndef INCLUDE_UNODE
#define INCLUDE_UNODE

#include <list>
#include <sstream>

using namespace std;

class unode;

class unode {
	private:
	int label;
	list<unode *> neighbors;
	int num_neighbors;

	public:
	unode() {
		label = -1;
		neighbors = list<unode *>();
		num_neighbors = 0;
	}
	unode(int l) {
		label = l;
		neighbors = list<unode *>();
		num_neighbors = 0;
	}

	void add_neighbor(unode *n) {
		neighbors.push_back(n);
		num_neighbors++;
	}

	bool remove_neighbor(unode *n) {
		list<unode *>::iterator i;
		bool result = false;
		for(i = neighbors.begin(); i != neighbors.end(); i++) {
			if ((*i) == n) {
				neighbors.remove(*i);
				result = true;
				num_neighbors--;
				break;
			}
		}
		return result;
	}

	string str() const {
		stringstream ss;
		ss << label;
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

	const list<unode *> &const_neighbors() const {
		return neighbors;
	}

};

#endif
