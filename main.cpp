// Jarno Leppänen 13.11.2013

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <tuple>
#include <array>
#include <stack>
#include <cstring>
#include <set>
#include <random>

using namespace std;
static const int inf = numeric_limits<int>::max();

struct progress {
	int n;
	int p;
	int i;
	progress(int n) {
		reset(n);
	}
	~progress() {
		stop();
	}
	void stop() {
		cerr << "\r         \r";
	}
	void reset(int n) {
		this->n = n;
		this->p = -1;
		this->i = -1;
	}
	void operator()() {
		int q = 100 * ++i / n;
		if (q == p)
			return;
		p = q;
		cerr << '\r' << p << "% done";
	}
};

struct shortest_paths {
	int n;
	vector<int> m;

	bool load_file(char const *filename) {
		// read adjacency matrix
		ifstream f(filename);
		if (!f.is_open()) {
			cerr << "can't open file " << filename << endl;
			return false;
		}
		cerr << "reading graph data from file " << filename << "..." << endl;
		int en;
		f >> n;
		f >> en;
		if (f.fail()) {
			cerr << "can't read number of vertices/edges" << endl;
			return false;
		}

		// read to distance matrix
		m = vector<int>(n * n, inf);
		// set diagonal to zero
		for (int i = 0, j = 0; i < n; ++i, j += n + 1) {
			m[j] = 0;
		}
		// set upper triangular to -1
		int k = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = i + 1; j < n; ++j) {
				m[i*n + j] = -1;
			}
		}
		progress p(en);
		for (int i = 0; i < en; ++i) {
			p();
			int u, v, w;
			f >> u >> v >> w;
			if (u > v) {
				m[u*n + v] = w;
			} else {
				m[v*n + u] = w;
			}
		}
		if (f.fail()) {
			cerr << "can't read edge data" << endl;
		}
		return true;
	}

	inline void update(int i, int j, int ir, int jr, int k, int w) {
		// check for overflow (max + max == -2)
		if (w >= 0 && m[ir + j] > w) {
			m[ir + j] = w;
			m[jr + i] = k;
		}
	}

	bool calculate_shortest_paths() {
		// symmetric floyd warshall with path reconstruction
		// shortest paths graph satisfies triangle equality!
		//
		// lower triangular matrix will consist of shortest path 
		// lenghts between vertices while upper triangular matrix
		// will be the path reconstruction matrix
		//
		// premature optimization is the root of all evil...
		cerr << "calculating shortest paths..." << endl;
		progress p(n);
		for (int k = 0, kr = 0; k < n; ++k, kr += n) {
			p();
			// i < k
			int i = 0, ir = 0;
			for (; i < k; ++i, ir += n) {
				// j < i < k
				for (int j = 0, jr = 0; j <= i; ++j, jr += n) {
					int w = m[kr + i] + m[kr + j];
					update(i, j, ir, jr, k, w);
				}
			}
			// i >= k
			for (; i < n; ++i, ir += n) {
				int j = 0, jr = 0;
				// j < k
				for (; j < k; ++j, jr += n) {
					int w = m[ir + k] + m[kr + j];
					update(i, j, ir, jr, k, w);
				}
				// j >= k
				for (; j < i; ++j, jr += n) {
					int w = m[ir + k] + m[jr + k];
					update(i, j, ir, jr, k, w);
				}
			}
		}
		// check negative cycles
		for (int i = 0, j = 0; i < n; ++i, j += n + 1) {
			if (m[j] != 0) {
				cerr << "graph contains negative cycles" << endl;
				return false;
			}
		}
		return true;
	}

	void print_matrix() {
		int k = 0;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (m[k] == inf) {
					cerr << '~';
				}
				else {
					cerr << m[k];
				}
				++k;
				cerr << ' ';
			}
			cerr << endl;
		}
	}

	int dist(int i, int j) const {
		if (i < j) {
			return m[j*n + i];
		}
		return m[i*n + j];
	}

	int next(int i, int j) const {
		if (i < j) {
			return m[i*n + j];
		}
		return m[j*n + i];
	}

	// floyd warshall path reconstruction
	void print_path_(int i, int j) const {
		if (dist(i, j) == inf) {
			cerr << "error finding path";
			return;
		}
		int k = next(i, j);
		if (k != -1) {
			print_path_(i, k);
			cout << k << endl;
			print_path_(k, j);
		}
	}

	void print_path(int i, int j) const {
		if (i == j)
			return;
		cout << i << endl;
		print_path_(i, j);
	}
};

// helper structs
struct edge {
	int w;
	int i;
	int j;
	edge(int w, int i, int j) : w(w), i(i), j(j) {}
};
struct edge_comp {
	bool operator()(edge const& a, edge const& b) {
		return a.w > b.w;
	}
};
struct neighbors {
	int d;
	array<int, 2> n;
	neighbors() : d(0), n{{-1, -1}} {} // zomg c++11
	inline void add(int k) {
		n[d++] = k;
	}
	inline int& operator[](int i) {
		return n[i];
	}
	inline int const& operator[](int i) const {
		return n[i];
	}
	inline void flip() {
		swap(n[0], n[1]);
	}
};

struct greedy_tsp {

	// local variables
	shortest_paths const& sp;
	priority_queue<edge, vector<edge>, edge_comp> q;
	vector<neighbors> ns;

	greedy_tsp(shortest_paths const& sp) : sp(sp), ns(sp.n) {}

	void sort_edges() {
		cerr << "sorting edges..." << endl;
		progress p(sp.n);
		for (int i = 0; i < sp.n; ++i) {
			p();
			for (int j = 0; j < i; ++j) {
				q.push(edge(sp.dist(i, j), i, j));
			}
		}
	}

	void find_greedy() {
		cerr << "constructing greedy tsp..." << endl;
		progress p(q.size());
		while (!q.empty()) {
			p();
			edge e(q.top());
			q.pop();
			// check degree
			if (ns[e.i].d == 2 || ns[e.j].d == 2)
				continue;
			// check cycles smaller than n
			int p = e.j;
			int i = e.i;
			int k = 0;
			for (;;) {
				if (++k == sp.n || i == -1) {
					ns[e.i].add(e.j);
					ns[e.j].add(e.i);
					break;
				}
				if (i == e.j) {
					break;
				}
				// next
				if (ns[i][0] == p) {
					p = i;
					i = ns[i][1];
				}
				else {
					p = i;
					i = ns[i][0];
				}
			}
		}
	}

	void align_edges() {
		int i = 0;
		int j = ns[0][1];
		while (j != 0) {
			if (ns[j][0] != i) {
				ns[j].flip();
			}
			i = j;
			j = ns[j][1];
		}
	}

	// assuming that ns is a cycle

	bool incident(int i, int j) {
		return i == j || i == ns[j][1] || i == ns[j][0];
	}

	void flip_path(int i, int j) {
		while (i != j) {
			ns[i].flip();
			i = ns[i][0];
		}
	}

	bool find_opt2_(int i, int e) {
		int j = ns[ns[i][1]][1];
		while (!incident(j, e)) {
			if (sp.dist(i, ns[i][1]) + sp.dist(j, ns[j][1]) >
				sp.dist(i, j) + sp.dist(ns[i][1], ns[j][1])) {
				int i2 = ns[i][1];
				ns[ns[i][1]][0] = ns[j][1];
				ns[ns[j][1]][0] = ns[i][1];
				ns[i][1] = j;
				ns[j][1] = i;
				flip_path(i2, i);
				return true;
			}
			j = ns[j][1];
		}
		return false;
	}

	void find_opt2() {
		int e = 0;
		int i = ns[0][1];
		while (i != e) {
			if (find_opt2_(i, e)) {
				e = i;
			}
			i = ns[i][1];
		}
		// went full
	}

	void calculate_tsp() {
		sort_edges();
		find_greedy();
		align_edges();
	}

	void print_tsp() {
		// travel path and print to stdout
		cerr << "printing output..." << endl;
		int length = 0;
		int i = 0;
		int j = ns[0][1];
		while (j != 0) {
			length += sp.dist(i, j);
			sp.print_path(i, j);
			if (ns[j][0] == i) {
				i = j;
				j = ns[j][1];
			}
			else {
				i = j;
				j = ns[j][0];
			}
		}
		length += sp.dist(i, j);
		sp.print_path(i, j);
		cout << j << endl;

		cerr << "total length: " << length << endl;
	}

	int length() const {
		int length = 0;
		int i = 0;
		int j = ns[0][1];
		while (j != 0) {
			length += sp.dist(i, j);
			i = j;
			j = ns[j][1];
		}
		length += sp.dist(i, j);
		return length;
	}
};

struct mst_preorder_tsp {
	// local variables
	shortest_paths const& sp;
	// adjacency list
	vector<vector<int>> mst;
	priority_queue<edge, vector<edge>, edge_comp> q;

	mst_preorder_tsp(shortest_paths const& sp) : sp(sp), mst(sp.n) {}

	void queue_children(int v, int prev) {
		mst[v].push_back(prev);
		mst[prev].push_back(v);
		for (int i = 0; i < sp.n; ++i) {
			if (mst[i].size() == 0) {
				q.push(edge(sp.dist(v, i), v, i));
			}
		}
	}

	void find_mst() {
		// prim
		for (int i = 1; i < sp.n; ++i) {
			q.push(edge(sp.dist(0, i), 0, i));
		}
		while (!q.empty()) {
			edge e(q.top());
			q.pop();
			if (mst[e.i].size() == 0) {
				queue_children(e.i, e.j);
			}
			else if (mst[e.j].size() == 0) {
				queue_children(e.j, e.i);
			}
		}
		// mst is minimum spanning tree now
	}

	void calculate_tsp() {
		find_mst();
	}

	int length;
	
	int print_preorder_(int v, int parent) {
		int prev = v;
		for (int i = 0; i < mst[v].size(); ++i) {
			if (mst[v][i] != parent) {
				length += sp.dist(prev, mst[v][i]);
				sp.print_path(prev, mst[v][i]);
				prev = print_preorder_(mst[v][i], v);
			}
		}
		return prev;
	}

	void print_tsp() {
		// print mst in preorder
		length = 0;
		int last = print_preorder_(0, -1);
		length += sp.dist(last, 0);
		sp.print_path(last, 0);
		cout << 0 << endl;

		cerr << "total length: " << length << endl;
	}
};

struct cw_edge_comp {
	bool operator()(edge const& a, edge const& b) {
		return a.w < b.w;
	}
};

struct cw_tsp {

	// local variables
	shortest_paths const& sp;
	priority_queue<edge, vector<edge>, cw_edge_comp> q;
	vector<neighbors> ns;

	cw_tsp(shortest_paths const& sp) : sp(sp), ns(sp.n) {}

	void calculate_tsp() {
		for (int i = 1; i < sp.n; ++i) {
			for (int j = 1; j < i; ++j) {
				q.push(edge(sp.dist(0, i) + sp.dist(0, j) - sp.dist(i, j), i, j));
			}
		}
		int n = 2;
		while (n < sp.n && !q.empty()) {
			edge e(q.top());
			q.pop();
			// check degree
			if (ns[e.i].d == 2 || ns[e.j].d == 2)
				continue;
			// check cycles
			int p = e.j;
			int i = e.i;
			for (;;) {
				// no cycles, add
				if (i == -1) {
					ns[e.i].add(e.j);
					ns[e.j].add(e.i);
					++n;
					break;
				}
				// cycle, discard
				if (i == e.j) {
					break;
				}
				// next
				if (ns[i][0] == p) {
					p = i;
					i = ns[i][1];
				}
				else {
					p = i;
					i = ns[i][0];
				}
			}
		}
		// extract last two
		while (!q.empty()) {
			edge e(q.top());
			q.pop();
			if (ns[e.i].d != 2 && ns[e.j].d != 2) {
				ns[0].add(e.i);
				ns[0].add(e.j);
				ns[e.i].add(0);
				ns[e.j].add(0);
			}
		}
	}

	void print_tsp() {
		// travel path and print to stdout
		cerr << "printing output..." << endl;
		int length = 0;
		int i = 0;
		int j = ns[0][1];
		while (j != 0) {
			length += sp.dist(i, j);
			sp.print_path(i, j);
			if (ns[j][0] == i) {
				i = j;
				j = ns[j][1];
			}
			else {
				i = j;
				j = ns[j][0];
			}
		}
		length += sp.dist(i, j);
		sp.print_path(i, j);
		cout << j << endl;

		cerr << "total length: " << length << endl;
	}
};

struct random_greedy_tsp {

	// local variables
	shortest_paths const& sp;
	typedef multiset<edge, cw_edge_comp> queue;
	queue q;

	vector<neighbors> ns;

	// random generator
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution;

	random_greedy_tsp(shortest_paths const& sp, int variance) : sp(sp), ns(sp.n), generator(1), distribution(0, variance - 1) {}

	void sort_edges() {
		cerr << "sorting edges..." << endl;
		progress p(sp.n);
		for (int i = 0; i < sp.n; ++i) {
			p();
			for (int j = 0; j < i; ++j) {
				q.insert(edge(sp.dist(i, j), i, j));
			}
		}
	}

	void find_random_greedy() {
		cerr << "constructing greedy tsp..." << endl;
		progress p(q.size());
		while (!q.empty()) {
			p();

			// choose rth best
			int r = distribution(generator);
			queue::iterator ei = q.begin();
			for (int i = 0; i < r; ++i) {
				if (++ei == q.end()) {
					ei = q.begin();
				}
			}
			edge e(*ei);
			q.erase(ei);

			// check degree
			if (ns[e.i].d == 2 || ns[e.j].d == 2)
				continue;
			// check cycles smaller than n
			int p = e.j;
			int i = e.i;
			int k = 0;
			for (;;) {
				if (++k == sp.n || i == -1) {
					ns[e.i].add(e.j);
					ns[e.j].add(e.i);
					break;
				}
				if (i == e.j) {
					break;
				}
				// next
				if (ns[i][0] == p) {
					p = i;
					i = ns[i][1];
				}
				else {
					p = i;
					i = ns[i][0];
				}
			}
		}
	}

	void align_edges() {
		int i = 0;
		int j = ns[0][1];
		while (j != 0) {
			if (ns[j][0] != i) {
				ns[j].flip();
			}
			i = j;
			j = ns[j][1];
		}
	}

	// assuming that ns is a cycle

	bool incident(int i, int j) {
		return i == j || i == ns[j][1] || i == ns[j][0];
	}

	void flip_path(int i, int j) {
		while (i != j) {
			ns[i].flip();
			i = ns[i][0];
		}
	}

	bool find_opt2_(int i, int e) {
		int j = ns[ns[i][1]][1];
		while (!incident(j, e)) {
			if (sp.dist(i, ns[i][1]) + sp.dist(j, ns[j][1]) >
				sp.dist(i, j) + sp.dist(ns[i][1], ns[j][1])) {
				int i2 = ns[i][1];
				ns[ns[i][1]][0] = ns[j][1];
				ns[ns[j][1]][0] = ns[i][1];
				ns[i][1] = j;
				ns[j][1] = i;
				flip_path(i2, i);
				return true;
			}
			j = ns[j][1];
		}
		return false;
	}

	void find_opt2() {
		int e = 0;
		int i = ns[0][1];
		while (i != e) {
			if (find_opt2_(i, e)) {
				e = i;
			}
			i = ns[i][1];
		}
		// went full
	}

	void clear_path() {
		for (int i = 0; i < sp.n; ++i) {
			ns[i] = neighbors();
		}
	}

	void calculate_tsp(int repeat) {
		int best = inf;
		vector<neighbors> best_path;
		for (int i = 0; i < repeat; ++i) {
			sort_edges();
			find_random_greedy();
			align_edges();
			find_opt2();
			int l = length();
			cerr << l << endl;
			if (l < best) {
				best = l;
				best_path = ns;
			}
			clear_path();
		}
		ns = best_path;
	}

	void print_tsp() {
		// travel path and print to stdout
		cerr << "printing output..." << endl;
		int length = 0;
		int i = 0;
		int j = ns[0][1];
		while (j != 0) {
			length += sp.dist(i, j);
			sp.print_path(i, j);
			if (ns[j][0] == i) {
				i = j;
				j = ns[j][1];
			}
			else {
				i = j;
				j = ns[j][0];
			}
		}
		length += sp.dist(i, j);
		sp.print_path(i, j);
		cout << j << endl;

		cerr << "total length: " << length << endl;
	}

	int length() const {
		int length = 0;
		int i = 0;
		int j = ns[0][1];
		while (j != 0) {
			length += sp.dist(i, j);
			i = j;
			j = ns[j][1];
		}
		length += sp.dist(i, j);
		return length;
	}
};

void usage()
{
	cerr << "usage: tsp input method [variance repeats]" << endl;
}

int main(int argc, char const *argv[])
{
	if (argc < 3) {
		usage();
		return 0;
	}

	shortest_paths sp;
	if (!sp.load_file(argv[1]) || 
		!sp.calculate_shortest_paths()) {
		return 0;
	}

	if (strcmp(argv[2], "greedy") == 0) {
		greedy_tsp greedy(sp);
		greedy.calculate_tsp();
		greedy.print_tsp();
	}
	else if (strcmp(argv[2], "greedyopt2") == 0) {
		greedy_tsp greedy(sp);
		greedy.calculate_tsp();
		greedy.find_opt2();
		greedy.print_tsp();
	}
	else if (strcmp(argv[2], "rndgreedyopt2") == 0) {
		if (argc < 5) {
			usage();
			return 0;
		}
		random_greedy_tsp greedy(sp, atoi(argv[3]));
		greedy.calculate_tsp(atoi(argv[4]));
		greedy.print_tsp();
	}
	else if (strcmp(argv[2], "mstp") == 0) {
		mst_preorder_tsp mstp(sp);
		mstp.calculate_tsp();
		mstp.print_tsp();
	}
	else if (strcmp(argv[2], "cw") == 0) {
		cw_tsp cw(sp);
		cw.calculate_tsp();
		cw.print_tsp();
	}

	return 0;
}
