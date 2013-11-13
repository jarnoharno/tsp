// Jarno Leppänen 13.11.2013

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <tuple>

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

void usage()
{
	cerr << "usage: tsp input" << endl;
}

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
		if (w >= 0 && w < inf && m[ir + j] > w) {
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

struct tsp_greedy {

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
		int i;
		int j;
		neighbors() : i(-1), j(-1) {}
		int degree() {
			return (i == -1 ? 0 : 1) + (j == -1 ? 0 : 1);
		}
		void add(int k) {
			if (i == -1) {
				i = k;
			}
			else if (j == -1) {
				j = k;
			}
		}
	};

	// local variables
	shortest_paths const& sp;
	priority_queue<edge, vector<edge>, edge_comp> q;
	vector<neighbors> ns;

	tsp_greedy(shortest_paths const& sp) : sp(sp), ns(sp.n) {}

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
			if (ns[e.i].degree() == 2 || ns[e.j].degree() == 2)
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
				if (ns[i].i == p) {
					p = i;
					i = ns[i].j;
				}
				else {
					p = i;
					i = ns[i].i;
				}
			}
		}
	}

	void calculate_tsp() {
		sort_edges();
		find_greedy();
	}

	void print_tsp() {
		// travel path and print to stdout
		cerr << "printing output..." << endl;
		int length = 0;
		int i = 0;
		int j = ns[0].j;
		while (j != 0) {
			length += sp.dist(i, j);
			sp.print_path(i, j);
			if (ns[j].i == i) {
				i = j;
				j = ns[j].j;
			}
			else {
				i = j;
				j = ns[j].i;
			}
		}
		length += sp.dist(i, j);
		sp.print_path(i, j);
		cout << j << endl;

		cerr << "total length: " << length << endl;
	}
};

int main(int argc, char const *argv[])
{
	if (argc < 2) {
		usage();
		return 0;
	}

	shortest_paths sp;
	if (!sp.load_file(argv[1]) || 
		!sp.calculate_shortest_paths()) {
		return 0;
	}

	tsp_greedy greedy(sp);
	greedy.calculate_tsp();
	greedy.print_tsp();

	return 0;
}
