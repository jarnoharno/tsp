// Jarno Leppänen 13.11.2013

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <queue>
#include <tuple>

using namespace std;
static const int inf = numeric_limits<int>::max();

struct sq_matrix {
	vector<int> m;
	int n;
	sq_matrix(int n, int v):
		m(n * n, v),
		n(n)
	{
	}

	int& operator()(int row, int col) {
		return m[row * n + col];
	}

	int const& operator()(int row, int col) const {
		return m[row * n + col];
	}

	void print() {
		int i = 0;
		for (int row = 0; row < n; ++row) {
			for (int col = 0; col < n; ++col) {
				if (m[i] == inf) {
					cerr << '~';
				}
				else {
					cerr << m[i];
				}
				cerr << ' ';
				++i;
			}
			cerr << endl;
		}
	}
};

// floyd warshall path reconstruction
void print_path_(sq_matrix const& dist, sq_matrix const& next, int i, int j) {
	if (dist(i, j) == inf) {
		cerr << "error finding path";
		return;
	}
	int k = next(i, j);
	if (k != -1) {
		print_path_(dist, next, i, k);
		cout << k << endl;
		print_path_(dist, next, k, j);
	}
}

void print_path(sq_matrix const& dist, sq_matrix const& next, int i, int j) {
	cout << i << endl;
	print_path_(dist, next, i, j);
}

struct progress {
	int n;
	int p;
	int i;
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

int main(int argc, char *argv[])
{
	if (argc < 2) {
		usage();
		return 0;
	}

	// read adjacency matrix
	ifstream f(argv[1]);
	if (!f.is_open()) {
		cerr << "can't open file " << argv[1] << argv << endl;
		return 0;
	}
	cerr << "reading graph data from file " << argv[1] << "..." << endl;
	int n, en;
	f >> n;
	f >> en;
	if (f.fail()) {
		cerr << "can't read number of vertices/edges" << endl;
	}
	sq_matrix adj(n, inf);
	progress p;
	p.reset(en);
	for (int i = 0; i < en; ++i) {
		p();
		int u, v, w;
		f >> u >> v >> w;
		adj(u, v) = w;
		adj(v, u) = w;
	}
	p.stop();
	if (f.fail()) {
		cerr << "can't read edge data" << endl;
	}

	// symmetric floyd warshall with path reconstruction
	// optimization ideas: johnson's algorithm, memory optimization
	cerr << "calculating shortest paths..." << endl;
	sq_matrix dist(adj);
	sq_matrix next(n, -1);
	for (int i = 0; i < n; ++i) {
		dist(i, i) = 0;
	}
	p.reset(n);
	for (int k = 0; k < n; ++k) {
		p();
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j <= i; ++j) {
				int w = dist(i, k) + dist(k, j);
				// check for overflow (max + max == -2)
				if (w > 0 && w < inf && dist(i, j) > w) {
					dist(i, j) = dist(j, i) = w;
					next(i, j) = next(j, i) = k;
				}
			}
		}
	}
	p.stop();
	// check negative cycles
	for (int i = 0; i < n; ++i) {
		if (dist(i, i) != 0) {
			cerr << "graph contains negative cycles" << endl;
			return 0;
		}
	}
	// shortest paths graph satisfies triangle equality!

	// greedy
	cerr << "sorting edges..." << endl;
	struct edge {
		int w;
		int i;
		int j;
		edge(int w, int i, int j) : w(w), i(i), j(j) {}
	};
	struct edge_comp {
		bool operator()(edge const& a, edge const& b) {
			return a.w < b.w;
		}
	};
	p.reset(n);
	priority_queue<edge, vector<edge>, edge_comp> q;
	for (int i = 0; i < n; ++i) {
		p();
		for (int j = 0; j < i; ++j) {
			q.push(edge(dist(i, j), i, j));
		}
	}
	p.stop();

	cerr << "constructing greedy tsp..." << endl;
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
	vector<neighbors> ns(n);
	p.reset(q.size());
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
			if (++k == n || i == -1) {
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
			} else {
				p = i;
				i = ns[i].i;
			}
		}
	}
	p.stop();

	// travel path and print to stdout
	cerr << "printing output..." << endl;
	int length = 0;
	int i = 0;
	int j = ns[0].j;
	while (j != 0) {
		length += dist(i, j);
		print_path(dist, next, i, j);
		if (ns[j].i == i) {
			i = j;
			j = ns[j].j;
		}
		else {
			i = j;
			j = ns[j].i;
		}
	}
	length += dist(i, j);
	print_path(dist, next, i, j);
	cout << j << endl;

	cerr << "total length: " << length << endl;

	return 0;
}
