#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

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

	void print() {
		int i = 0;
		for (int row = 0; row < n; ++row) {
			for (int col = 0; col < n; ++col) {
				if (m[i] == inf) {
					cout << '~';
				}
				else {
					cout << m[i];
				}
				cout << ' ';
				++i;
			}
			cout << endl;
		}
	}
};

struct progress {
	int n;
	int p;
	int i;
	void stop() {
		cout << "\r         \r";
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
		cout << '\r' << p << "% done";
	}
};

void usage()
{
	cout << "usage: tsp input" << endl;
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
		cout << "can't open file " << argv << endl;
		return 0;
	}
	cout << "reading graph data..." << endl;
	int n, en;
	f >> n;
	f >> en;
	if (f.fail()) {
		cout << "can't read number of vertices/edges" << endl;
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
		cout << "can't read edge data" << endl;
	}

	// symmetric floyd warshall with path reconstruction
	// optimization ideas: johnson's algorithm, memory optimization
	cout << "calculating shortest paths..." << endl;
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
	// check negative cycles
	for (int i = 0; i < n; ++i) {
		if (dist(i, i) != 0) {
			cout << "graph contains negative cycles" << endl;
			return 0;
		}
	}
	// shortest paths graph satisfies triangle equality!
	next.print();

	return 0;
}
