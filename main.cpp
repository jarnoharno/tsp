#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

using namespace std;

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
				if (m[i] == numeric_limits<int>::max()) {
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
	ifstream f(argv[1]);
	if (!f.is_open()) {
		cout << "can't open file " << argv << endl;
		return 0;
	}
	int vertices, edges;
	f >> vertices;
	f >> edges;
	if (f.fail()) {
		cout << "can't read number of vertices/edges" << endl;
	}
	sq_matrix adj_matrix(vertices, numeric_limits<int>::max());
	for (int i = 0; i < edges; ++i) {
		int u, v, w;
		f >> u >> v >> w;
		adj_matrix(u, v) = w;
		adj_matrix(v, u) = w;
	}
	if (f.fail()) {
		cout << "can't read edge data" << endl;
	}
	adj_matrix.print();
	return 0;
}
