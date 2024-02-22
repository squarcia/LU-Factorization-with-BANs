#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <numeric>
#include <algorithm>
#include "ban.h"

typedef unsigned int uint;
using namespace std;
using namespace Eigen;

bool na_simplex(const Matrix<ban, Dynamic, Dynamic> &A, const Vector<ban, Dynamic> &b, const Vector<ban, Dynamic> &c, vector<uint> &B, T tol,
                Vector<ban, Dynamic> const &x, ban &optimal_value);

bool na_simplex_test(const Matrix<ban, Dynamic, Dynamic> &A, const Vector<ban, Dynamic> &b, const Vector<ban, Dynamic> &c, vector<uint> &B, T tol,
                Vector<ban, Dynamic> const &x, ban &optimal_value, bool th);

void modify(Matrix<ban, Dynamic, Dynamic> &A, Vector<ban, Dynamic> &b, Vector<ban, Dynamic> &c, const vector<int> &t, vector<uint> &B);

bool i_big_m(Matrix<ban, Dynamic, Dynamic> A, Vector<ban, Dynamic> b, Vector<ban, Dynamic> c, const vector<int> &t, T tol,
             Vector<ban, Dynamic> const &x, ban &optimal_value);

bool i_big_m_test(Matrix<ban, Dynamic, Dynamic> A, Vector<ban, Dynamic> b, Vector<ban, Dynamic> c, const vector<int> &t, T tol,
             Vector<ban, Dynamic> const &x, ban &optimal_value, bool th);