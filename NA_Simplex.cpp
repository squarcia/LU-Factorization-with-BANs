#include "NA_Simplex.h"
#include <Eigen/Eigenvalues>
#define EALL Eigen::placeholders::all

void debug(const Vector<ban,Dynamic> &A, const char* s){
    cout << s << ": ";
    for(auto x : A)
        cout << x << " ";
    cout << endl << endl;
}

void debug(const vector<uint> &v, const char* s){
    cout << s << ": ";
    for(auto x : v)
        cout << x << " ";
    cout << endl << endl;
}

//print a squared matrix of bans
void debugFact(Matrix<ban, Dynamic, Dynamic> &LU, int dim){
    cout << "LU:" << endl;
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            cout << "LU(" << i << "," << j << "): " << LU(i,j); 
        }
        cout << endl;
    }
}

//print a non-squared matrix of bans
void debugMat(const Matrix<ban, Dynamic, Dynamic> &matrix, int rows, int columns){
    cout << "matrix:" << endl;
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            cout << "matrix(" << i << "," << j << "): " << matrix(i,j); 
        }
        cout << endl;
    }
}

/// @brief compute matrix rank given the finite components of the matrix
/// @param A input matrix
/// @param n dimension of the matrix
void computeRank(const Matrix<ban,Dynamic,Dynamic> &A, int n){
    // Calcola la decomposizione LU
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix(i, j) = A(i,j).lead_mon();
        }
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    int rank = svd.rank();
    cout << "Matrix dim: " << n << endl;
    cout << "Matrix rank: " << rank << endl;
}

/// @brief compute the condition number of the matrix given the finite components of the matrix
/// @param LU_ input matrix
/// @param n dimension of the matrix
void computeCondNumber(const Matrix<ban,Dynamic,Dynamic> &LU_, int n){
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix(i, j) = LU_(i,j).lead_mon();
        }
    }
    //cout << "matrix: " << endl << matrix << endl;
    double det = matrix.determinant();
    cout << "Determinant: " << det << endl;
    Eigen::MatrixXd A_inverse = matrix.transpose() * matrix; // Calcola A^T * A
    double detInv = A_inverse.determinant();
    cout << "Inverse determinant: " << detInv << endl;
    //cout << "inverse: " << endl << A_inverse << endl;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A_inverse, Eigen::ComputeFullU | Eigen::ComputeFullV);
    double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
    std::cout << "Condition number: " << cond << endl;
}

/// @brief check if input matrix is singular given the finite components of the matrix
/// @param LU input matrix
/// @param n dimension of the matrix
void checkSingularity(const Matrix<ban,Dynamic,Dynamic> &LU, int n){
    // Calcola la decomposizione LU
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix(i, j) = LU(i,j).lead_mon();
        }
    }
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
    //cout << lu_decomp.matrixLU() << endl;
    // Ottieni il rango della matrice
    int rank = lu_decomp.rank();
    //cout << "matrix: " << endl << matrix << endl;
    // Determina se la matrice è prossima alla singolarità
    bool isNearSingularity = rank < std::min(LU.rows(), LU.cols());
    if (isNearSingularity) {
        std::cout << "Matrix is close to singularity" << std::endl;
    } else {
        std::cout << "Matrix is not close to singularity." << std::endl;
    }
}

/// @brief check if matrix is lu factorizable
/// @param A input matrix
/// @param n dimension of the matrix
void checkFactorizable(const Matrix<ban,Dynamic,Dynamic> &A, int n){
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix(i, j) = A(i,j).lead_mon();
        }
    }
    //cout << "Matrix: " << endl << matrix << endl;
    for(int i = 0; i < n; i++) {
        int dim = i + 1;
        MatrixXd minMat = matrix.topLeftCorner(dim, dim);
        //cout << "MinMat: " << endl << minMat << endl;
        int max_index = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(matrix(k, i)) > abs(matrix(max_index, i))) {
                max_index = k;
            }
        }
        matrix.row(max_index).swap(matrix.row(i));
        double det = minMat.determinant();
        //cout << "MinMat determinant: " << det << endl;
        if(det == 0) {
            cout << "Matrix not factorizable" << endl;
            return;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

void LU_fact_tot(const Matrix<ban, Dynamic, Dynamic> &A, Matrix<ban, Dynamic, Dynamic> const &LU, vector<uint> &pivot, Vector<ban, Dynamic> &b){
    T tol = 1e-4;
    const uint n_constraints = A.rows();
    cout << "TEST A" << endl;
    checkFactorizable(A,n_constraints);
    computeCondNumber(A,n_constraints);
    checkSingularity(A,n_constraints);
    computeRank(A,n_constraints);
    cout << endl;
    Matrix<ban, Dynamic, Dynamic> LU_ = A; 
    debugFact(LU_, n_constraints);
    iota(pivot.begin(), pivot.end(), 0);
    for(auto k = 0; k < n_constraints; ++k){
        uint kp = k;
        uint kj = k;
        ban amax = abs(LU_(k,k));
        ban amin = abs(LU_(k,k));
        for(auto j = k+1; j < n_constraints; ++j){
            for(auto i = k+1; i < n_constraints; ++i){
                ban absi = abs(LU_(i,j));
                if(absi > amax) {
                    kp = i;
                    kj = j;
                    amax = absi;
                }
                if(absi > 0 && absi < amin)
                    amin = absi;
            }
        }
        cout << amax;
        //if(amax.lead_mon() == 0) continue;
        //compute the threshold as the division of the monosemnia of the minimum value and the monosemnia of the maximum of the column
        //auto th = amin/amax;
        //cout << "Threshold: " << th << endl;
        cout << kp <<  "," << kj << "= " << LU_(kp,kj) <<  endl;
        if(LU_(kp,kj) != 0){
            //cout << "elem not zero: " << kp <<  "," << k << "= " << LU_(kp,k) <<  endl;
            //cout << kp << endl;
            if((kp != k) && (kj!=k)){
                auto tmp_pivot = pivot[k];
                pivot[k] = pivot[kp];
                pivot[kp] = tmp_pivot;
                //debug(pivot, "pivot");

                auto tmp_b = b[k];
                b[k] = b[kj];
                b[kp] = tmp_b;

                for(auto i = 0; i < n_constraints; ++i){
                    auto tmp = LU_(k,i);
                    LU_(k,i) = LU_(kp, i);
                    LU_(kp,i) = tmp;
                }

                for(auto j = 0; j < n_constraints; ++j){
                    auto tmp = LU_(j,k);
                    LU_(j,k) = LU_(j,kj);
                    LU_(j,kj) = tmp;
                }
            }
            
            //cout << "Akkinv: 1/" << LU_(k,k) << endl;
            auto Akkinv = 1/LU_(k,k);
            //cout << "Akkinv: " << Akkinv << endl;
            for(auto i = k+1; i < n_constraints; ++i) {
                LU_(i, k) *= Akkinv;
            }
        }
        else{
            throw domain_error("pivot equal to zero detected");
        }
        for(auto j = k+1; j < n_constraints; ++j){
            for(auto i = k+1; i < n_constraints; ++i){
                LU_(i,j) -= LU_(i,k) * LU_(k,j);
                LU_(i,j).denoise(tol);
            }
        }
        debugFact(LU_, n_constraints);
    }
    //cout << "TEST LU" << endl;
    //computeCondNumber(LU_, n_constraints);
    //checkSingularity(LU_, n_constraints);
    //computeRank(LU_, n_constraints);
    debugFact(LU_, n_constraints);
    debug(pivot, "pivot");
    (const_cast< Matrix<ban, Dynamic, Dynamic> &>(LU)) = LU_;
}

void LU_fact(const Matrix<ban, Dynamic, Dynamic> &A, Matrix<ban, Dynamic, Dynamic> const &LU, vector<uint> &pivot){
    T tol = 1e-4;
    const uint n_constraints = A.rows();
    cout << "TEST A" << endl;
    checkFactorizable(A,n_constraints);
    computeCondNumber(A,n_constraints);
    checkSingularity(A,n_constraints);
    computeRank(A,n_constraints);
    cout << endl;
    Matrix<ban, Dynamic, Dynamic> LU_ = A; 
    iota(pivot.begin(), pivot.end(), 0);
    for(auto k = 0; k < n_constraints; ++k){
        uint kp = k;
        ban amax = abs(LU_(k,k));
        ban amin = abs(LU_(k,k));
        for(auto i = k+1; i < n_constraints; ++i){
            ban absi = abs(LU_(i,k));
            if(absi > amax) {
                kp = i;
                amax = absi;
            }
            if(absi > 0 && absi < amin)
                amin = absi;
        }
        //compute the threshold as the division of the monosemnia of the minimum value and the monosemnia of the maximum of the column
        auto th = amin/amax;
        //cout << "Threshold: " << th << endl;
        if(LU_(kp,k) != 0){
            //cout << "elem not zero: " << kp <<  "," << k << "= " << LU_(kp,k) <<  endl;
            //cout << kp << endl;
            if(kp != k){
                auto tmp_pivot = pivot[k];
                pivot[k] = pivot[kp];
                pivot[kp] = tmp_pivot;
                //debug(pivot, "pivot");

                for(auto i = 0; i < n_constraints; ++i){
                    auto tmp = LU_(k,i);
                    LU_(k,i) = LU_(kp, i);
                    LU_(kp,i) = tmp;
                }
            }
            auto Akkinv = 1/LU_(k,k);
            //cout << "Akkinv: " << Akkinv << endl;
            for(auto i = k+1; i < n_constraints; ++i) {
                LU_(i, k) *= Akkinv;
            }
        }
        else{
            throw domain_error("pivot equal to zero detected");
        }
        for(auto j = k+1; j < n_constraints; ++j){
            for(auto i = k+1; i < n_constraints; ++i){
                LU_(i,j) -= LU_(i,k) * LU_(k,j);
                LU_(i,j).denoise(tol);
            }
        }
        //debugFact(LU_, n_constraints);
    }
    //cout << "TEST LU" << endl;
    //computeCondNumber(LU_, n_constraints);
    //checkSingularity(LU_, n_constraints);
    //computeRank(LU_, n_constraints);
    debugFact(LU_, n_constraints);
    debug(pivot, "pivot");
    (const_cast< Matrix<ban, Dynamic, Dynamic> &>(LU)) = LU_;
}

void LU_fact_th(const Matrix<ban, Dynamic, Dynamic> &A, Matrix<ban, Dynamic, Dynamic> const &LU, vector<uint> &pivot){
    T tol = 1e-4;
    const uint n_constraints = A.rows();
    cout << "TEST A" << endl;
    checkFactorizable(A,n_constraints);
    computeCondNumber(A,n_constraints);
    checkSingularity(A,n_constraints);
    computeRank(A,n_constraints);
    cout << endl;
    Matrix<ban, Dynamic, Dynamic> LU_ = A; 
    iota(pivot.begin(), pivot.end(), 0);
    for(auto k = 0; k < n_constraints; ++k){
        uint kp = k;
        ban amax = abs(LU_(k,k));
        ban amin = abs(LU_(k,k));
        for(auto i = k+1; i < n_constraints; ++i){
            ban absi = abs(LU_(i,k));
            if(absi > amax) {
                kp = i;
                amax = absi;
            }
            if(absi > 0 && absi < amin)
                amin = absi;
        }
        //cout << "MAX: " << amax << endl;
        //cout << "MIN: " << amin << endl;
        //compute the threshold as the division of the monosemnia of the minimum value and the monosemnia of the maximum of the column
        auto th = amin;
        if(amax == 0){
            th = 0;
            LU_(k,k) = ETA;
            //continue;
        }
        else{
            th = amin/amax;
        }
        cout << "Threshold: " << th << endl;
        if(LU_(kp,k) != 0){
            //cout << "elem not zero: " << kp <<  "," << k << "= " << LU_(kp,k) <<  endl;
            //cout << kp << endl;
            if(kp != k){
                auto tmp_pivot = pivot[k];
                pivot[k] = pivot[kp];
                pivot[kp] = tmp_pivot;
                //debug(pivot, "pivot");

                for(auto i = 0; i < n_constraints; ++i){
                    auto tmp = LU_(k,i);
                    LU_(k,i) = LU_(kp, i);
                    LU_(kp,i) = tmp;
                }
            }
            auto Akkinv = 1/LU_(k,k);
            /*ban Akkinv = 1;
            if(LU_(k,k).lead_mon()== 0) Akkinv = ALPHA;
            else
                Akkinv = 1/LU_(k,k);
            */
            cout << "Akkinv: " << Akkinv << endl;
            //Check if the value for which we divide the column monosemnia is greater or equal the threshold;
            //if not, don't scale the coefficients
            if(abs(Akkinv)>=th){
                for(auto i = k+1; i < n_constraints; ++i) {
                    LU_(i, k) *= Akkinv;
                }
            }
        }
        else{
            throw domain_error("pivot equal to zero detected");
        }
        for(auto j = k+1; j < n_constraints; ++j){
            for(auto i = k+1; i < n_constraints; ++i){
                LU_(i,j) -= LU_(i,k) * LU_(k,j);
                //new control
                /*if(LU_(i,j) < tol){
                    LU_(i,j).denoise(tol);
                }*/
                LU_(i,j).denoise(tol);
            }
        }
        //debugFact(LU_, n_constraints);
    }
    //cout << "TEST LU" << endl;
    //computeCondNumber(LU_, n_constraints);
    //checkSingularity(LU_, n_constraints);
    //computeRank(LU_, n_constraints);
    //debugFact(LU_, n_constraints);
    
    debugFact(LU_, n_constraints);
    debug(pivot, "pivot");
    (const_cast< Matrix<ban, Dynamic, Dynamic> &>(LU)) = LU_;
}

void LU_solve(const Matrix<ban, Dynamic, Dynamic> &LU, const Vector<ban, Dynamic> &b, const vector<uint> &pivot, Vector<ban, Dynamic> const &x){
    T tol = 1e-6; //tolerance for testing how much reducing the precision will impact the solution
    const uint dim = LU.cols();
    Vector<ban, Dynamic> xx(dim);
    Vector<ban, Dynamic> yy(dim);
    //the vectors are randomly initialized, so we force them to zero
    yy.setZero();
    xx.setZero();

    // solve for L
    for(auto i = 0; i < dim; ++i){
        auto bi = b(pivot[i]);
        //cout<< "bi: " << bi;
        for(auto j = 0; j < i; ++j){
            //cout << "factor: LU(" <<i<<","<<j<<")=" << LU(i,j) << "*" << yy(j);
            bi -= LU(i,j)*yy(j);
            //cout << "bi: " << bi;
        }
        yy(i) = bi;
        //cout << "yi: " << yy(i);
        //yy(i).denoise(tol);
    }
    for(auto& yyi : yy){
        yyi.denoise(tol);
    }
    //debug(yy, "y");
    // solve for U
    for(int i = dim-1; i >= 0; --i){
        auto yi = yy(i);
        //cout << "yi: " << yi;
        for(int j = dim-1; j > i; --j){
            //cout << "factor: LU(" <<i<<","<<j<<")=" << LU(i,j) << "*" << xx(j);
            yi -= LU(i,j)*xx(j);
            //cout << "yi: " << yi;
        }
        //cout << "yi after: " << yi;
        //cout << "LU(" << i << "," << i << ")=" << LU(i,i);
        xx(i) = yi/LU(i,i);
        //xx(i).denoise(tol);
        //cout << "xi: " << xx(i) << endl;
        //debug(xx, "xx");
    }
    for(auto& xxi : xx){
        xxi.denoise(tol);
    }
    debug(xx, "x");

    (const_cast< Vector<ban, Dynamic> &>(x)) = xx;
}


void LU_trans_solve(const Matrix<ban, Dynamic, Dynamic> &LU, const Vector<ban, Dynamic> &b, const vector<uint> &pivot, Vector<ban, Dynamic> const &x){
    T tol = 1e-6; //tolerance for testing how much reducing the precision will impaxt the solution
    const uint dim = LU.cols();
    Vector<ban, Dynamic> xx(dim);
    Vector<ban, Dynamic> yy(dim);
    //the vectors are randomly initialized, so we force them to zero
    yy.setZero();
    xx.setZero();
    
    // solve for U^T
    for(auto i = 0; i < dim; ++i){
        auto bi = b(i);
        //cout << "bi: " << bi;
        for(auto j = 0; j < i; ++j){
            //cout << "factor: LU(" <<j<<","<<i<<")=" << LU(j,i) << "*" << yy(j);
            bi -= LU(j,i)*yy(j);
            //cout << "bi: " << bi;
        }
        //cout << "bi after: " << bi;
        //cout << "LU(" << i << "," << i << ")=" << LU(i,i);
        yy(i) = bi/LU(i,i);
        //yy(i).denoise(tol);
        //cout << "yi: " << yy(i) << endl;
    }
    for(auto& yyi : yy){
        yyi.denoise(tol);
    }

    // solve for L^T
    for(int i = dim-1; i >= 0; --i){
        auto yi = yy(i);
        //cout << "yi: " << yi;
        for(int j = dim-1; j > i; --j){
            //cout << "factor: LU(" <<j<<","<<i<<")=" << LU(j,i) << "*" << xx(j);
            yi -= LU(j,i)*xx(j);
            //cout << "yi: " << yi;
        }
        xx(i) = yi;
        //cout << "xi: " << xx(i)<<endl;
        //xx(i).denoise(tol);
        //debug(xx, "xx");
    }
    for(auto& xxi : xx){
        xxi.denoise(tol);
    }
    debug(xx, "x_transpose");

    xx(pivot) = xx;
    (const_cast< Vector<ban, Dynamic> &>(x)) = xx;
}

/* The method implements the revised simplex method in Box 7.1 on page 103 of Chvatal

	Revised Simplex
	
	max  c'*x
	s.t. Ax = b
	      x >= 0

*/

bool na_simplex_test(const Matrix<ban, Dynamic, Dynamic> &A, const Vector<ban, Dynamic> &b, const Vector<ban, Dynamic> &c, vector<uint> &B, T tol,
                Vector<ban, Dynamic> const &x, ban &optimal_value, bool th){

	// TODO control dimensions coherence
    const uint n_constraints = A.rows();
    const uint n_variables = A.cols();
    if(n_constraints != b.size() || n_variables != c.size()){
        cout << "Dimensions not coherent" << endl;
        return false;
    }
    //debugMat(A, n_constraints, n_variables);
    // Assume rank non-deficient initial base matrix
    vector<uint> variables;
    variables.reserve(n_variables);
    for(auto i = 0; i < n_variables; ++i)
        variables.push_back(i);

    vector<uint> N;
    N.reserve(n_variables-n_constraints);
    set_difference(variables.begin(), variables.end(), B.begin(), B.end(), inserter(N, N.begin()));
    debug(B,"B");
    debug(N,"N");
    vector<uint> pivot(n_constraints);
    Matrix<ban, Dynamic, Dynamic> LU(n_constraints, n_constraints);
    if(th)
        LU_fact_th(A(EALL, B), LU, pivot);
    else
        LU_fact(A(EALL, B), LU, pivot);
    Vector<ban, Dynamic> xB(n_constraints);
    LU_solve(LU, b, pivot, xB);
    for(auto& xBi : xB) {
        xBi.denoise(tol);
    }
    debug(xB, "xB");
    Vector<ban, Dynamic> xx = Vector<ban, Dynamic>::Zero(n_variables);
    xx(B) = xB;

    Vector<ban, Dynamic> y(n_constraints), sN(n_variables-n_constraints), d(n_constraints), quality(n_constraints);
    int k;
    uint l;
    ban k_val;
    uint tmp;
    vector<uint> zz;
    bool one_positive;
    Index ii;
    int steps = 0;
    while(true){
        steps++;
        if(steps==20) break;
        cout << "################################################################STEP " << steps << "################################################################" << endl << endl;
        LU_trans_solve(LU, c(B), pivot, y);
        sN = c(N) - A(EALL, N).transpose() * y;
        //debug(A(EALL, N).transpose() * y, "diff");
        //debug(c,"c");
        //exit(1);
        debug(sN,"sN");

        // entering index
        k = -1;
        for(auto i=0; i < n_variables - n_constraints; ++i){
            k_val = sN(i).lead_mon();
            if(k_val > tol) {
                k = i;
                break;
            }
        }
        cout<<"k: " << k << endl;

        if(k == -1) {
            // solution found
            (const_cast< Vector<ban, Dynamic> &>(x)) = xx;
            debug(x, "Optimal point");
            optimal_value = (c.transpose() * x)(0);
            cout << "Solution found in " << steps << " steps"<<endl;
            return true;
        }
        cout << "Find d"<<endl;
        LU_solve(LU, A(EALL, N[k]), pivot, d);
        debug(d,"d");

        zz.clear();
        one_positive = false;
        for(int i = 0; i < n_constraints; ++i)
            if (d(i) > 0){
                zz.push_back(i);
                one_positive = true;
            }

        if(!one_positive){
            // problem unbounded
            cout << "Problem unbounded" << endl;
            return false;
        }
        quality = (xB(zz)).array() / (d(zz)).array();
        debug(quality, "quality");
        quality.minCoeff(& ii);

        debug(zz, "zz");
        //cout << "i: " << ii<< endl;
        //cout << "zz[i]: " << zz[ii]<< endl;
        
        //Efficient but unstable update when using bans
        /*
        auto t = quality(ii);
        xB -= t * d;
        xx(N[k]) = t;
        */

        // more stable update
        cout << "Update base" << endl;
        l = zz[ii];
        tmp = B[l];
        B[l] = N[k];
        debug(B,"B");
        N[k] = tmp;
        debug(N,"N");
        cout << "Find new base solution" << endl;
        if(th)
            LU_fact_th(A(EALL, B), LU, pivot);
        else
            LU_fact(A(EALL, B), LU, pivot);
        LU_solve(LU, b, pivot, xB);
        xx(N[k]) = 0;
        for(auto& elem : xB) {
            elem.denoise(tol);
        }
        debug(xB,"xB");
        xx(B) = xB;
    }

}

bool na_simplex(const Matrix<ban, Dynamic, Dynamic> &A, const Vector<ban, Dynamic> &b, const Vector<ban, Dynamic> &c, vector<uint> &B, T tol,
                Vector<ban, Dynamic> const &x, ban &optimal_value){

	// TODO control dimensions coherence
    const uint n_constraints = A.rows();
    const uint n_variables = A.cols();
    if(n_constraints != b.size() || n_variables != c.size()){
        cout << "Dimensions not coherent" << endl;
        return false;
    }
    //debugMat(A, n_constraints, n_variables);
    // Assume rank non-deficient initial base matrix
    vector<uint> variables;
    variables.reserve(n_variables);
    for(auto i = 0; i < n_variables; ++i)
        variables.push_back(i);

    vector<uint> N;
    N.reserve(n_variables-n_constraints);
    set_difference(variables.begin(), variables.end(), B.begin(), B.end(), inserter(N, N.begin()));
    debug(B,"B");
    debug(N,"N");
    vector<uint> pivot(n_constraints);
    Matrix<ban, Dynamic, Dynamic> LU(n_constraints, n_constraints);
    //LU_fact(A(EALL, B), LU, pivot);
    Vector<ban, Dynamic> b_tmp = b;
    LU_fact_tot(A(EALL, B), LU, pivot, b_tmp);
    Vector<ban, Dynamic> xB(n_constraints);
    LU_solve(LU, b_tmp, pivot, xB);
    for(auto& xBi : xB) {
        xBi.denoise(tol);
    }
    debug(xB, "xB");
    Vector<ban, Dynamic> xx = Vector<ban, Dynamic>::Zero(n_variables);
    xx(B) = xB;

    Vector<ban, Dynamic> y(n_constraints), sN(n_variables-n_constraints), d(n_constraints), quality(n_constraints);
    int k;
    uint l;
    ban k_val;
    uint tmp;
    vector<uint> zz;
    bool one_positive;
    Index ii;
    int steps = 0;
    while(true){
        steps++;
        cout << "################################################################STEP " << steps << "################################################################" << endl << endl;
        LU_trans_solve(LU, c(B), pivot, y);
        sN = c(N) - A(EALL, N).transpose() * y;
        //debug(A(EALL, N).transpose() * y, "diff");
        //debug(c,"c");
        //exit(1);
        debug(sN,"sN");

        // entering index
        k = -1;
        for(auto i=0; i < n_variables - n_constraints; ++i){
            k_val = sN(i).lead_mon();
            if(k_val > tol) {
                k = i;
                break;
            }
        }
        cout<<"k: " << k << endl;

        if(k == -1) {
            // solution found
            (const_cast< Vector<ban, Dynamic> &>(x)) = xx;
            debug(x, "Optimal point");
            optimal_value = (c.transpose() * x)(0);
            cout << "Solution found in " << steps << " steps"<<endl;
            return true;
        }
        cout << "Find d"<<endl;
        LU_solve(LU, A(EALL, N[k]), pivot, d);
        debug(d,"d");

        zz.clear();
        one_positive = false;
        for(int i = 0; i < n_constraints; ++i)
            if (d(i) > 0){
                zz.push_back(i);
                one_positive = true;
            }

        if(!one_positive){
            // problem unbounded
            cout << "Problem unbounded" << endl;
            return false;
        }
        quality = (xB(zz)).array() / (d(zz)).array();
        debug(quality, "quality");
        quality.minCoeff(& ii);

        debug(zz, "zz");
        //cout << "i: " << ii<< endl;
        //cout << "zz[i]: " << zz[ii]<< endl;
        
        //Efficient but unstable update when using bans
        /*
        auto t = quality(ii);
        xB -= t * d;
        xx(N[k]) = t;
        */

        // more stable update
        cout << "Update base" << endl;
        l = zz[ii];
        tmp = B[l];
        B[l] = N[k];
        debug(B,"B");
        N[k] = tmp;
        debug(N,"N");
        cout << "Find new base solution" << endl;
        //LU_fact(A(EALL, B), LU, pivot);
        Vector<ban, Dynamic> b_tmp = b;
        LU_fact_tot(A(EALL, B), LU, pivot, b_tmp);
        LU_solve(LU, b_tmp, pivot, xB);
        xx(N[k]) = 0;
        for(auto& elem : xB) {
            elem.denoise(tol);
        }
        debug(xB,"xB");
        xx(B) = xB;
    }

}

/*  problem form

    max c^T x
    s.t. Ax <= b if t < 0
		 Ax  = b if t = 0
		 Ax >= b if t > 0
		  x >= 0

    Assume b >= 0
*/

void modify(Matrix<ban, Dynamic, Dynamic> &A, Vector<ban, Dynamic> &b, Vector<ban, Dynamic> &c, const vector<int> &t, vector<uint> &B){

    uint n_constraints = A.rows();
    uint n_variables = A.cols();

    auto AA = A;
    auto bb = b;
    auto cc = c;

    vector<int> idx_smaller;
    vector<int> idx_equal;
    vector<int> idx_larger;

    for(auto i = 0; i < t.size(); ++i){
        if(t[i] < 0)
            idx_smaller.push_back(i);
        else if(t[i] == 0)
            idx_equal.push_back(i);
        else
            idx_larger.push_back(i);
    }

    uint n_less = idx_smaller.size();
    uint n_equal = idx_equal.size();
    uint n_larger = idx_larger.size();

    /* Make A as follows
             _					   _
			|  Ale	I	0	0	0	|
	   A = 	|  Aeq	0	I	0	0	|
			|_ Age  0	0	I  -I  _|
	*/


    A = Matrix<ban, Dynamic, Dynamic>::Zero(n_constraints, n_variables + n_less + n_equal + 2 * n_larger);
    A(seq(0, n_less - 1), seq(0, n_variables - 1)) = AA(idx_smaller, EALL);
    A(seq(0, n_less - 1), seq(n_variables, n_variables + n_less - 1)) = Matrix<ban, Dynamic, Dynamic>::Identity(n_less, n_less);
    A(seq(n_less, n_less + n_equal -1), seq(0, n_variables - 1)) = AA(idx_equal, EALL);
    A(seq(n_less, n_less + n_equal -1), seq(n_variables + n_less, n_variables + n_less + n_equal - 1)) = Matrix<ban, Dynamic, Dynamic>::Identity(n_equal, n_equal);
    A(seq(n_less + n_equal, n_less + n_equal + n_larger -1), seq(0, n_variables - 1)) = AA(idx_larger, EALL);
    A(seq(n_less + n_equal, n_less + n_equal + n_larger -1), seq(n_variables + n_less + n_equal, n_variables + n_less + n_equal + n_larger - 1)) =  Matrix<ban, Dynamic, Dynamic>::Identity(n_larger, n_larger);
    A(seq(n_less + n_equal, n_less + n_equal + n_larger -1), seq(n_variables + n_less + n_equal + n_larger, n_variables + n_less + n_equal + 2 * n_larger - 1)) = - Matrix<ban, Dynamic, Dynamic>::Identity(n_larger, n_larger);

    /* Make b as follows

        b =  |_ ble	beq	 bge_|^T

     */

    b(seq(0, n_less - 1)) = bb(idx_smaller);
    b(seq(n_less, n_less + n_equal - 1)) = bb(idx_equal);
    b(seq(n_less + n_equal, n_less + n_equal + n_larger -1)) = bb(idx_larger);


    /* Make c as follows

        c =  | _c 0 -α -α 0 |^T

     */

    c = Vector<ban, Dynamic>::Zero(n_variables + n_less + n_equal + 2 * n_larger);
    c(seq(0, n_variables-1)) = cc;
    //c(seq(n_variables + n_less, n_variables + n_less + n_equal + n_larger - 1), 0) = ALPHA * Matrix<T, Dynamic, 1>::Ones(n_equal + n_larger);
    for(auto &elem : c(seq(n_variables + n_less, n_variables + n_less + n_equal + n_larger - 1)))
        elem = -ALPHA;
    //debug(c, "c");
    B.reserve(n_constraints);
    for(auto i = n_variables; i < n_constraints + n_variables; ++i)
        B.push_back(i);
}

bool i_big_m(Matrix<ban, Dynamic, Dynamic> A, Vector<ban, Dynamic> b, Vector<ban, Dynamic> c, const vector<int> &t, T tol,
             Vector<ban, Dynamic> const &x, ban &optimal_value){

    // TODO check positivity of b
    
    cout << "b is positive" << endl;
    // TODO check non-infinity of c
     for(int i = 0; i < c.size(); i++){
        float v = c(i).lead_mon();
        //cout << v << endl;
        if (isinf(v)){
            cout << "c is infinite" << endl;
            return false;
        }
    }
    cout << "c is finite" << endl;

    vector<uint> B;
    modify(A, b, c, t, B);

    return na_simplex(A, b, c, B, tol, x, optimal_value);
}

bool i_big_m_test(Matrix<ban, Dynamic, Dynamic> A, Vector<ban, Dynamic> b, Vector<ban, Dynamic> c, const vector<int> &t, T tol,
             Vector<ban, Dynamic> const &x, ban &optimal_value, bool th){

    // TODO check positivity of b
    
    cout << "b is positive" << endl;
    // TODO check non-infinity of c
     for(int i = 0; i < c.size(); i++){
        float v = c(i).lead_mon();
        //cout << v << endl;
        if (isinf(v)){
            cout << "c is infinite" << endl;
            return false;
        }
    }
    cout << "c is finite" << endl;

    vector<uint> B;
    modify(A, b, c, t, B);

    return na_simplex_test(A, b, c, B, tol, x, optimal_value, th);
}

