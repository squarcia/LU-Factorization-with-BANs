#include "NA_Simplex.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

//template <typename Derived>
//template class Eigen::internal::traits<ban>;
void debug1(const Vector<ban,Dynamic> &A, const char* s){
    cout << s << ": ";
    for(auto x : A)
        cout << x << " ";
    cout << endl << endl;
}
int main(int argc, char *argv[]) {
    while(true){
        cout <<"RUN TEST" << endl;
        cout <<"Test 0: exit" << endl;
        cout <<"Test 1: standard" << endl;
        cout <<"Test 2: i-big-m standard" << endl;
        cout <<"Test 3: na standard" << endl;
        cout <<"Test 4: na i-big-m" << endl;
        cout <<"Test 5: full na" << endl;
        cout <<"Test 6: full na i-big-m (c=0)" << endl;
        cout <<"Test 7: full na i-big-m" << endl;
        cout <<"Test 8: full na with infinite coefficients" << endl;
        cout <<"Test 9: full with na infinitesimal coefficients" << endl;
        cout <<"Test 10: full na infinite eigenvalues" << endl;
        cout <<"Test 13: full na with singular matrix" << endl;
        cout <<"Test 11: standard (test 5 without BANs)" << endl;
        cout <<"Test 12: i-big-m standard (test 6 without BANs)" << endl;
        int info = 0;
        cout << "Insert test to run" << endl;
        cin >> info;
        if((info < 0) || (info > 20)) continue;
        if(info == 0) break;
        bool th = 0;
        cout<<"Select method: threshold active:1, non active:0" << endl;
        cin >> th;
        if(th!= 0 && th != 1) continue;

        ban optimal_value;
        T tol = 1e-4;
        bool flag;

        if(info == 1){
            // standard
            Matrix<ban, 6, 9> A;
            Matrix<ban, 6, 1> b;
            Matrix<ban, 9, 1> c;

            A <<  2,  1, -3, 1, 0, 0, 0, 0, 0,
                2,  3, -2, 0, 1, 0, 0, 0, 0,
                4,  3,  3, 0, 0, 1, 0, 0, 0,
                0,  0,  1, 0, 0, 0, 1, 0, 0,
                0,  0, -1, 0, 0, 0, 0, 1, 0,
                -1, -2, -1, 0, 0, 0, 0, 0, 1;

            b << 90, 190, 300, 10, -10, -70;

            c << 8, 12, 7, 0, 0, 0, 0, 0, 0;

            vector<uint> B = { 1, 2, 3, 4, 5, 7};
            Vector<ban, 9> x;

            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value,th);
            cout << "Standard Kite" << endl;
            //debug1(x, "x");
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 2){
            // i-big-m standard
            Matrix<ban, 5, 3> A;
            Matrix<ban, 5, 1> b;
            Matrix<ban, 3, 1> c;

            A << 2, 1, -3,
                2, 3, -2,
                4, 3,  3,
                0, 0,  1,
                1, 2,  1;

            b << 90, 190, 300, 10, 70;

            c << 8, 12, 7;

            vector<int> t = {-1, -1, -1, 0, 1};

            Matrix<ban, 5, 1> x;

            flag = i_big_m_test(A, b, c, t, tol, x, optimal_value,th);

            //cout<<x<<endl;
            //cout<< A(0) << endl;
            //cout<< b[0] << endl;
            //cout<< c[0] << endl;
            cout << "Standard I-Big-M" << endl;
            cout<<optimal_value<<endl;
            cout<<flag<<endl;
        }

        if(info == 3){
            // na
            Matrix<ban, 6, 9> A;
            Matrix<ban, 6, 1> b;
            Matrix<ban, 9, 1> c;

            A <<  2,  1, -3, 1, 0, 0, 0, 0, 0,
                2,  3, -2, 0, 1, 0, 0, 0, 0,
                4,  3,  3, 0, 0, 1, 0, 0, 0,
                0,  0,  1, 0, 0, 0, 1, 0, 0,
                0,  0, -1, 0, 0, 0, 0, 1, 0,
                -1, -2, -1, 0, 0, 0, 0, 0, 1;

            b << 90, 190, 300, 10, -10, -70;

            c << 8 + 14 * ETA, 12 +  10 * ETA, 7 + 2 * ETA, 0, 0, 0, 0, 0, 0;

            vector<uint> B = {1, 2, 3, 4, 5, 7};
            Vector<ban, 9> x;

            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value, th);

            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 4){
            // i-big-m na
            Matrix<ban, 5, 3> A;
            Matrix<ban, 5, 1> b;
            Matrix<ban, 3, 1> c;

            A << 2, 1, -3,
                2, 3, -2,
                4, 3,  3,
                0, 0,  1,
                1, 2,  1;

            b << 90, 190, 300, 10, 70;

            c << 8 + 14 * ETA, 12 +  10 * ETA, 7 + 2 * ETA;
            //debug1(c,"c: ");
            //exit(1);
            vector<int> t = {-1, -1, -1, 0, 1};

            Matrix<ban, 5, 1> x;


            flag = i_big_m_test(A, b, c, t, tol, x, optimal_value,th);

            //cout<<x<<endl;
            cout<<optimal_value<<endl;
            cout<<flag<<endl;
        }

        if(info == 5){
            // full na
            Matrix<ban, 9, 11> A;
            Vector<ban, 9> b;
            Vector<ban, 11> c;

            A <<  -1,            -1,            1, 0, 0, 0, 0, 0, 0, 0, 0,
                -15,           -30 - ETA,      0, 1, 0, 0, 0, 0, 0, 0, 0,
                1,            -1,             0, 0, 1, 0, 0, 0, 0, 0, 0,
                30 - ETA,     15,             0, 0, 0, 1, 0, 0, 0, 0, 0,
                1,             1,             0, 0, 0, 0, 1, 0, 0, 0, 0,
                5 - 2 * ETA,  10 - 2 * ETA,   0, 0, 0, 0, 0, 1, 0, 0, 0,
                1,             1,             0, 0, 0, 0, 0, 0, 1, 0, 0,
                35,           35 - 2 * ETA,   0, 0, 0, 0, 0, 0, 0, 1, 0,
                -1,            1,             0, 0, 0, 0, 0, 0, 0, 0, 1;

            b << -30,
                -900 + 45 * ETA - ETA * ETA,
                60,
                1800 - 75 * ETA - 2 * ETA * ETA,
                75,
                525 - 145 * ETA,
                70,
                2450 - 70 * ETA - 2 * ETA * ETA,
                70;

            c = Vector<ban, Dynamic>::Zero(11);
            c(0) = 8 + 14 * ETA;
            c(1) = 12 +  10 * ETA;

            vector<uint> B = {0, 1, 4, 5, 6, 7, 8, 9, 10};
            Vector<ban, 9> x;

            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value, th);
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 6){
            // full na i-big-m (c=0)
            Matrix<ban, 9, 2> A;
            Matrix<ban, 9, 1> b;
            Matrix<ban, 2, 1> c;

            A <<  1,             1,
                    15,            30 + ETA,
                    1,            -1,
                    30 + ETA,     15,
                    1,             1,
                    5 - 2 * ETA,  10 - 2 * ETA,
                    1,             1,
                    35,            35 - 2 * ETA,
                    -1,            1;

            b << 30,
                900 -45 * ETA + ETA * ETA,
                60,
                1800 - 75 * ETA - 2 * ETA * ETA,
                75,
                525 - 145 * ETA,
                70,
                2450 - 70 * ETA - 2 * ETA * ETA,
                70;

            //This c is for testing
            //c << 8 + 14 * ETA, 12 +  10 * ETA;
            vector<int> t = {1, 1, -1, -1, -1, -1, -1, -1, -1};
            Matrix<ban, 9,1> x;
            flag = i_big_m_test(A, b, c, t, tol, x, optimal_value, th);
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 7){
            // full na i-big-m
            Matrix<ban, 9, 2> A;
            Matrix<ban, 9, 1> b;
            Matrix<ban, 2, 1> c;

            A <<  1,             1,
                    15,            30 + ETA,
                    1,            -1,
                    30 + ETA,     15,
                    1,             1,
                    5 - 2 * ETA,  10 - 2 * ETA,
                    1,             1,
                    35,            35 - 2 * ETA,
                    -1,            1;

            b << 30,
                900 -45 * ETA + ETA * ETA,
                60,
                1800 - 75 * ETA - 2 * ETA * ETA,
                75,
                525 - 145 * ETA,
                70,
                2450 - 70 * ETA - 2 * ETA * ETA,
                70;

            //This c is for testing
            c << 8 + 14 * ETA, 12 +  10 * ETA;
            vector<int> t = {1, 1, -1, -1, -1, -1, -1, -1, -1};
            Matrix<ban, 9,1> x;
            flag = i_big_m_test(A, b, c, t, tol, x, optimal_value, th);
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 8){
            // full na with infinite coefficients --> the threshold is not enough because we consider only the value of monosemnia
            //and not the infinite coefficient!!
            Matrix<ban, 9, 11> A;
            Vector<ban, 9> b;
            Vector<ban, 11> c;

            A <<  -ALPHA,            -1,   1, 0, 0, 0, 0, 0, 0, 0, 0,
                -15,           -30 - ETA,      0, 1, 0, 0, 0, 0, 0, 0, 0,
                1,            -1,             0, 0, 1, 0, 0, 0, 0, 0, 0,
                30 - ETA,     15,             0, 0, 0, 1, 0, 0, 0, 0, 0,
                1,             1,             0, 0, 0, 0, 1, 0, 0, 0, 0,
                5 - 2 * ETA,  10 - 2 * ETA,   0, 0, 0, 0, 0, 1, 0, 0, 0,
                1,             ALPHA,             0, 0, 0, 0, 0, 0, 1, 0, 0,
                35,           35 - 2 * ETA,   0, 0, 0, 0, 0, 0, 0, 1, 0,
                -1,            1,             0, 0, 0, 0, 0, 0, 0, 0, 1;

            b << -30,
                -900 + 45 * ETA - ETA * ETA,
                60,
                1800 - 75 * ETA - 2 * ETA * ETA,
                75,
                525 - 145 * ETA,
                70,
                2450 - 70 * ETA - 2 * ETA * ETA,
                70;

            c = Vector<ban, Dynamic>::Zero(11);
            c(0) = 8 + 14 * ETA;
            c(1) = 12 +  10 * ETA;

            vector<uint> B = {0, 1, 4, 5, 6, 7, 8, 9, 10};
            Vector<ban, 9> x;

            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value, th);
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 9){
            // full na with infinitesimal coefficients
            Matrix<ban, 9, 11> A;
            Vector<ban, 9> b;
            Vector<ban, 11> c;

            A <<  -1,            -1,            1, 0, 0, 0, 0, 0, 0, 0, 0,
                -15,           - ETA,      0, 1, 0, 0, 0, 0, 0, 0, 0,
                1,            -1,             0, 0, 1, 0, 0, 0, 0, 0, 0,
                30 - ETA,     15,             0, 0, 0, 1, 0, 0, 0, 0, 0,
                1,             1,             0, 0, 0, 0, 1, 0, 0, 0, 0,
                -2 * ETA,  10 - 2 * ETA,   0, 0, 0, 0, 0, 1, 0, 0, 0,
                1,             1,             0, 0, 0, 0, 0, 0, 1, 0, 0,
                35,           35 - 2 * ETA,   0, 0, 0, 0, 0, 0, 0, 1, 0,
                -1,            1,             0, 0, 0, 0, 0, 0, 0, 0, 1;

            b << -30,
                -900 + 45 * ETA - ETA * ETA,
                60,
                1800 - 75 * ETA - 2 * ETA * ETA,
                75,
                525 - 145 * ETA,
                70,
                2450 - 70 * ETA - 2 * ETA * ETA,
                70;

            c = Vector<ban, Dynamic>::Zero(11);
            c(0) = 8 + 14 * ETA;
            c(1) = 12 +  10 * ETA;

            vector<uint> B = {0, 1, 4, 5, 6, 7, 8, 9, 10};
            Vector<ban, 9> x;

            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value, th);
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 10){
            // standard
            Matrix<ban, 6, 9> A;
            Matrix<ban, 6, 1> b;
            Matrix<ban, 9, 1> c;

            A <<  2*ALPHA,  1, -3, 1, 0, 0, 0, 0, 0,
                2,  3, -2, 0, 1, 0, 0, 0, 0,
                4,  3*ALPHA,  3, 0, 0, 1, 0, 0, 0,
                0,  0,  1, 0, 0, 0, 1, 0, 0,
                0,  0, -1, 0, 0, 0, 0, 1, 0,
                -1, -2, -1, 0, 0, 0, 0, 0, 1;

            b << 90, 190, 300, 10, -10, -70;

            c << 8, 12, 7, 0, 0, 0, 0, 0, 0;

            vector<uint> B = { 1, 2, 3, 4, 5, 7};
            Vector<ban, 9> x;

            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value,th);
            cout << "Standard Kite" << endl;
            //debug1(x, "x");
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 11){
            // full na i-big-m
            Matrix<ban, 9, 2> A;
            Matrix<ban, 9, 1> b;
            Matrix<ban, 2, 1> c;

            A <<  1,                1,
                    15,            30,
                    1,            -1,
                    30,             15,
                    1,             1,
                    5,              10,
                    1,             1,
                    35,            35,
                    -1,            1;

            b << 30,
                900,
                60,
                1800,
                75,
                525,
                70,
                2450,
                70;

            //This c is for testing
            //c << 8 + 14 * ETA, 12 +  10 * ETA;
            vector<int> t = {1, 1, -1, -1, -1, -1, -1, -1, -1};
            Matrix<ban, 9,1> x;
            flag = i_big_m(A, b, c, t, tol, x, optimal_value);
            cout << optimal_value << endl;
            cout << flag << endl;
            return 0;
        }

        if(info == 12){
            //case 5 standard (without bans)
            Matrix<ban, 9, 11> A;
            Vector<ban, 9> b;
            Vector<ban, 11> c;

            A <<  -1,            -1,            1, 0, 0, 0, 0, 0, 0, 0, 0,
                -15,           -30,            0, 1, 0, 0, 0, 0, 0, 0, 0,
                1,            -1,             0, 0, 1, 0, 0, 0, 0, 0, 0,
                30,           15,             0, 0, 0, 1, 0, 0, 0, 0, 0,
                1,             1,             0, 0, 0, 0, 1, 0, 0, 0, 0,
                5,            10,             0, 0, 0, 0, 0, 1, 0, 0, 0,
                1,             1,             0, 0, 0, 0, 0, 0, 1, 0, 0,
                35,           35,             0, 0, 0, 0, 0, 0, 0, 1, 0,
                -1,            1,             0, 0, 0, 0, 0, 0, 0, 0, 1;

            b << -30,
                -900,
                60,
                1800,
                75,
                525,
                70,
                2450,
                70;

            c = Vector<ban, Dynamic>::Zero(11);
            c(0) = 8;// + 14 * ETA;
            c(1) = 12;// +  10 * ETA;

            vector<uint> B = {0, 1, 4, 5, 6, 7, 8, 9, 10};
            Vector<ban, 9> x;
        
            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value, th);
            cout << optimal_value << endl;
            cout << flag << endl;
        }

        if(info == 13){
            // standard
            Matrix<ban, 6, 9> A;
            Matrix<ban, 6, 1> b;
            Matrix<ban, 9, 1> c;

            A <<  2*ETA,  1, -3, 1, 0, 0, 0, 0, 0,
                2,  3, -2, 0, 1, 0, 0, 0, 0,
                4,  3*ETA,  3, 0, 0, 1, 0, 0, 0,
                0,  0,  1, 0, 0, 0, 1, 0, 0,
                0,  0, -1, 0, 0, 0, 0, 1, 0,
                -1, -2, -1, 0, 0, 0, 0, 0, 1;

            b << 90, 190, 300, 10, -10, -70;

            c << 8, 12, 7, 0, 0, 0, 0, 0, 0;

            vector<uint> B = { 1, 2, 3, 4, 5, 7};
            Vector<ban, 9> x;

            flag = na_simplex_test(A, b, c, B, tol, x, optimal_value,th);
            cout << "Standard Kite" << endl;
            //debug1(x, "x");
            cout << optimal_value << endl;
            cout << flag << endl;
        }
    }
    /*
    Matrici da testare:
    -autovalori finiti/non finiti
    -singolari/non singolari
    -simmetrie/non simmetrie
    -coefficienti finiti/infiniti (alcuni) / infinitesimi (alcuni)
    Print:
    -fattorizzazione intermedia
    -fattorizzazione finale
    -eliminazione gaussiana (divisione E sottrazione)
    Testare con e senza threshold e documentare quando funziona e quando no e perché
    I test si fanno con i monosegni perché l'errore è introdotto da quelli e non dagli eta nei casi non-standard 
    */
    return 0;
}
