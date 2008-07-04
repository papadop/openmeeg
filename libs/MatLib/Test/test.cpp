/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include "MatLibConfig.h"
#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"
#include "chrono.h"

#include <cmath>
#include <iostream>

using namespace std;

template<class T>
void genericTest(T &M)
{
    cout<<" Generic Test "<<endl;
    cout<<"   nlin  = " << (int)M.nlin() << endl;
    cout<<"   ncol  = " << (int)M.ncol() << endl;
    vecteur v(M.ncol());
    v.set(1);
    v = M*v;

    cout << endl << "BASE :" << endl;
    cout << M;

    // Test IO
    cout << endl << "BIN :" << endl;
    M.saveBin("tmp.bin");
    M.loadBin("tmp.bin");
    cout << M;

    cout << endl << "TXT :" << endl;
    M.saveTxt("tmp.txt");
    M.loadTxt("tmp.txt");
    cout << M;

    cout << "   operator * OK" << endl;
    cout.flush();
}

int main ()
{
    // section vecteur
    cout<<endl<<"========== vectors =========="<<endl;
    vecteur v(8);
    v.set(0);
    v.saveBin("tmp.bin");
    for(int i=0;i<8;i++) v(i)=i;
    v.saveTxt("tmp.txt");
    v.loadBin("tmp.bin");
    cout<<"v= "<<endl<<v<<endl;
    v.loadTxt("tmp.txt");
    cout<<"v= "<<endl<<v<<endl;

    // section matrice
    cout<<endl<<"========== matrices =========="<<endl;
    matrice M(5,5);

    for(size_t i=0;i<M.nlin();i++)
        for(size_t j=0;j<M.ncol();j++)
            M(i,j)=pow(2.0,(double)i)+pow(2.0,(double)j);

    genericTest(M);

    matrice Q = M.submat(3,1,2,3); // select submatrice
    cout<<"Q= "<<endl<<Q<<endl;

    matrice P(3,3);
    P(0,0) = 25 ; P(0,1) = 3 ; P(0,2) = 6 ;
    P(1,0) = 12 ; P(1,1) = 5 ; P(1,2) = 32 ;
    P(2,0) = 4 ; P(2,1) = 10 ; P(2,2) = 4 ;
    cout << "Matrice P : " << endl;
    cout << P;

    matrice Pinv = P.inverse();
    cout << "P Inverse matrice : " << endl;
    cout << Pinv;

    matrice unit = P*Pinv;
    double eps = 0.01;
    for(unsigned int i = 0; i<unit.nlin(); i++)
        for(unsigned int j = 0; j<unit.ncol(); j++){
            if(i == j){
                if(abs(unit(i,j)-1) > eps){
                    cerr << "Error: inverse is WRONG-1" << endl;
                    exit(1);
                }
            }
            else
                if(abs(unit(i,j)) > eps){
                    cerr << "Error: inverse is WRONG-2 " << "unit(" << i << "," << j << ") = " << unit(i,j) << endl;
                    exit(1);
                }
        }

#ifdef USE_MATIO
    cout << "MAT :" << endl;
    M.saveMat("tmp_matrice.mat");
    M.loadMat("tmp_matrice.mat");
    cout << M;
#endif

    // section symmatrice
    cout<<endl<<"========== symmetric matrices =========="<<endl;
    symmatrice S(4);
    for(unsigned int i=0; i<4; i++)
        for(unsigned int j=i; j<4; j++)
            S(i,j)=pow(2.0,(double)i)+pow(3.0,(double)j);

    genericTest(S);
    matrice R = S(1,2,0,2); // extract submatrice
    cout << "R= " << endl << R << endl;

    // section sparse_matrice
    cout<<endl<<"========== sparse matrices =========="<<endl;
    sparse_matrice spM(10,10);
    size_t _n=0;
    for(size_t i=0;i<5;i++)
    {
        _n=(_n*1237+1493)%1723;
        int _p=(_n*1237+1493)%1723;
        spM(_n%10,_p%10)=_n;
    }
    genericTest(spM);

#ifdef USE_MATIO
    cout << "MAT :" << endl;
    spM.saveMat("tmp_sparse_matrice.mat");
    spM.loadMat("tmp_sparse_matrice.mat");
    cout << spM;
#endif

    matrice U(10,10);
    U.set(1.0);
    matrice T = spM*U;
    cout << T;
    T = U*spM;
    cout << T;

    cout<<endl<<"========== fast sparse matrices =========="<<endl;
    cout << spM;
    fast_sparse_matrice fspM(spM);
    cout << fspM;
    genericTest(fspM);

    return 0;
}
