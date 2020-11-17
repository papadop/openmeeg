/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

/// \file
/// \brief File containing the integral operators.

#pragma once

#include <iostream>
#include <iomanip>

#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <symm_block_matrix.h>
#include <sparse_matrix.h>
#include <geometry.h>
#include <integrator.h>
#include <analytics.h>

#include <progressbar.h>

namespace OpenMEEG {

    // TODO: Use overloading and remove the internal suffix.

    void operatorFerguson(const Vect3&,const Mesh&,Matrix&,const unsigned&,const double&);
    void operatorDipolePotDer(const Vect3&,const Vect3&,const Mesh&,Vector&,const double&,const unsigned,const bool);
    void operatorDipolePot(const Vect3&,const Vect3&,const Mesh&,Vector&,const double&,const unsigned,const bool);

    namespace Details {
        // #define ADAPT_LHS

        template <template <typename,typename> class Integrator>
        void operatorDipolePot(const Vect3& r0,const Vect3& q,const Mesh& m,Vector& rhs,const double& coeff,const unsigned gauss_order) {
            static analyticDipPot anaDP;

            anaDP.init(q,r0);
            Integrator<double,analyticDipPot> gauss(0.001);
            gauss->setOrder(gauss_order);

            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const auto& triangle : m.triangles()) {
            #elif defined OPENMP_ITERATOR
            for (Triangles::const_iterator tit=m.triangles().begin();tit<m.triangles().end();++tit) {
                const Triangle& triangle = *tit;
            #else
            for (int i=0;i<m.triangles().size();++i) {
                const Triangle& triangle = *(m.triangles().begin()+i);
            #endif
                const double d = gauss->integrate(anaDP,triangle);
                #pragma omp critical
                rhs(triangle.index()) += d*coeff;
            }
        }

        inline Vect3 operatorFerguson(const Vect3& x,const Vertex& V,const Mesh& m) {
            Vect3 result;
            result = 0.0;

            //  Loop over triangles of which V is a vertex

            for (const auto& tp : m.triangles(V)) {
                const Triangle& T    = *tp;
                const Edge&     edge = T.edge(V);

                // A, B are the two opposite vertices to V (triangle A, B, V)
                const Vertex& A = edge.vertex(0);
                const Vertex& B = edge.vertex(1);
                const Vect3 AB = (A-B)/(2*T.area());
                
                analyticS analyS(V,A,B);
                const double opS = analyS.f(x);

                result += (AB*opS);
            }

            return result;
        }
    }

    class DiagonalBlock {
    public:

        DiagonalBlock(const Mesh& m): mesh(m) { }

        template <typename T>
        void addId(const double& coeff,T& matrix) const {
            // The Matrix is incremented by the identity P1P0 operator
            std::cout << "OPERATOR P1P0... (arg : mesh " << mesh.name() << " )" << std::endl;
            for (const auto& triangle : mesh.triangles())
                for (const auto& vertex : triangle)
                    matrix(triangle.index(),vertex->index()) += Id(triangle,*vertex)*coeff;
        }

    private:

        static double Id(const Triangle& T,const Vertex& V) {
            return (T.contains(V)) ? T.area()/3.0 : 0.0;
        }

        const Mesh& mesh;
    };

    class PartialBlock {
    public:

        PartialBlock(const Mesh& m): mesh(m) { }

        void addD(const double& coeff,const Vertices& points,Matrix& matrix) const {
            std::cout << "PARTAL OPERATOR D..." << std::endl;
            for (const auto& triangle : mesh.triangles()) {
                analyticD3 analyD(triangle);
                for (const auto& vertex : points) {
                    const Vect3& integrals = analyD.f(vertex);
                    for (unsigned i=0;i<3;++i)
                        matrix(vertex.index(),triangle.vertex(i).index()) += integrals(i)*coeff;
                }
            }
        }

        void S(const double& coeff,const Vertices& points,Matrix& matrix) const {
            std::cout << "PARTIAL OPERATOR S..." << std::endl;
            for (const auto& triangle : mesh.triangles()) {
                const analyticS analyS(triangle);
                for (const auto& vertex : points)
                    matrix(vertex.index(),triangle.index()) = coeff*analyS.f(vertex);
            }
        }


    private:

        const Mesh& mesh;
    };

    template <bool VERBOSE=true>
    class Operators {

        class Bloc: public Matrix {

            typedef Matrix base;

        public:

            Bloc(const unsigned r0,const unsigned c0,const unsigned n,const unsigned m): i0(r0),j0(c0),base(n,m) { }

            double& operator()(const unsigned i,const unsigned j)       { return base::operator()(i-i0,j-j0); }
            double  operator()(const unsigned i,const unsigned j) const { return base::operator()(i-i0,j-j0); }

        private:

            const unsigned i0;
            const unsigned j0;
        };
        
        class SymBloc: public SymMatrix {

            typedef SymMatrix base;

        public:

            SymBloc(const unsigned off,const unsigned sz): offset(off),base(sz) { }

            double& operator()(const unsigned i,const unsigned j)       { return base::operator()(i-offset,j-offset); }
            double  operator()(const unsigned i,const unsigned j) const { return base::operator()(i-offset,j-offset); }

        private:

            const unsigned offset;
        };
        
        void message(const char* op_name) const {
            if (VERBOSE)
                std::cout << "OPERATOR " << std::left << std::setw(2) << op_name << "... (arg : mesh " << mesh1.name() << " , mesh " << mesh2.name() << " )" << std::endl;
        }

    public:

        // This constructor takes the following arguments:
        //  - The 2 interacting meshes.
        //  - The gauss order parameter (for adaptive integration).
        //  - A verbosity parameters (for printing the action on the terminal).

        Operators(const Mesh& m1,const Mesh& m2,const unsigned order): mesh1(m1),mesh2(m2),gauss_order(order),same_mesh(&mesh1==&mesh2) { }

        // SymMatrix is initialized at once.

        static void init(SymMatrix& matrix) { matrix.set(0.0); }
        void set_blocks(SymMatrix&) const   { }

        // SymmetricBlockMatrix is initialized blockwise.

        static void init(maths::SymmetricBlockMatrix&) { }
        void set_blocks(maths::SymmetricBlockMatrix& matrix) const {
            const Ranges& vrange1 = mesh1.vertices_ranges();
            const Ranges& vrange2 = mesh2.vertices_ranges();
            const Ranges& trange1 = { mesh1.triangles_range() };
            const Ranges& trange2 = { mesh2.triangles_range() };

            matrix.add_blocks(trange1,trange2); // S blocks.
            matrix.add_blocks(vrange1,vrange2); // N blocks.
            matrix.add_blocks(trange1,vrange2); // D blocks.
            if (!same_mesh)
                matrix.add_blocks(trange2,vrange1); // D* blocks when they are not the transpose of D blocks.
        }
        
        // Various operators take the following arguments:
        //  - The coefficient to be applied to each matrix element (depending on conductivities, ...)
        //  - The storage Matrix for the result

        template <typename T>
        void S(const double& coeff,T& matrix) const {
            message("S");
            ProgressBar pb(mesh1.triangles().size());

            // Operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(l,t) a P0 test function on layer l and triangle t.

            if (same_mesh) {

                // When meshes are equal, optimized computation for a symmetric matrix.

                const Triangles& triangles = mesh1.triangles();
                for (Triangles::const_iterator tit1=triangles.begin(); tit1!=triangles.end(); ++tit1,++pb) {
                    const Triangle& triangle1 = *tit1;
                    const analyticS analyS(triangle1);
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_ITERATOR
                    for (Triangles::const_iterator tit2=tit1; tit2!=triangles.end(); ++tit2) {
                        const Triangle& triangle2 = *tit2;
                    #else
                    for (int i2=tit1-triangles.begin(); i2<triangles.size(); ++i2) {
                        const Triangle& triangle2 = *(triangles.begin()+i2);
                    #endif
                        matrix(triangle1.index(),triangle2.index()) = S(analyS,triangle2)*coeff;
                    }
                }
            } else {
                // TODO check the symmetry of S. 
                // if we invert tit1 with tit2: results in HeadMat differs at 4.e-5 which is too big.
                // using ADAPT_LHS with tolerance at 0.000005 (for S) drops this at 6.e-6 (but increase the computation time).

                for (const auto& triangle1 : mesh1.triangles()) {
                    const analyticS analyS(triangle1);
                    const Triangles& m2_triangles = mesh2.triangles();
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_RANGEFOR
                    for (const auto& triangle2 : m2_triangles) {
                    #elif defined OPENMP_ITERATOR
                    for (Triangles::const_iterator tit2=m2_triangles.begin();tit2<m2_triangles.end();++tit2) {
                        const Triangle& triangle2 = *tit2;
                    #else
                    for (int i2=0;i2<m2_triangles.size();++i2) {
                        const Triangle& triangle2 = *(m2_triangles.begin()+i2);
                    #endif
                        matrix(triangle1.index(),triangle2.index()) = S(analyS,triangle2)*coeff;
                    }
                    ++pb;
                }
            }
        }

        template <typename T1,typename T2>
        void N(const double& coeff,const T1& S,T2& matrix) const {

            message("N");

            ProgressBar pb(mesh1.vertices().size());
            if (same_mesh) {

                // When meshes are equal, optimized computation for a symmetric matrix.

                for (auto vit1=mesh1.vertices().begin(); vit1!=mesh1.vertices().end(); ++vit1) {
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_ITERATOR
                    for (auto vit2=vit1; vit2<mesh1.vertices().end(); ++vit2) {
                    #else
                    for (int i2=0;i2<=vit1-mesh1.vertices().begin();++i2) {
                        const auto vit2 = mesh1.vertices().begin()+i2;
                    #endif
                        matrix((*vit1)->index(),(*vit2)->index()) += N(**vit1,**vit2,mesh1,S)*coeff;
                    }
                    ++pb;
                }
            } else {
                const VerticesRefs& m2_vertices = mesh2.vertices();
                for (const auto& vertex1 : mesh1.vertices()) {
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_RANGEFOR
                    for (const auto& vertex2 : m2_vertices) {
                    #elif defined OPENMP_ITERATOR
                    for (auto vit2=m2_vertices.begin(); vit2<m2_vertices.end(); ++vit2) {
                        const Vertex* vertex2 = *vit2;
                    #else
                    for (int i2=0; i2<m2_vertices.size(); ++i2) {
                        const Vertex* vertex2 = *(m2_vertices.begin()+i2);
                    #endif
                        matrix(vertex1->index(),vertex2->index()) += N(*vertex1,*vertex2,mesh1,mesh2,S)*coeff;
                    }
                    ++pb;
                }
            }
        }

        //  Operator N in case one of the mesh is a current barrier.

        template <typename T>
        void N(const double& coeff,T& matrix) const {
            if (same_mesh) { // When meshes are equal, optimized computation for a symmetric matrix.
                SymBloc Sbloc(mesh1.triangles().front().index(),mesh1.triangles().size());
                S(1.0,Sbloc);
                N(coeff,Sbloc,matrix);
            } else {
                Bloc Sbloc(mesh1.triangles().front().index(),mesh2.triangles().front().index(),mesh1.triangles().size(),mesh2.triangles().size());
                S(1.0,Sbloc);
                N(coeff,Sbloc,matrix);
            }
        }

        template <typename T>
        void D(const double& coeff,T& matrix) const {
            // This function (OPTIMIZED VERSION) has the following arguments:
            //    the coefficient to be applied to each matrix element (depending on conductivities, ...)
            //    the storage Matrix for the result

            message("D");
            D(mesh1,mesh2,coeff,matrix);
        }

        template <typename T>
        void Dstar(const double& coeff,T& matrix) const {
            // This function (OPTIMIZED VERSION) has the following arguments:
            //    the coefficient to be applied to each matrix element (depending on conductivities, ...)
            //    the storage Matrix for the result

            message("D*");
            D(mesh2,mesh1,coeff,matrix);
        }

    private:

        double S(const analyticS& analyS,const Triangle& triangle) const {
        #ifdef ADAPT_LHS
            AdaptiveIntegrator<double,analyticS> gauss(0.005);
        #else
            STATIC_OMP Integrator<double,analyticS> gauss;
        #endif
            gauss.setOrder(gauss_order);
            return gauss.integrate(analyS,triangle);
        }

        template <typename T>
        double N(const Vertex& V1,const Vertex& V2,const Mesh& m,const T& matrix) const {
            return N(0.25,V1,V2,m,m,matrix);
        }

        // This function assumes that the two meshes are different.

        template <typename T>
        double N(const Vertex& V1,const Vertex& V2,const Mesh& m1,const Mesh& m2,const T& matrix) const {
            const double coeff = (V1==V2) ? 0.5 : 0.25;
            return N(coeff,V1,V2,m1,m2,matrix);
        }

        template <typename T>
        double N(const double& factor,const Vertex& V1,const Vertex& V2,const Mesh& m1,const Mesh& m2,const T& matrix) const {

            double result = 0.0;
            for (const auto& tp1 : m1.triangles(V1)) {
                const Edge& edge1 = tp1->edge(V1);
                const Vect3& CB1 = edge1.vertex(0)-edge1.vertex(1);
                for (const auto& tp2 : m2.triangles(V2)) {

                    const Edge& edge2 = tp2->edge(V2);
                    const Vect3& CB2 = edge2.vertex(0)-edge2.vertex(1);

                    result -= factor*dotprod(CB1,CB2)*matrix(tp1->index(),tp2->index())/(tp1->area()*tp2->area());
                }
            }
            return result;
        }

        template <typename T>
        void D(const Triangle& T1,const Triangle& T2,T& mat,const double& coeff) const {
            //this version of operatorD add in the Matrix the contribution of T2 on T1
            // for all the P1 functions it gets involved
            // consider varying order of quadrature with the distance between T1 and T2
            analyticD3 analyD(T2);

        #ifdef ADAPT_LHS
            AdaptiveIntegrator<Vect3,analyticD3> gauss(0.005);
            gauss.setOrder(gauss_order);
        #else
            STATIC_OMP Integrator<Vect3, analyticD3> gauss(gauss_order);
        #endif
            const Vect3 total = gauss.integrate(analyD,T1);

            for (unsigned i=0; i<3; ++i)
                mat(T1.index(),T2.vertex(i).index()) += total(i)*coeff;
        }

        template <typename T>
        void D(const Mesh& m1,const Mesh& m2,const double& coeff,T& mat) const {
            // This function (OPTIMIZED VERSION) has the following arguments:
            //    the 2 interacting meshes
            //    the storage Matrix for the result
            //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
            //    the gauss order parameter (for adaptive integration)

            // In this version of the function, in order to skip multiple computations of the same quantities
            //    loops are run over triangles but the Matrix cannot be filled in this function anymore
            //    That's why the filling is done is function Details::operatorD
            //

            ProgressBar pb(m1.triangles().size());
            const Triangles& m1_triangles = m1.triangles();
            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const auto& triangle1 : m1_triangles) {
            #elif defined OPENMP_ITERATOR
            for (Triangles::const_iterator tit1=m1_triangles.begin(); tit1<m1_triangles.end(); ++tit1) {
                const Triangle& triangle1 = *tit1;
            #else
            for (int i1=0; i1<m1_triangles.size(); ++i1) {
                const Triangle& triangle1 = *(m1_triangles.begin()+i1);
            #endif
                for (const auto& triangle2 : m2.triangles())
                    D(triangle1,triangle2,mat,coeff);
                ++pb;
            }
        }

        const Mesh&    mesh1;
        const Mesh&    mesh2;
        const unsigned gauss_order;
        const bool     same_mesh;
    };
}
