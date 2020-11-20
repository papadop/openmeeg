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

#pragma once

#include <cmath>
#include <iostream>

#include <vertex.h>
#include <triangle.h>
#include <mesh.h>

namespace OpenMEEG {

    // light class containing d Vect3

    template <unsigned d>
    class OPENMEEG_EXPORT Vect3array {

        Vect3 t[d];

    public:

        Vect3array() {};

        inline Vect3array(const double x) {
            for (unsigned i=0;i<d;++i)
                t[i] = Vect3(x);
        }

        inline Vect3array<d> operator*(const double x) const {
            Vect3array<d> r;
            for (unsigned i=0;i<d;++i)
                r.t[i] = t[i]*x;
            return r;
        }

        inline Vect3  operator()(const int i) const { return t[i]; }
        inline Vect3& operator()(const int i)       { return t[i]; }
    };

    // Quadrature rules are from Marc Bonnet's book: Equations integrales..., Appendix B.3

    static const double cordBars[4][16][4] = {
        //parameters for N=3
        {
            {0.166666666666667, 0.166666666666667, 0.666666666666667, 0.166666666666667},
            {0.166666666666667, 0.666666666666667, 0.166666666666667, 0.166666666666667},
            {0.666666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        }
        ,
        // parameters for N=6
        {
            {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.111690794839005},
            {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.111690794839005},
            {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.111690794839005},
            {0.091576213509771, 0.091576213509771, 0.816847572980458, 0.054975871827661},
            {0.091576213509771, 0.816847572980458, 0.091576213509771, 0.054975871827661},
            {0.816847572980458, 0.091576213509771, 0.091576213509771, 0.054975871827661},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        }
        ,
            // parameters for N=7
        {
            {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.1125},
            {0.470142064105115, 0.470142064105115, 0.059715871789770, 0.066197076394253},
            {0.470142064105115, 0.059715871789770, 0.470142064105115, 0.066197076394253},
            {0.059715871789770, 0.470142064105115, 0.470142064105115, 0.066197076394253},
            {0.101286507323456, 0.101286507323456, 0.797426985353088, 0.062969590272414},
            {0.101286507323456, 0.797426985353088, 0.101286507323456, 0.062969590272414},
            {0.797426985353088, 0.101286507323456, 0.101286507323456, 0.062969590272414},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        }
        ,

            // parameters for N=16
        {
            {0.333333333333333, 0.333333333333333, 0.333333333333333, 0.072157803838893},
            {0.081414823414554, 0.459292588292722, 0.459292588292722, 0.047545817133642},
            {0.459292588292722, 0.081414823414554, 0.459292588292722, 0.047545817133642},
            {0.459292588292722, 0.459292588292722, 0.081414823414554, 0.047545817133642},
            {0.898905543365937, 0.050547228317031, 0.050547228317031, 0.016229248811599},
            {0.050547228317031, 0.898905543365937, 0.050547228317031, 0.016229248811599},
            {0.050547228317031, 0.050547228317031, 0.898905543365937, 0.016229248811599},
            {0.658861384496479, 0.170569307751760, 0.170569307751761, 0.051608685267359},
            {0.170569307751760, 0.658861384496479, 0.170569307751761, 0.051608685267359},
            {0.170569307751760, 0.170569307751761, 0.658861384496479, 0.051608685267359},
            {0.008394777409957, 0.728492392955404, 0.263112829634639, 0.013615157087217},
            {0.728492392955404, 0.008394777409957, 0.263112829634639, 0.013615157087217},
            {0.728492392955404, 0.263112829634639, 0.008394777409957, 0.013615157087217},
            {0.008394777409957, 0.263112829634639, 0.728492392955404, 0.013615157087217},
            {0.263112829634639, 0.008394777409957, 0.728492392955404, 0.013615157087217},
            {0.263112829634639, 0.728492392955404, 0.008394777409957, 0.013615157087217}
        }

    }; // end of gaussTriangleParams

    static const unsigned nbPts[4] = {3, 6, 7, 16};

    class OPENMEEG_EXPORT Integrator {

        static unsigned safe_order(const unsigned n) {
            if (n>0 && n<4)
                return n;
            std::cout << "Unavailable Gauss order " << n << ": min is 1, max is 3" << std::endl;
            return (n<1) ? 1 : 3;
        }

    protected:

        typedef Vect3 Point;
        typedef Point TrianglePoints[3];

    public:

        Integrator(const unsigned ord=3,const double tol=0.0): order(safe_order(ord)) {  }

        template <typename T,typename Function>
        T integrate(const Function& function,const Triangle& triangle) const {
            const TrianglePoints tripts = { triangle.vertex(0), triangle.vertex(1), triangle.vertex(2) };
            return triangle_integration(function,tripts);
        }

    protected:

        template <typename T,typename Function>
        T triangle_integration(const Function& function,const TrianglePoints& triangle) const {
            T result = 0.0;
            for (unsigned i=0;i<nbPts[order];++i) {
                Vect3 v(0.0,0.0,0.0);
                for (unsigned j=0;j<3;++j)
                    v.multadd(cordBars[order][i][j],triangle[j]);
                result += cordBars[order][i][3]*function(v);
            }

            // compute double area of triangle defined by points

            const double area2 = crossprod(triangle[1]-triangle[0],triangle[2]-triangle[0]).norm();
            return result*area2;
        }

        const unsigned order;
    };

    class OPENMEEG_EXPORT AdaptiveIntegrator: public Integrator {

        typedef Integrator base;

    public:

        AdaptiveIntegrator(const unsigned ord): base(ord),tolerance(0.0),max_depth(0) { }
        AdaptiveIntegrator(const unsigned ord,const double tol): base(ord),tolerance(tol),max_depth(10) { }
        AdaptiveIntegrator(const unsigned ord,const unsigned levels,const double tol=0.0001): base(ord),tolerance(tol),max_depth(levels) { }

        double norm(const double a) const { return fabs(a);  }
        double norm(const Vect3& a) const { return a.norm(); }

        // TODO: T can be deduced from Function.

        template <typename T,typename Function>
        T integrate(const Function& function,const Triangle& triangle) const {
            const TrianglePoints tripts = { triangle.vertex(0), triangle.vertex(1), triangle.vertex(2) };
            const T& coarse = base::triangle_integration<T>(function,tripts);
            return (max_depth==0) ? coarse : adaptive_integration<T>(function,tripts,coarse,max_depth);
        }

    private:

        using typename base::Point;
        using typename base::TrianglePoints;

        template <typename T,typename Function>
        T adaptive_integration(const Function& function,const TrianglePoints& triangle,const T& coarse,const unsigned level) const {
            const Point midpoints[] = { 0.5*(triangle[1]+triangle[2]), 0.5*(triangle[2]+triangle[0]), 0.5*(triangle[0]+triangle[1]) };
            const TrianglePoints new_triangles[] = {
                { triangle[0],  midpoints[1], midpoints[2] }, { midpoints[0], triangle[1],  midpoints[2] },
                { midpoints[0], midpoints[1], triangle[2]  }, { midpoints[0], midpoints[1], midpoints[2] }
            };

            T refined = 0.0;
            T integrals[4];
            for (unsigned i=0; i<4; ++i) {
                integrals[i] = base::triangle_integration<T>(function,new_triangles[i]);
                refined += integrals[i];
            }

            if (norm(coarse-refined)<=tolerance*norm(coarse) || level==0)
                return coarse;

            T sum = 0.0;
            for (unsigned i=0; i<4; ++i)
                sum += adaptive_integration<T>(function,new_triangles[i],integrals[i],level-1);
            return sum;
        }

        const double   tolerance;
        const unsigned max_depth;
    };
}
