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

#include <operators.h>

namespace OpenMEEG {

    // General routine for applying Details::operatorFerguson (see this function for further comments)
    // to an entire mesh, and storing coordinates of the output in a Matrix.

    void operatorFerguson(const Vect3& x,const Mesh& m,Matrix& mat,const unsigned& offsetI,const double& coeff) {
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& vertexp : m.vertices()) {
        #elif defined OPENMP_ITERATOR
        for (auto vit=m.vertices().begin();vit<m.vertices().end();++vit) {
            const Vertex* vertexp = *vit;
        #else
        for (int i=0;i<m.vertices().size();++i) {
            const Vertex* vertexp = *(m.vertices().begin()+i);
        #endif
            const unsigned vindex = vertexp->index();
            Vect3 v = Details::operatorFerguson(x,*vertexp,m);
            mat(offsetI+0,vindex) += v.x()*coeff;
            mat(offsetI+1,vindex) += v.y()*coeff;
            mat(offsetI+2,vindex) += v.z()*coeff;
        }
    }

    void operatorDipolePotDer(const Dipole& dipole,const Mesh& m,Vector& rhs,const double& coeff,const unsigned gauss_order,const bool adapt_rhs) {
        Integrator<Vect3,analyticDipPotDer>* gauss = (adapt_rhs) ? new AdaptiveIntegrator<Vect3,analyticDipPotDer>(gauss_order,0.001) :
                                                                   new Integrator<Vect3,analyticDipPotDer>(gauss_order);

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
            const analyticDipPotDer anaDPD(dipole,triangle);

            const Vect3& v = gauss->integrate(anaDPD,triangle);
            #pragma omp critical
            {
                for (unsigned i=0; i<3; ++i)
                    rhs(triangle.vertex(i).index()) += v(i)*coeff;
            }
        }
        delete gauss;
    }

    void operatorDipolePot(const Dipole& dipole,const Mesh& m,Vector& rhs,const double& coeff,const unsigned gauss_order,const bool adapt_rhs) {
        const analyticDipPot anaDP(dipole);

        AdaptiveIntegrator<double,analyticDipPot> gauss(gauss_order,0.001);

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
            const double d = gauss.integrate(anaDP,triangle);
            rhs(triangle.index()) += d*coeff;
        }
    }
}
