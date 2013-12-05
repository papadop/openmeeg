/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
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

#ifndef OPENMEEG_CGAL_MESH_INCLUDE_H
#define OPENMEEG_CGAL_MESH_INCLUDE_H

// for verbosity
#define CGAL_MESH_3_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Image_3.h>

namespace OpenMEEG {

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;

    template <typename C3t3>
    Mesh CGAL_to_OM(C3t3 c3t3) {
        typedef typename C3t3::Triangulation Tr;
        Mesh m(c3t3.triangulation().number_of_vertices(), c3t3.number_of_facets());

        std::map<typename C3t3::Vertex_handle, unsigned> V;
        unsigned inum = 0;
        for ( typename Tr::Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin(), end = c3t3.triangulation().finite_vertices_end(); vit != end; ++vit)
        {
            const typename Tr::Point& p = vit->point();
            m.add_vertex(Vertex(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())));
            V[vit] = inum++;
        }

        for ( typename C3t3::Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
            const typename Tr::Cell_handle cell = fit->first;
            const int& index = fit->second;
            const int index1 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 0))];
            const int index2 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 1))];
            const int index3 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 2))];
            Triangle t(m.vertices()[index1], m.vertices()[index2], m.vertices()[index3] );
            m.push_back(t);
        }

        m.update();
        return m;
    }
}
#endif  //! OPENMEEG_CGAL_MESH_INCLUDE_H
