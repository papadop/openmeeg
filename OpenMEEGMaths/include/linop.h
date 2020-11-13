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

#include <cstdlib>
#include <cmath>

#ifdef HAVE_SHARED_PTR_ARRAY_SUPPORT
#include <memory>
template <typename T>
using SharedPtr = std::shared_ptr<T>;
#else
#include <boost/shared_ptr.hpp>
template <typename T>
using SharedPtr = boost::shared_ptr<T>;
#endif

#include "OpenMEEGMathsConfig.h"
#include <OMassert.H>

namespace OpenMEEG {

    namespace maths {
        struct OPENMEEGMATHS_EXPORT MathsIO;
    }

    // Properly convert a size_t int to an int

    OPENMEEGMATHS_EXPORT inline BLAS_INT sizet_to_int(const size_t& num) {
        BLAS_INT num_out = static_cast<BLAS_INT>(num);
        om_assert(num_out>=0);
        return num_out;
    }

    class OPENMEEGMATHS_EXPORT LinOpInfo {
    public:

        typedef maths::MathsIO* IO;

        typedef enum { FULL, SYMMETRIC, BLOCK, BLOCK_SYMMETRIC, SPARSE } StorageType;
        typedef unsigned Dimension;

        LinOpInfo() { }
        LinOpInfo(const size_t m,const size_t n,const StorageType st,const Dimension d):
            num_lines(m),num_cols(n),storage(st),dim(d)  { }

        virtual ~LinOpInfo() {};

        size_t  nlin() const { return num_lines; }
        size_t& nlin()       { return num_lines; }

        virtual size_t  ncol() const { return num_cols; }
                size_t& ncol()       { return num_cols; }

        StorageType  storageType() const { return storage; }
        StorageType& storageType()       { return storage; }

        Dimension   dimension()   const { return dim;     }
        Dimension&  dimension()         { return dim;     }

        IO& default_io() { return DefaultIO; }

    protected:

        size_t      num_lines;
        size_t      num_cols;
        StorageType storage;
        Dimension   dim;
        IO          DefaultIO;
    };

    class OPENMEEGMATHS_EXPORT LinOp: public LinOpInfo {

        typedef LinOpInfo base;

    public:

        LinOp() { }
        LinOp(const size_t m,const size_t n,const StorageType st,const Dimension d): base(m,n,st,d) { }

        virtual size_t size() const = 0;
        virtual void   info() const = 0;
    };

    typedef enum { DEEP_COPY } DeepCopy;

    struct OPENMEEGMATHS_EXPORT LinOpValue: public SharedPtr<double[]> {
        typedef SharedPtr<double[]> base;

        LinOpValue(): base(0) { }
        LinOpValue(const size_t n): base(new double[n]) { }
        LinOpValue(const size_t n,const double* initval): LinOpValue(n) { std::copy(initval,initval+n,&(*this)[0]); }
        LinOpValue(const size_t n,const LinOpValue& v):   LinOpValue(n,&(v[0])) { }

        ~LinOpValue() { }

        bool empty() const { return static_cast<bool>(*this); }
    };
}
