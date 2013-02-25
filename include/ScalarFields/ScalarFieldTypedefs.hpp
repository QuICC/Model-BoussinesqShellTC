/** \file ScalarFieldTypedefs.hpp
 *  \brief Definition of some typedefs for the different scalar fields
 */

#ifndef SCALARFIELDTYPEDEFS_HPP
#define SCALARFIELDTYPEDEFS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "ScalarFields/FlatLayout.hpp"
#include "ScalarFields/SlicedLayout.hpp"
#include "ScalarFields/ScalarField1D.hpp"
#include "ScalarFields/ScalarField2D.hpp"
#include "ScalarFields/ScalarField3D.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @name Scalar field typedefs in 1D
    */
   //@{
   /// Typedef for a real flat 1D scalar field
   typedef ScalarField1D<EPMFloat, FlatLayout>  FlatDScalar1D;

   /// Typedef for a real sliced 1D scalar field
   typedef ScalarField1D<EPMFloat, SlicedLayout>  SlicedDScalar1D;

   /// Typedef for a complex flat 1D scalar field
   typedef ScalarField1D<EPMComplex, FlatLayout>  FlatZScalar1D;

   /// Typedef for a complex flat 1D scalar field
   typedef ScalarField1D<EPMComplex, SlicedLayout>  SlicedZScalar1D;
   //@}

   /**
    * @name Scalar field typedefs in 2D
    */
   //@{
   /// Typedef for a real flat 2D scalar field
   typedef ScalarField2D<EPMFloat, FlatLayout>  FlatDScalar2D;

   /// Typedef for a real sliced 2D scalar field
   typedef ScalarField2D<EPMFloat, SlicedLayout>  SlicedDScalar2D;

   /// Typedef for a complex flat 2D scalar field
   typedef ScalarField2D<EPMComplex, FlatLayout>  FlatZScalar2D;

   /// Typedef for a complex flat 2D scalar field
   typedef ScalarField2D<EPMComplex, SlicedLayout>  SlicedZScalar2D;
   //@}

   /**
    * @name Scalar field typedefs in 3D
    */
   //@{
   /// Typedef for a real flat 3D scalar field
   typedef ScalarField3D<EPMFloat, FlatLayout>  FlatDScalar3D;

   /// Typedef for a real sliced 3D scalar field
   typedef ScalarField3D<EPMFloat, SlicedLayout>  SlicedDScalar3D;

   /// Typedef for a complex flat 3D scalar field
   typedef ScalarField3D<EPMComplex, FlatLayout>  FlatZScalar3D;

   /// Typedef for a complex flat 3D scalar field
   typedef ScalarField3D<EPMComplex, SlicedLayout>  SlicedZScalar3D;
   //@}
}
}

#endif // SCALARFIELDTYPEDEFS_HPP
