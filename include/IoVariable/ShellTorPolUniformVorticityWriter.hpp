/** 
 * @file ShellTorPolUniformVorticityWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHELLTORPOLUNIFORMVORTICITYWRITER_HPP
#define SHELLTORPOLUNIFORMVORTICITYWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
    */
   class ShellTorPolUniformVorticityWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellTorPolUniformVorticityWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellTorPolUniformVorticityWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute energy for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Write State to file
          */
         virtual void write();

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const;
         
      protected:

      private:
         /**
          * @brief Spherical volume to normalize vorticity
          */
         MHDFloat mVolume;

         /*
          * @brief Thickness of 10* boundary layer
          */
         MHDFloat mDelta;

         /**
          * @brief Storage for the x component of UV
          */
         MHDFloat mUVx;

         /**
          * @brief Storage for the y component of UV
          */
         MHDFloat mUVy;

         /**
          * @brief Storage for the z component of UV
          */
         MHDFloat mUVz;

         /**
          * @brief Chebyshev operator to integrate in radius
          */
         SparseMatrix mIntgOp;

         /**
          * @brief Chebyshev operator for spherical integral in radius (include r^2 factor)
          */
         SparseMatrix mSphIntgOp;
   };

   inline bool ShellTorPolUniformVorticityWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<ShellTorPolUniformVorticityWriter> SharedShellTorPolUniformVorticityWriter;

}
}

#endif // SHELLTORPOLUNIFORMVORTICITYWRITER_HPP
