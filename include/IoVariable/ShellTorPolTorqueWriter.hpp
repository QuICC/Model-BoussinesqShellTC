/** 
 * @file ShellTorPolTorqueWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHELLTORPOLTORQUEWRITER_HPP
#define SHELLTORPOLTORQUEWRITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics inner sphere torque calculation
    * for a Toroidal/Poloidal field in a spherical shell
    */
   class ShellTorPolTorqueWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellTorPolTorqueWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellTorPolTorqueWriter();

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
          * @brief Torque value to be computed and printed
          */
         MHDFloat mTorque;

         /*
          * @brief bool flag to signal different cores whether they have the T_1^0 field or not
          */
         bool mComputeFlag;

         /*
          * @brief array storing the evaluated basis of r d/dr (1/r T) de factor beying a projection from spectral
          * to physical value at the boundary
          */
         Array mProj;

         /*
          * @brief precomputed E-dependent scaling for the torque
          */
         MHDFloat mFactor;
   };

   inline bool ShellTorPolTorqueWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<ShellTorPolTorqueWriter> SharedShellTorPolTorqueWriter;

}
}

#endif // SHELLTORPOLTORQUEWRITER_HPP
