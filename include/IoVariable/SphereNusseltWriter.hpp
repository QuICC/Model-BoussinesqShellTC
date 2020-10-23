/** 
 * @file SphereNusseltWriter.hpp
 * @brief Implementation of the Nusselt number in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_SPHERENUSSELTWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERENUSSELTWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the Nusselt number in a sphere
    */
   class SphereNusseltWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereNusseltWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereNusseltWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Write State to file
          */
         virtual void write();

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const; 
         
      protected:
         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

      private:
         /**
          * @brief Nusselt number
          */
         MHDFloat mNusselt;

         /*
          * @brief Spherical volume to normalize energy to energy density
          */
         MHDFloat mTb;

         /**
          * @brief Origin projector
          */
         Matrix mOrigin;
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphereNusseltWriter> SharedSphereNusseltWriter;

   inline bool SphereNusseltWriter::isHeavy() const
   {
      return false;
   }

}
}

#endif // QUICC_IO_VARIABLE_SPHERENUSSELTWRITER_HPP
