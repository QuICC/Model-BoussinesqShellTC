/** 
 * @file NusseltWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef NUSSELTWRITER_HPP
#define NUSSELTWRITER_HPP

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

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII Nusselt number writer
    */
   class NusseltWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         NusseltWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~NusseltWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<NusseltWriter> SharedNusseltWriter;

}
}

#endif // NUSSELTWRITER_HPP
