/** 
 * @file Cartesian1DNusseltZWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer through Z boundary
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DNUSSELTZWRITER_HPP
#define CARTESIAN1DNUSSELTZWRITER_HPP

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
   class Cartesian1DNusseltZWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DNusseltZWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DNusseltZWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian1DNusseltZWriter> SharedCartesian1DNusseltZWriter;

}
}

#endif // CARTESIAN1DNUSSELTZWRITER_HPP
