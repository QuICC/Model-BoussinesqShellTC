/** 
 * @file Cartesian3DNusseltZWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer for a 3D box through Z
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN3DNUSSELTZWRITER_HPP
#define CARTESIAN3DNUSSELTZWRITER_HPP

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
#include "IoVariable/IVariableAsciiWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII Nusselt number writer for a 3D box through Z
    */
   class Cartesian3DNusseltZWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         Cartesian3DNusseltZWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian3DNusseltZWriter();

         /**
          * @brief Initialise the operator and file
          */
         virtual void init();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Nusselt calculation operator
          */
         SparseMatrix   mNusseltOp;

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian3DNusseltZWriter> SharedCartesian3DNusseltZWriter;

}
}

#endif // CARTESIAN3DNUSSELTZWRITER_HPP
