/** 
 * @file Cartesian1DNusseltDZWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer through the Z boundary extracted from temperature field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DNUSSELTDZWRITER_HPP
#define CARTESIAN1DNUSSELTDZWRITER_HPP

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
    * @brief Implementation of the ASCII Nusselt number writer through the Z boundary extracted from temperature field
    */
   class Cartesian1DNusseltDZWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DNusseltDZWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DNusseltDZWriter();

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
   typedef SharedPtrMacro<Cartesian1DNusseltDZWriter> SharedCartesian1DNusseltDZWriter;

}
}

#endif // CARTESIAN1DNUSSELTDZWRITER_HPP
