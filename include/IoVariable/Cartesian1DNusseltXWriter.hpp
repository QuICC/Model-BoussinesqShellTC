/** 
 * @file Cartesian1DNusseltXWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer through the X boundary
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DNUSSELTXWRITER_HPP
#define CARTESIAN1DNUSSELTXWRITER_HPP

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
    * @brief Implementation of the ASCII Nusselt number writer through the X boundary
    */
   class Cartesian1DNusseltXWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DNusseltXWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DNusseltXWriter();

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
   typedef SharedPtrMacro<Cartesian1DNusseltXWriter> SharedCartesian1DNusseltXWriter;

}
}

#endif // CARTESIAN1DNUSSELTXWRITER_HPP
