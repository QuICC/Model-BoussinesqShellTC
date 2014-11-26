/** 
 * @file NusseltBeta3DQGPerWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer for the Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef NUSSELTBETA3DQGPERWRITER_HPP
#define NUSSELTBETA3DQGPERWRITER_HPP

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
    * @brief Implementation of the ASCII Nusselt number writer for the Beta 3DQG model
    */
   class NusseltBeta3DQGPerWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         NusseltBeta3DQGPerWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~NusseltBeta3DQGPerWriter();

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
   typedef SharedPtrMacro<NusseltBeta3DQGPerWriter> SharedNusseltBeta3DQGPerWriter;

}
}

#endif // NUSSELTBETA3DQGPERWRITER_HPP
