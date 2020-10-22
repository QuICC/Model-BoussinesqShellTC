/** 
 * @file ContinuityWriter.hpp
 * @brief Implementation of the ASCII maximal continuity writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CONTINUITYWRITER_HPP
#define CONTINUITYWRITER_HPP

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
    * @brief Implementation of the ASCII maximal continuity writer
    */
   class ContinuityWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         ContinuityWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~ContinuityWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<ContinuityWriter> SharedContinuityWriter;

}
}

#endif // CONTINUITYWRITER_HPP
