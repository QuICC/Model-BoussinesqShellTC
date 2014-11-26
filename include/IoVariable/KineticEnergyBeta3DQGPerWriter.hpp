/** 
 * @file KineticEnergyBeta3DQGPerWriter.hpp
 * @brief Implementation of the ASCII kinetic energy writer for the Beta 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef KINETICENERGYBETA3DQGPERWRITER_HPP
#define KINETICENERGYBETA3DQGPERWRITER_HPP

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
   class KineticEnergyBeta3DQGPerWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         KineticEnergyBeta3DQGPerWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~KineticEnergyBeta3DQGPerWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<KineticEnergyBeta3DQGPerWriter> SharedKineticEnergyBeta3DQGPerWriter;

}
}

#endif // KINETICENERGYBETA3DQGPERWRITER_HPP
