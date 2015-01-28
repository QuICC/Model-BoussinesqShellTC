/** 
 * @file IVariableHeavyAsciiEWriter.hpp
 * @brief Implementation of a generic variable to ASCII file writer with heavy calculations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IVARIABLEHEAVYASCIIEWRITER_HPP
#define IVARIABLEHEAVYASCIIEWRITER_HPP

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of a generic variable to ASCII file writer with heavy calculations
    */
   class IVariableHeavyAsciiEWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name       Filename
          * @param ext        File extension
          * @param header     Header string of file
          * @param type       Type string of file
          * @param version    Version string of file
          * @param id         ID of the dimension space
          */
         IVariableHeavyAsciiEWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id);

         /**
          * @brief Destructor
          */
         virtual ~IVariableHeavyAsciiEWriter();

         /**
          * @brief Perform heavy calculations
          */
         virtual void compute(Transform::TransformCoordinatorType& coord) = 0;

      protected:

      private:
   };

   /// Typedef for a smart reference counting pointer of a Variable HDF5 numbering writer
   typedef SharedPtrMacro<IVariableHeavyAsciiEWriter>   SharedIVariableHeavyAsciiEWriter;

}
}

#endif // IVARIABLEHEAVYASCIIEWRITER_HPP
