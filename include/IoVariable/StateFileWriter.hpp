/** \file StateFileWriter.hpp
 *  \brief Implementation of the HDF5 state file writer
 */

#ifndef STATEFILEWRITER_HPP
#define STATEFILEWRITER_HPP

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
#include "IoVariable/IVariableHdf5NWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * \brief Implementation of the HDF5 state file writer
    */
   class StateFileWriter: public IVariableHdf5NWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type       Type of the file (typically scheme name)
          * @param isRegular  Is data regular?
          */
         StateFileWriter(std::string type, bool isRegular);

         /**
          * @brief Destructor
          */
         virtual ~StateFileWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:
         /**
          * @brief Create group for scalar field
          *
          * @param name    Name of the field
          * @param scalar  Scalar field values
          */
         void writeSpectralScalar(const std::string& name, const Datatypes::SpectralScalarType& scalar);

         /**
          * @brief Create group for vector field
          *
          * @param name    Name of the field
          * @param vector  Vector of components
          */
         void writeSpectralVector(const std::string& name, const std::vector<Datatypes::SpectralScalarType>& vector);

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<StateFileWriter> SharedStateFileWriter;

}
}

#endif // STATEFILEWRITER_HPP
