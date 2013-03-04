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
#include "Enums/PhysicalNames.hpp"
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
         * @param spRes  Shared resolution information
         * @param type  Type of the file (typically scheme name)
         */
         StateFileWriter(std::string type);

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
          * @brief Write Physical parameters to file
          */
         void writePhysical();

         /**
          * @brief Write run information to file
          *
          * @param time Reached simulation time
          * @param step Current timestep
          */
         void writeRun(const MHDFloat time, const MHDFloat step);

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
