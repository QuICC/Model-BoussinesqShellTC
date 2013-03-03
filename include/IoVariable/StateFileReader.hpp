/** \file StateFileReader.hpp
 *  \brief Implementation of HDF5 state file reader
 */

#ifndef STATEFILEREADER_HPP
#define STATEFILEREADER_HPP

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
#include "IoVariable/IVariableHdf5Reader.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * \brief Implementation of HDF5 state file reader
    */
   class StateFileReader: public IVariableHdf5Reader
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name  Name of the file
         * @param type  Type of the file (typically scheme name)
         */
         StateFileReader(std::string name, std::string type);

         /**
         * @brief Destructor
         */
         virtual ~StateFileReader();

         /**
          * @brief Read State from file
          */
         virtual void read();

         /**
          * @brief Read resolution and data information from state file
          */
         void readSetup();

         /**
          * @brief Read Partial data from state file
          */
         template <typename TFilter> void readPartial();

         /**
          * @brief Get state file time
          */
         MHDFloat time() const;

         /**
          * @brief Get state file timestep
          */
         MHDFloat timestep() const;
         
      protected:
         /**
          * @brief Time read from file
          */
         MHDFloat mTime;

         /**
          * @brief Timestep read from file
          */
         MHDFloat mTimestep;

         /**
          * @brief Read run information to file
          */
         void readRun();

         /**
          * @brief Read scalar field values from file
          *
          * @param name    Name of the scalar field
          * @param rScalar Storage for the field
          */
         void readSpectralScalar(const std::string& name, Datatypes::SpectralScalarType& rScalar);

         /**
          * @brief Read vector field values from file
          *
          * @param name    Name of the vector field
          * @param rVector Storage for the field
          */
         void readSpectralVector(const std::string& name, std::vector<Datatypes::SpectralScalarType>& rVector);

         /**
          * @brief Read vector field component from file
          *
          * @param name    Name of the vector field
          * @param id      ID of the component
          * @param rComp   Storage for field component
          */
         void readSpectralComponent(const std::string& name, FieldComponents::Spectral::Id id, Datatypes::SpectralScalarType& rComp);

      private:
   };

   inline MHDFloat StateFileReader::time() const
   {
      return this->mTime;
   }

   inline MHDFloat StateFileReader::timestep() const
   {
      return this->mTimestep;
   }

}
}

#endif // STATEFILEREADER_HPP
