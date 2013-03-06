/** \file StateFileWriter.cpp
 *  \brief Source of the implementation of the HDF5 state file writer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/StateFileWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldComponents.hpp"
#include "IoVariable/StateFileTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   StateFileWriter::StateFileWriter(std::string type, const bool isRegular)
      : IVariableHdf5NWriter(StateFileTags::BASENAME, StateFileTags::EXTENSION, StateFileTags::HEADER, type, StateFileTags::VERSION, Dimensions::Space::SPECTRAL, isRegular)
   {
   }

   StateFileWriter::~StateFileWriter()
   {
   }

   void StateFileWriter::write()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Write the Physical parameters
      this->writePhysical();

      // Write the truncation information
      this->writeTruncation();

      // Write the run information
      this->writeRun(0.0, 0.0);

      // Write all the scalars
      StateFileWriter::scalar_iterator_range sRange = this->scalarRange();
      StateFileWriter::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         this->writeSpectralScalar(IoTools::IdToHuman::toString(sit->first), sit->second->dom(0).perturbation());
      }

      // Write all the vectors
      StateFileWriter::vector_iterator_range vRange = this->vectorRange();
      StateFileWriter::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         this->writeSpectralVector(IoTools::IdToHuman::toString(vit->first), vit->second->dom(0).perturbation().data());
      }

      // Close file
      this->postWrite();
   }

   void StateFileWriter::writePhysical()
   {
      // Create the Physical parameters group
      hid_t group = H5Gcreate(this->file(), StateFileTags::PHYSICAL.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /// \mhdBug This has not yet been implemented!
      
      // close group
      H5Gclose(group);
   }

   void StateFileWriter::writeRun(const MHDFloat time, const MHDFloat step)
   {
      // Create the Run parameters group
      hid_t group = H5Gcreate(this->file(), StateFileTags::RUN.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write reached simulation time to file
      this->writeScalar(group, StateFileTags::RUNTIME, time);

      // Write last timestep to file
      this->writeScalar(group, StateFileTags::RUNSTEP, step);
      
      // close group
      H5Gclose(group);
   }

   void StateFileWriter::writeSpectralScalar(const std::string& name, const Datatypes::SpectralScalarType& scalar)
   {
      // Create the scalar group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Check for data regularity
      if(this->mIsRegular)
      {
         /// \mhdBug Need to to produce a "sliced" data type out of the flat storage
         // Write the scalar values
         //this->writeRegularField(group, name, scalar.sliced());
      } else
      {
         /// \mhdBug Irregular grid is not working, requires modification of simulation resolution
         throw Exception("Irregular data is not yet implemented in HDF5 storage");
         // Write the scalar values
         //this->writeIrregularField(group, name, scalar.sliced());
      }
      
      // close group
      H5Gclose(group);
   }

   void StateFileWriter::writeSpectralVector(const std::string& name, const std::vector<Datatypes::SpectralScalarType>& vector)
   {
      // Create the vector field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Check for data regularity
      if(this->mIsRegular)
      {
         for(size_t i = 0; i < vector.size(); i++)
         {
            /// \mhdBug Need to to produce a "sliced" data type out of the flat storage
            // Write vector component
            //this->writeRegularField(group, name+"_"+IoTools::IdToHuman::toTag(static_cast<FieldComponents::Spectral::Id>(i)), vector.at(i).sliced());
         }
      } else
      {
         /// \mhdBug Irregular grid is not working, requires modification of simulation resolution
         throw Exception("Irregular data is not yet implemented in HDF5 storage");
         for(size_t i = 0; i < vector.size(); i++)
         {
            // Write vector component
            //this->writeIrregularField(group, name+"_"+IoTools::IdToHuman::toTag(static_cast<FieldComponents::Spectral::Id>(i)), vector.at(i).sliced());
         }
      }
      
      // close group
      H5Gclose(group);
   }

}
}
