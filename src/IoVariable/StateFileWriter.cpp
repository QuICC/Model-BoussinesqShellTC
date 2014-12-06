/** 
 * @file StateFileWriter.cpp
 * @brief Source of the implementation of the HDF5 state file writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
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
      this->writeRun();

      // Write all the scalars
      StateFileWriter::scalar_iterator_range sRange = this->scalarRange();
      StateFileWriter::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         this->writeSpectralScalar(IoTools::IdToHuman::toTag(sit->first), sit->second->dom(0).perturbation());
      }

      // Write all the vectors
      StateFileWriter::vector_iterator_range vRange = this->vectorRange();
      StateFileWriter::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         this->writeSpectralVector(IoTools::IdToHuman::toTag(vit->first), vit->second->dom(0).perturbation().data());
      }

      // Close file
      this->postWrite();
   }

   void StateFileWriter::writeSpectralScalar(const std::string& name, const Datatypes::SpectralScalarType& scalar)
   {
      // Create the scalar group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, const Datatypes::SpectralScalarType::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(scalar);

      // Check for data regularity
      if(this->mIsRegular)
      {
         // Write the scalar values
         this->writeRegularField(group, name, fieldInfo);
      } else
      {
         // Write the scalar values
         this->writeIrregularField(group, name, fieldInfo);
      }
      
      // close group
      H5Gclose(group);
   }

   void StateFileWriter::writeSpectralVector(const std::string& name, const std::map<FieldComponents::Spectral::Id,Datatypes::SpectralScalarType>& vector)
   {
      // Create the vector field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, const Datatypes::SpectralScalarType::PointType *> > fieldInfo;

      // Check for data regularity
      std::map<FieldComponents::Spectral::Id,Datatypes::SpectralScalarType>::const_iterator it;
      if(this->mIsRegular)
      {
         for(it = vector.begin(); it != vector.end(); ++it)
         {
            // create component field information
            fieldInfo = Datatypes::FieldTools::createInfo(it->second);

            // Write vector component
            this->writeRegularField(group, name+"_"+IoTools::IdToHuman::toTag(it->first), fieldInfo);
         }
      } else
      {
         for(it = vector.begin(); it != vector.end(); ++it)
         {
            // create component field information
            fieldInfo = Datatypes::FieldTools::createInfo(it->second);

            // Write vector component
            this->writeIrregularField(group, name+"_"+IoTools::IdToHuman::toTag(it->first), fieldInfo);
         }
      }
      
      // close group
      H5Gclose(group);
   }

}
}
