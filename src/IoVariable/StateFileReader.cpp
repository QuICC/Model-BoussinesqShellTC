/** 
 * @file StateFileReader.cpp
 * @brief Source of the implementation of HDF5 state file reader
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
#include "IoVariable/StateFileReader.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/StateFileTags.hpp"
#include "IoVariable/VariableHdf5Tags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   StateFileReader::StateFileReader(std::string name, std::string type, const bool isRegular)
      : IVariableHdf5Reader(StateFileTags::BASENAME + name, StateFileTags::EXTENSION, StateFileTags::HEADER, type, StateFileTags::VERSION, Dimensions::Space::TRANSFORM, isRegular), mTime(-1.0), mTimestep(-1.0)
   {
   }

   StateFileReader::~StateFileReader()
   {
   }

   void StateFileReader::read()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read the run information
      this->readRun();

      // Read all the scalars
      StateFileReader::scalar_iterator_range sRange = this->scalarRange();
      StateFileReader::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         // Make sure full field is zero
         sit->second->setZeros();

         // Read field values
         this->readSpectralScalar(IoTools::IdToHuman::toTag(sit->first), sit->second->rDom(0).rPerturbation());
      }

      // Read all the vectors
      StateFileReader::vector_iterator_range vRange = this->vectorRange();
      StateFileReader::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         // Make sure full field is zero
         vit->second->setZeros();

         // Read field values
         this->readSpectralVector(IoTools::IdToHuman::toTag(vit->first), vit->second->rDom(0).rPerturbation().rData());
      }
   }

   void StateFileReader::readSetup()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation();

      // Set Read arguments
      this->setReadArguments();

      // Read the run information
      this->readRun();
   }

   void StateFileReader::readRun()
   {
      // Open the run paramters group
      hid_t group = H5Gopen(this->file(), VariableHdf5Tags::RUN.c_str(), H5P_DEFAULT);

      // Read the reached simulation time from file
      this->readScalar(group, VariableHdf5Tags::RUNTIME, this->mTime);

      // Read the last used timestep from file
      this->readScalar(group, VariableHdf5Tags::RUNSTEP, this->mTimestep);
      
      // close group
      H5Gclose(group);
   }

   void StateFileReader::readSpectralScalar(const std::string& name, Datatypes::SpectralScalarType& rScalar)
   {
      // Open the codensity scalar group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, Datatypes::SpectralScalarType::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(rScalar);

      // Check for data regularity
      if(this->mIsRegular)
      {
         this->readRegularField(group, name, fieldInfo);
      } else
      {
         this->readIrregularField(group, name, fieldInfo);
      }
      
      // close group
      H5Gclose(group);
   }

   void StateFileReader::readSpectralVector(const std::string& name, std::vector<Datatypes::SpectralScalarType>& rVector)
   {
      // Open the magnetic field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, Datatypes::SpectralScalarType::PointType *> > fieldInfo;

      // Check for data regularity
      if(this->mIsRegular)
      {
         for(size_t i = 0; i < rVector.size(); i++)
         { 
            // create component field information
            fieldInfo = Datatypes::FieldTools::createInfo(rVector.at(i));

            // Read component from file 
            this->readRegularField(group,name+"_"+IoTools::IdToHuman::toTag(static_cast<FieldComponents::Spectral::Id>(i)), fieldInfo);
         }
      } else
      {
         for(size_t i = 0; i < rVector.size(); i++)
         {
            // create component field information
            fieldInfo = Datatypes::FieldTools::createInfo(rVector.at(i));

            // Read component from file 
            this->readIrregularField(group,name+"_"+IoTools::IdToHuman::toTag(static_cast<FieldComponents::Spectral::Id>(i)), fieldInfo);
         }
      }
      
      // close group
      H5Gclose(group);
   }

   void StateFileReader::readSpectralComponent(const std::string& name, FieldComponents::Spectral::Id id, Datatypes::SpectralScalarType& rComp)
   {
      // Open the vector field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, Datatypes::SpectralScalarType::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(rComp);

      // Check for data regularity
      if(this->mIsRegular)
      {
         // Read the field component
         this->readRegularField(group, name+"_"+IoTools::IdToHuman::toTag(id), fieldInfo);
      } else
      {
         // Read the field component
         this->readIrregularField(group, name+"_"+IoTools::IdToHuman::toTag(id), fieldInfo);
      }
      
      // close group
      H5Gclose(group);
   }

}
}
