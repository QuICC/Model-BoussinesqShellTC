/** \file StateFileReader.cpp
 *  \brief Source of the implementation of HDF5 state file reader
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
#include "Enums/FieldComponents.hpp"
#include "IoVariable/StateFileTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   StateFileReader::StateFileReader(std::string name, std::string type)
      : IVariableHdf5Reader(StateFileTags::BASENAME + name, StateFileTags::EXTENSION, StateFileTags::HEADER, type, StateFileTags::VERSION), mTime(-1.0), mTimestep(-1.0)
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
      this->checkTruncation(Dimensions::Space::SPECTRAL);

      // Set Read arguments
      this->setReadArguments(Dimensions::Space::SPECTRAL);

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
         this->readSpectralScalar(IoTools::IdToHuman::toString(sit->first), sit->second->rDom(0).rPerturbation());
      }

      // Read all the vectors
      StateFileReader::vector_iterator_range vRange = this->vectorRange();
      StateFileReader::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         // Make sure full field is zero
         vit->second->setZeros();

         // Read field values
         this->readSpectralVector(IoTools::IdToHuman::toString(vit->first), vit->second->rDom(0).rPerturbation().rData());
      }
   }

   void StateFileReader::readSetup()
   {
      // Read the truncation information
      this->readTruncation();

      // Check file compatibility with data truncation
      this->checkTruncation(Dimensions::Space::SPECTRAL);

      // Set Read arguments
      this->setReadArguments(Dimensions::Space::SPECTRAL);

      // Read the run information
      this->readRun();
   }

   void StateFileReader::readRun()
   {
      // Open the run paramters group
      hid_t group = H5Gopen(this->file(), StateFileTags::RUN.c_str(), H5P_DEFAULT);

      // Read the reached simulation time from file
      this->readScalar(group, StateFileTags::RUNTIME, this->mTime);

      // Read the last used timestep from file
      this->readScalar(group, StateFileTags::RUNSTEP, this->mTimestep);
      
      // close group
      H5Gclose(group);
   }

   void StateFileReader::readSpectralScalar(const std::string& name, Datatypes::SpectralScalarType& rScalar)
   {
      // Open the codensity scalar group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      /// \mhdBug Irregular grid writing is not supported yet
      // Read the codensity expansion
      //this->readIrregularField(group, name, rScalar.rSliced());
      
      // close group
      H5Gclose(group);
   }

   void StateFileReader::readSpectralVector(const std::string& name, std::vector<Datatypes::SpectralScalarType>& rVector)
   {
      // Open the magnetic field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      /// \mhdBug Irregular grid writing is not supported yet
      //for(int i = 0; i < rVector.size(); i++)
      //{
      //   // Read component from file 
      //   this->readIrregularField(group,name+"_"+IoTools::IdToHuman::toTag(static_cast<FieldComponents::Spectral::Component>(i)), rVector.at(i).rSliced());
      //}
      
      // close group
      H5Gclose(group);
   }

   void StateFileReader::readSpectralComponent(const std::string& name, FieldComponents::Spectral::Id id, Datatypes::SpectralScalarType& rComp)
   {
      // Open the vector field group
      hid_t group = H5Gopen(this->file(), name.c_str(), H5P_DEFAULT);

      /// \mhdBug Irregular grid writing is not supported yet
      // Read the field component
      //this->readIrregularField(group, name+"_"+IoTools::IdToHuman::toTag(id), rComp.rSliced());
      
      // close group
      H5Gclose(group);
   }

}
}
