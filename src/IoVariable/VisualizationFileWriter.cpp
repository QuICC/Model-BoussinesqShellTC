/** \file VisualizationFileWriter.cpp
 *  \brief Source of the implementation of the HDF5 visualisation file writer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/VisualizationFileWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldComponents.hpp"
#include "IoVariable/VisualizationFileTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   VisualizationFileWriter::VisualizationFileWriter(std::string type, const std::vector<Array>& mesh)
      : IVariableHdf5NWriter(VisualizationFileTags::BASENAME, VisualizationFileTags::EXTENSION, VisualizationFileTags::HEADER, type, VisualizationFileTags::VERSION, Dimensions::Space::PHYSICAL, true), mMesh(mesh)
   {
   }

   VisualizationFileWriter::~VisualizationFileWriter()
   {
   }

   void VisualizationFileWriter::write()
   {
      // Create file
      this->preWrite();

      // Create the header and version information
      this->createFileInfo();

      // Write the Physical parameters
      this->writeMesh();

      // Write the truncation information
      this->writeTruncation();

      // Write the run information
      this->writeRun(0.0, 0.0);

      // Write all the scalars
      VisualizationFileWriter::scalar_iterator_range sRange = this->scalarRange();
      VisualizationFileWriter::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         this->writePhysicalScalar(IoTools::IdToHuman::toString(sit->first), sit->second->dom(0).phys());
      }

      // Write all the vectors
      VisualizationFileWriter::vector_iterator_range vRange = this->vectorRange();
      VisualizationFileWriter::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         this->writePhysicalVector(IoTools::IdToHuman::toString(vit->first), vit->second->dom(0).phys().data());
      }

      // Close file
      this->postWrite();
   }

   void VisualizationFileWriter::writeMesh()
   {
      // Create the Physical parameters group
      hid_t group = H5Gcreate(this->file(), VisualizationFileTags::MESH.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for(int i = 0; i < this->mspRes->cpu()->nDim(); i ++)
      {
         this->writeArray(group, VisualizationFileTags::GRID+"_"+IoTools::IdToHuman::toTag(static_cast<FieldComponents::Physical::Id>(i)), this->mMesh.at(i));
      }
      
      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writeRun(const MHDFloat time, const MHDFloat step)
   {
      // Create the Run parameters group
      hid_t group = H5Gcreate(this->file(), VisualizationFileTags::RUN.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write reached simulation time to file
      this->writeScalar(group, VisualizationFileTags::RUNTIME, time);

      // Write last timestep to file
      this->writeScalar(group, VisualizationFileTags::RUNSTEP, step);
      
      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalScalar(const std::string& name, const Datatypes::PhysicalScalarType& scalar)
   {
      // Create the Codensity scalar group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the scalar field
      /// \mhdBug need to restructure general HDF5 data writer
      //this->writeRegularField(group, name, scalar.sliced());
      
      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalVector(const std::string& name, const std::vector<Datatypes::PhysicalScalarType>& vector)
   {
      // Create the Magnetic Field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /// \mhdBug need to restructure general HDF5 data writer
      //for(int i = 0; i < vector.size(); i++)
      //{
      //   // Write the vector field
      //   this->writeRegularField(group, name+"_"+IoTools::IdToHuman::toTag(static_cast<FieldComponents::Physical::Id>(i)), vector.at(i).sliced());
      //}
      
      // close group
      H5Gclose(group);
   }

}
}
