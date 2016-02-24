/** 
 * @file VisualizationFileWriter.cpp
 * @brief Source of the implementation of the HDF5 visualisation file writer
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
#include "IoVariable/VisualizationFileWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/VisualizationFileTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   VisualizationFileWriter::VisualizationFileWriter(std::string type)
      : IVariableHdf5NWriter(VisualizationFileTags::BASENAME, VisualizationFileTags::EXTENSION, VisualizationFileTags::HEADER, type, VisualizationFileTags::VERSION, Dimensions::Space::PHYSICAL, true)
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
      this->writePhysical();

      // Write the Physical parameters
      this->writeMesh();

      // Write the truncation information
      this->writeTruncation();

      // Write the run information
      this->writeRun();

      // Write all the scalars
      VisualizationFileWriter::scalar_iterator_range sRange = this->scalarRange();
      VisualizationFileWriter::scalar_iterator sit;
      for(sit = sRange.first; sit != sRange.second; ++sit)
      {
         if(sit->second->dom(0).hasPhys())
         {
            this->writePhysicalScalar(IoTools::IdToHuman::toTag(sit->first), sit->second->dom(0).phys());
         }

         if(sit->second->dom(0).hasGrad())
         {
            this->writePhysicalVector(IoTools::IdToHuman::toTag(sit->first)+"_grad", sit->second->dom(0).grad().data());
         }

         if(sit->second->dom(0).hasGrad2())
         {
            this->writePhysicalTensor(IoTools::IdToHuman::toTag(sit->first)+"_grad2", sit->second->dom(0).grad2().data());
         }
      }

      // Write all the vectors
      VisualizationFileWriter::vector_iterator_range vRange = this->vectorRange();
      VisualizationFileWriter::vector_iterator vit;
      for(vit = vRange.first; vit != vRange.second; ++vit)
      {
         if(vit->second->dom(0).hasPhys())
         {
            this->writePhysicalVector(IoTools::IdToHuman::toTag(vit->first), vit->second->dom(0).phys().data());
         }

         if(vit->second->dom(0).hasGrad())
         {
            std::vector<FieldComponents::Spectral::Id> fId;
            fId.push_back(FieldComponents::Spectral::ONE);
            fId.push_back(FieldComponents::Spectral::TWO);
            fId.push_back(FieldComponents::Spectral::THREE);
            for(unsigned int i = 0; i < fId.size(); i ++)
            {
               if(fId.at(i) != FieldComponents::Spectral::NOTUSED)
               {
                  this->writePhysicalVector(IoTools::IdToHuman::toTag(vit->first) + "_" + IoTools::IdToHuman::toTag(fId.at(i)) + "_grad", vit->second->dom(0).grad(fId.at(i)).data());
               }
            }
         }

         if(vit->second->dom(0).hasCurl())
         {
            this->writePhysicalVector(IoTools::IdToHuman::toTag(vit->first)+"_curl", vit->second->dom(0).curl().data());
         }
      }

      // Close file
      this->postWrite();
   }

   void VisualizationFileWriter::writeMesh()
   {
      // Create the Physical parameters group
      hid_t group = H5Gcreate(this->file(), VisualizationFileTags::MESH.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      std::vector<FieldComponents::Physical::Id> gridId;
      gridId.push_back(FieldComponents::Physical::ONE);
      gridId.push_back(FieldComponents::Physical::TWO);
      gridId.push_back(FieldComponents::Physical::THREE);
      for(int i = 0; i < this->mspRes->cpu()->nDim(); i ++)
      {
         if(gridId.at(i) != FieldComponents::Physical::NOTUSED)
         {
            this->writeArray(group, VisualizationFileTags::GRID+"_"+IoTools::IdToHuman::toTag(gridId.at(i)), this->mMesh.at(i));
         }
      }
      
      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalScalar(const std::string& name, const Datatypes::PhysicalScalarType& scalar)
   {
      // Create the scalar group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, const Datatypes::PhysicalScalarType::PointType *> > fieldInfo = Datatypes::FieldTools::createInfo(scalar);

      // Write the scalar field
      this->writeRegularField(group, name, fieldInfo);
      
      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalVector(const std::string& name, const std::map<FieldComponents::Physical::Id,Datatypes::PhysicalScalarType>& vector, const std::string& joint)
   {
      // Create the Field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, const Datatypes::PhysicalScalarType::PointType *> > fieldInfo;

      std::map<FieldComponents::Physical::Id,Datatypes::PhysicalScalarType>::const_iterator it;
      for(it = vector.begin(); it != vector.end(); ++it)
      {
         // create component field information
         fieldInfo = Datatypes::FieldTools::createInfo(it->second);

         // Write the vector field
         this->writeRegularField(group, name + joint + IoTools::IdToHuman::toTag(it->first), fieldInfo);
      }
      
      // close group
      H5Gclose(group);
   }

   void VisualizationFileWriter::writePhysicalTensor(const std::string& name, const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,Datatypes::PhysicalScalarType>& tensor, const std::string& joint)
   {
      // Create the Field group
      hid_t group = H5Gcreate(this->file(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Storage for the field information
      std::vector<std::tr1::tuple<int,int, const Datatypes::PhysicalScalarType::PointType *> > fieldInfo;

      std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,Datatypes::PhysicalScalarType>::const_iterator it;
      for(it = tensor.begin(); it != tensor.end(); ++it)
      {
         // create component field information
         fieldInfo = Datatypes::FieldTools::createInfo(it->second);

         // Write the tensor field
         this->writeRegularField(group, name + joint + IoTools::IdToHuman::toTag(it->first.first) + IoTools::IdToHuman::toTag(it->first.second), fieldInfo);
      }
      
      // close group
      H5Gclose(group);
   }

}
}
