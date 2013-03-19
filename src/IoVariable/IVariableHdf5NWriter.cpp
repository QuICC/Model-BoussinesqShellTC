/** \file IVariableHdf5NWriter.cpp 
 *  \brief Source of the implementation of the generic variable to file writer
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/IVariableHdf5NWriter.hpp"

// Project includes
//
#include "IoVariable/VariableHdf5Tags.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   IVariableHdf5NWriter::IVariableHdf5NWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const bool isRegular)
      : IHdf5NWriter(name, ext, header, type, version), mTime(-1.0), mTimestep(-1.0), mIsRegular(isRegular), mSpaceId(id)
   {
   }

   IVariableHdf5NWriter::~IVariableHdf5NWriter()
   {
   }

   Dimensions::Space::Id IVariableHdf5NWriter::space() const
   {
      return this->mSpaceId;
   }

   void IVariableHdf5NWriter::setPhysical(const std::map<std::string,MHDFloat>& parameters)
   {
      this->mPhysical = parameters;
   }

   void IVariableHdf5NWriter::setMesh(const std::vector<Array>& mesh)
   {
      this->mMesh = mesh;
   }

   void IVariableHdf5NWriter::setSimTime(const MHDFloat time, const MHDFloat timestep)
   {
      this->mTime = time;

      this->mTimestep = timestep;
   }

   void IVariableHdf5NWriter::setResolution(SharedResolution spRes)
   {
      // Store resolution object
      this->mspRes = spRes;

      // Set dataset dimensions
      this->setDatasetSize();

      // Set dataset offsets
      this->setDatasetOffsets();

      // This set collective IO operations
      this->setCollIo();
   }

   void IVariableHdf5NWriter::expect(const PhysicalNames::Id id)
   {
      this->mExpected.insert(id);
   }

   bool IVariableHdf5NWriter::isFull() const
   {
      bool status = true;

      // Check that all expected scalars and vectors are present
      status = status && (this->mScalars.size() + this->mVectors.size() == this->mExpected.size());

      // Check that the resolution has been set
      status = status && this->mspRes;

      if(this->mspRes && (this->mSpaceId == Dimensions::Space::PHYSICAL))
      {
         status = status && (this->mMesh.size() == static_cast<size_t>(this->mspRes->cpu()->nDim()));
      }

      return status;
   }

   void IVariableHdf5NWriter::addScalar(const std::pair<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalar)
   {
      // Only add the variable if it is listed as expected
      if(this->mExpected.count(scalar.first))
      {
         this->mScalars.insert(scalar);

         // Resolution is not yet set extract from scalar
         if(!this->mspRes)
         {
            this->setResolution(scalar.second->dom(0).spRes());
         }
      }
   }
 
   void IVariableHdf5NWriter::addVector(const std::pair<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vector)
   {
      // Only add the variable if it is listed as expected
      if(this->mExpected.count(vector.first))
      {
         this->mVectors.insert(vector);

         // Resolution is not yet set extract from vector
         if(!this->mspRes)
         {
            this->setResolution(vector.second->dom(0).spRes());
         }
      }
   }

   IVariableHdf5NWriter::scalar_iterator_range  IVariableHdf5NWriter::scalarRange()
   {
      return std::make_pair(this->mScalars.begin(), this->mScalars.end());
   }

   IVariableHdf5NWriter::vector_iterator_range  IVariableHdf5NWriter::vectorRange()
   {
      return std::make_pair(this->mVectors.begin(), this->mVectors.end());
   }

   void IVariableHdf5NWriter::setDatasetSize()
   {
      if(this->mIsRegular)
      {
         // Get dimensions ordered by index access speed (fast -> slow)
         ArrayI oDims = this->mspRes->sim()->orderedDims(this->mSpaceId);

         int nDims = oDims.size();
         for(int i = 0; i < nDims; ++i)
         {
            // Set the dimension size in reverse order for HDF5
            this->mFileDims.push_back(oDims(nDims-1-i));
         }
      } else
      {
         /// \mhdBug Irregular grid is not working, requires modification of simulation resolution
         throw Exception("Irregular data is not yet implemented in HDF5 storage");
//         // Set the slow dimension size
//IS WRONG         this->mFileDims.push_back(this->mspRes->sim()->nSlow(this->mSpaceId));
//
//         // Set the second dimension size (fastest in C ordering)
//IS WRONG         this->mFileDims.push_back(this->mspRes->sim()->nFast(this->mSpaceId));
      }
   }

   void IVariableHdf5NWriter::setDatasetOffsets()
   {
      Dimensions::Transform::Id  transId;

      // Select transform dimension depending on dimension space
      if(this->mSpaceId == Dimensions::Space::SPECTRAL)
      {
         transId = Dimensions::Transform::TRA1D;
      } else //if(this->mSpaceId == Dimensions::Space::PHYSICAL)
      {
         transId = Dimensions::Transform::TRA3D;
      }

      if(this->mIsRegular)
      {
         // Loop over the stored indexes
         std::vector<hsize_t>  offV;
         offV.push_back(0);
         offV.push_back(0);
         for(int i=0; i < this->mspRes->cpu()->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++i)
         {
            // Compute offset for third dimension
            offV.at(0) = this->mspRes->cpu()->dim(transId)->idx<Dimensions::Data::DAT3D>(i);

            // Compute offset for second dimension
            offV.at(1) = this->mspRes->cpu()->dim(transId)->idx<Dimensions::Data::DAT2D>(0,i);

            // Add offset to vector
            this->mFileOffsets.push_back(offV);
         }
      } else
      {
         /// \mhdBug Irregular grid is not working, requires modification of simulation resolution
         throw Exception("Irregular data is not yet implemented in HDF5 storage");

         // offset HDF5 type
         hsize_t  offset;

         // Loop over the stored indexes
         std::vector<hsize_t>  offV;
         for(int i=0; i < this->mspRes->cpu()->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++i)
         {
            // Compute offset from previous index
            offset = 0;
            for(int j = 0; j < this->mspRes->cpu()->dim(transId)->idx<Dimensions::Data::DAT3D>(i); ++j)
            {
//IS WRONG              offset += this->mspRes->sim()->dim(Dimensions::Simulation::SIM3D, this->mSpaceId);
            }

            // Compute the offset for the local indexes
            for(int j = 0; j < this->mspRes->cpu()->dim(transId)->dim<Dimensions::Data::DAT2D>(); ++j)
            {
//IS WRONG               offV.push_back(offset + this->mspRes->cpu()->dim(transId)->idx<Dimensions::Data::DAT2D>(j,i));
            }

            // Add offset to vector
            this->mFileOffsets.push_back(offV);
            offV.clear();
         }
      }
   }

   void IVariableHdf5NWriter::setCollIo()
   {
      this->mCollIoWrite = 0;

      // Select transform dimension depending on dimension space
      Dimensions::Transform::Id  transId;
      if(this->mSpaceId == Dimensions::Space::SPECTRAL)
      {
         transId = Dimensions::Transform::TRA1D;
      } else
      {
         transId = Dimensions::Transform::TRA3D;
      }

      // Get the maximum number of slowest directions over all CPUs
      for(int i = 0; i < this->mspRes->nCpu(); ++i)
      {
         if(this->mCollIoWrite < this->mspRes->cpu(i)->dim(transId)->dim<Dimensions::Data::DAT3D>())
         {
            this->mCollIoWrite = this->mspRes->cpu(i)->dim(transId)->dim<Dimensions::Data::DAT3D>();
         }
      }
   }

   void IVariableHdf5NWriter::writeTruncation()
   {
      std::ostringstream   oss;

      // Create the truncation parameters group
      hid_t base = H5Gcreate(this->file(), VariableHdf5Tags::TRUNCATION.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Create the spectral truncation subgroup
      hid_t subGroup = H5Gcreate(base, VariableHdf5Tags::TRUNCSPECTRAL.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write spectral resolution information
      for(int i = 0; i < this->mspRes->cpu()->nDim(); i++)
      {
         oss << VariableHdf5Tags::TRUNCDIM << i+1 << "D";

         // Write truncation value to file
         this->writeScalar(subGroup, oss.str(), this->mspRes->sim()->dim(static_cast<Dimensions::Simulation::Id>(i),Dimensions::Space::SPECTRAL)-1);

         oss.str("");
      }
      
      // close group
      H5Gclose(subGroup);

      // Create the physical truncation subgroup
      subGroup = H5Gcreate(base, VariableHdf5Tags::TRUNCPHYSICAL.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write physical resolution information
      for(int i = 0; i < this->mspRes->cpu()->nDim(); i++)
      {
         oss << VariableHdf5Tags::TRUNCDIM << i+1 << "D";

         // Write truncation value to file
         this->writeScalar(subGroup, oss.str(), this->mspRes->sim()->dim(static_cast<Dimensions::Simulation::Id>(i),Dimensions::Space::PHYSICAL)-1);

         oss.str("");
      }
      
      // close group
      H5Gclose(subGroup);
      
      // close group
      H5Gclose(base);
   }

   void IVariableHdf5NWriter::writePhysical()
   {
      // Create the Physical parameters group
      hid_t group = H5Gcreate(this->file(), VariableHdf5Tags::PHYSICAL.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      std::map<std::string,MHDFloat>::const_iterator it;
      for(it = this->mPhysical.begin(); it != this->mPhysical.end(); ++it)
      {
         // Write reached simulation time to file
         this->writeScalar(group, it->first, it->second);
      }
      
      // close group
      H5Gclose(group);
   }

   void IVariableHdf5NWriter::writeRun()
   {
      // Create the Run parameters group
      hid_t group = H5Gcreate(this->file(), VariableHdf5Tags::RUN.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write reached simulation time to file
      this->writeScalar(group, VariableHdf5Tags::RUNTIME, this->mTime);

      // Write last timestep to file
      this->writeScalar(group, VariableHdf5Tags::RUNSTEP, this->mTimestep);
      
      // close group
      H5Gclose(group);
   }

}
}
