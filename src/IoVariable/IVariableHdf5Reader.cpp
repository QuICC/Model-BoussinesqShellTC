/** \file IVariableHdf5Reader.cpp
 *  \brief Source of the implementation of a generic variable data file reader
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
#include "IoVariable/IVariableHdf5Reader.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoVariable/VariableHdf5Tags.hpp"
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   IVariableHdf5Reader::IVariableHdf5Reader(std::string name, std::string ext, std::string header, std::string type, std::string version, const bool isRegular)
      : IHdf5Reader(name, ext, header, type, version), mIsRegular(isRegular)
   {
   }

   IVariableHdf5Reader::~IVariableHdf5Reader()
   {
   }

   void IVariableHdf5Reader::setResolution(SharedResolution spRes)
   {
      // Store resolution object
      this->mspRes = spRes;

      // This set collective IO operations
      this->setCollIo();
   }

   void IVariableHdf5Reader::expect(const PhysicalNames::Id id)
   {
      this->mExpected.insert(id);
   }

   bool IVariableHdf5Reader::isFull() const
   {
      // Check that all expected scalars and vectors are present
      bool sizeStatus = (this->mScalars.size() + this->mVectors.size() == this->mExpected.size());

      // Check that the resolution has been set
      bool resStatus = this->mspRes;

      return (sizeStatus && resStatus);
   }

   void IVariableHdf5Reader::addScalar(const std::pair<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalar)
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
 
   void IVariableHdf5Reader::addVector(const std::pair<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vector)
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

   IVariableHdf5Reader::scalar_iterator_range  IVariableHdf5Reader::scalarRange()
   {
      return std::make_pair(this->mScalars.begin(), this->mScalars.end());
   }

   IVariableHdf5Reader::vector_iterator_range  IVariableHdf5Reader::vectorRange()
   {
      return std::make_pair(this->mVectors.begin(), this->mVectors.end());
   }

   void IVariableHdf5Reader::setCollIo()
   {
      this->mCollIoRead = 0;
   }

   void IVariableHdf5Reader::setReadArguments(const Dimensions::Space::Id id)
   {
      /// \mhdBug Code assumes regular data

      Dimensions::Transform::Id  transId;

      // Select transform dimension depending on dimension space
      if(id == Dimensions::Space::SPECTRAL)
      {
         transId = Dimensions::Transform::TRA1D;
      } else
      {
         transId = Dimensions::Transform::TRA3D;
      }

      // Set block size (i.e. number of contiguous elements)
      this->mBlock = std::min(this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,id), this->mspFileRes->dim(Dimensions::Simulation::SIM1D,id));

      // offset HDF5 type
      hsize_t offset;

      // Loop over the stored indexes
      std::vector<hsize_t>  offV;
      for(int i=0; i < this->mspRes->cpu()->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         int idx3D = this->mspRes->cpu()->dim(transId)->idx<Dimensions::Data::DAT3D>(i);
         // Check if value is available in file
         if(idx3D < this->mspFileRes->dim(Dimensions::Simulation::SIM3D,id))
         {
            // Compute offset from previous index
            offset = 0;
            for(int j = 0; j < idx3D; ++j)
            {
               offset += this->mspFileRes->dim(Dimensions::Simulation::SIM2D,id);
            }

            // Compute the offset for the local indexes
            for(int j = 0; j < this->mspRes->cpu()->dim(transId)->dim<Dimensions::Data::DAT2D>(i); ++j)
            {
               if(this->mspRes->cpu()->dim(transId)->idx<Dimensions::Data::DAT2D>(j,i) < this->mspFileRes->dim(Dimensions::Simulation::SIM2D,id))
               {
                  offV.push_back(offset + this->mspRes->cpu()->dim(transId)->idx<Dimensions::Data::DAT2D>(j,i));
               }
            }

            // Add offset to vector
            this->mFileOffsets.push_back(offV);
            offV.clear();
         }
      }

      // Get the "global" local minimum for MPI code
      this->mCollIoRead = this->mFileOffsets.size();
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mCollIoRead, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      #endif // GEOMHDISCC_MPI
   }

   void IVariableHdf5Reader::readTruncation()
   {
      std::ostringstream   oss;

      // Open the truncation parameters group
      hid_t base = H5Gopen(this->file(), VariableHdf5Tags::TRUNCATION.c_str(), H5P_DEFAULT);

      // Open the spectral truncation subgroup
      hid_t subGroup = H5Gopen(base, VariableHdf5Tags::TRUNCSPECTRAL.c_str(), H5P_DEFAULT);

      ArrayI spec(this->mspRes->cpu()->nDim());

      // Read spectral resolution information
      for(int i = 0; i < this->mspRes->cpu()->nDim(); i++)
      {
         oss << VariableHdf5Tags::TRUNCDIM << i+1 << "D";

         // Read dimension from file
         this->readScalar(subGroup, oss.str(), spec(i));

         oss.str("");
      }
      spec.array() += 1;
      
      // close group
      H5Gclose(subGroup);

      // Open the physical truncation parameters group
      subGroup = H5Gopen(base, VariableHdf5Tags::TRUNCPHYSICAL.c_str(), H5P_DEFAULT);

      ArrayI phys(this->mspRes->cpu()->nDim());

      // Read spectral resolution information
      for(int i = 0; i < this->mspRes->cpu()->nDim(); i++)
      {
         oss << VariableHdf5Tags::TRUNCDIM << i+1 << "D";

         // Read dimension from file
         this->readScalar(subGroup, oss.str(), phys(i));

         oss.str("");
      }
      phys.array() += 1;
      
      // close group
      H5Gclose(subGroup);
      
      // close group
      H5Gclose(base);

      this->mspFileRes = SharedSimulationResolution(new SimulationResolution(phys,spec));
   }

   void IVariableHdf5Reader::checkTruncation(Dimensions::Space::Id id)
   {
      std::ostringstream   oss;

      for(int i = 0; i < this->mspRes->cpu()->nDim(); i++)
      {
         if(this->mspRes->sim()->dim(static_cast<Dimensions::Simulation::Id>(i),id) != this->mspFileRes->dim(static_cast<Dimensions::Simulation::Id>(i),id))
         {
            oss << "Dimension " << i+1 <<  " doesn't fit";
            IoTools::Formatter::printLine(std::cout, '-');
            IoTools::Formatter::printCentered(std::cout, oss.str(), '*');

            if(this->mspRes->sim()->dim(static_cast<Dimensions::Simulation::Id>(i),id) > this->mspFileRes->dim(static_cast<Dimensions::Simulation::Id>(i),id))
            {
               IoTools::Formatter::printCentered(std::cout, " ---> Zeros have been added!", ' ');
            } else
            {
               IoTools::Formatter::printCentered(std::cout, "---> File data has been truncated!", ' ');
            }
            IoTools::Formatter::printLine(std::cout, '-');
            IoTools::Formatter::printNewline(std::cout);
         }
      }
      IoTools::Formatter::printNewline(std::cout);
   }
}
}
