/** 
 * @file IVariableAsciiEWriter.cpp 
 * @brief Source of the implementation of the generic variable to ASCII file writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "IoVariable/IVariableAsciiEWriter.hpp"

// Project includes
//

namespace QuICC {

namespace IoVariable {

   IVariableAsciiEWriter::IVariableAsciiEWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id)
      : IAsciiEWriter(name, ext, header, type, version), mTime(-1.0), mTimestep(-1.0), mSpaceId(id)
   {
   }

   IVariableAsciiEWriter::~IVariableAsciiEWriter()
   {
   }

   Dimensions::Space::Id IVariableAsciiEWriter::space() const
   {
      return this->mSpaceId;
   }

   void IVariableAsciiEWriter::setPhysical(const std::map<std::string,MHDFloat>& parameters)
   {
      this->mPhysical = parameters;
   }

   void IVariableAsciiEWriter::setSimTime(const MHDFloat time, const MHDFloat timestep)
   {
      this->mTime = time;

      this->mTimestep = timestep;
   }

   void IVariableAsciiEWriter::setResolution(SharedResolution spRes)
   {
      // Store resolution object
      this->mspRes = spRes;
   }

   void IVariableAsciiEWriter::expect(const PhysicalNames::Id id)
   {
      this->mExpected.insert(id);
   }

   bool IVariableAsciiEWriter::isFull() const
   {
      bool status = true;

      // Check that all expected scalars and vectors are present
      status = status && (this->mScalars.size() + this->mVectors.size() == this->mExpected.size());

      // Check that the resolution has been set
      status = status && this->mspRes;

      return status;
   }

   void IVariableAsciiEWriter::addScalar(const std::pair<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalar)
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
 
   void IVariableAsciiEWriter::addVector(const std::pair<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vector)
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

   IVariableAsciiEWriter::scalar_iterator_range  IVariableAsciiEWriter::scalarRange()
   {
      return std::make_pair(this->mScalars.begin(), this->mScalars.end());
   }

   IVariableAsciiEWriter::vector_iterator_range  IVariableAsciiEWriter::vectorRange()
   {
      return std::make_pair(this->mVectors.begin(), this->mVectors.end());
   }

   void IVariableAsciiEWriter::compute(Transform::TransformCoordinatorType& coord)
   {
   }

}
}
