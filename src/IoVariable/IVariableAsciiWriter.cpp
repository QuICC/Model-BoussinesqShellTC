/** 
 * @file IVariableAsciiWriter.cpp 
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
#include "IoVariable/IVariableAsciiWriter.hpp"

// Project includes
//

namespace QuICC {

namespace IoVariable {

   IVariableAsciiWriter::IVariableAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const WriteMode mode)
      : IAsciiWriter(name, ext, header, type, version, mode), mTime(-1.0), mTimestep(-1.0), mSpaceId(id)
   {
   }

   IVariableAsciiWriter::~IVariableAsciiWriter()
   {
   }

   Dimensions::Space::Id IVariableAsciiWriter::space() const
   {
      return this->mSpaceId;
   }

   void IVariableAsciiWriter::setPhysical(const std::map<std::string,MHDFloat>& parameters)
   {
      this->mPhysical = parameters;
   }

   const Resolution& IVariableAsciiWriter::res() const
   {
      return *this->mspRes;
   }

   void IVariableAsciiWriter::setSimTime(const MHDFloat time, const MHDFloat timestep)
   {
      this->mTime = time;

      this->mTimestep = timestep;
   }

   void IVariableAsciiWriter::setResolution(SharedResolution spRes)
   {
      // Store resolution object
      this->mspRes = spRes;
   }

   void IVariableAsciiWriter::expect(const PhysicalNames::Id id)
   {
      this->mExpected.insert(id);
   }

   bool IVariableAsciiWriter::isFull() const
   {
      bool status = true;

      // Check that all expected scalars and vectors are present
      status = status && (this->mScalars.size() + this->mVectors.size() == this->mExpected.size());

      // Check that the resolution has been set
      status = status && this->mspRes;

      return status;
   }

   void IVariableAsciiWriter::addScalar(const std::pair<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalar)
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
 
   void IVariableAsciiWriter::addVector(const std::pair<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vector)
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

   IVariableAsciiWriter::scalar_iterator_range  IVariableAsciiWriter::scalarRange()
   {
      return std::make_pair(this->mScalars.cbegin(), this->mScalars.cend());
   }

   IVariableAsciiWriter::vector_iterator_range  IVariableAsciiWriter::vectorRange()
   {
      return std::make_pair(this->mVectors.cbegin(), this->mVectors.cend());
   }

   void IVariableAsciiWriter::compute(Transform::TransformCoordinatorType& coord)
   {
   }

}
}
