/** \file EquationData.cpp
 *  \brief Source of building block for the implementation of a time dependend evolution equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/EquationData.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   EquationData::EquationData(SharedIEquationParameters spEqParams)
      : mspEqParams(spEqParams)
   {
   }

   EquationData::~EquationData()
   {
   }

   void EquationData::setField(PhysicalNames::Id name, Datatypes::SharedScalarVariableType spField)
   {
      this->mScalars.insert(std::make_pair(name, spField));
   }

   void EquationData::setField(PhysicalNames::Id name, Datatypes::SharedVectorVariableType spField)
   {
      this->mVectors.insert(std::make_pair(name, spField));
   }

   const Datatypes::ScalarVariableType& EquationData::scalar(PhysicalNames::Id name) const
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 1);

      return *(this->mScalars.find(name)->second);
   }

   Datatypes::ScalarVariableType& EquationData::rScalar(PhysicalNames::Id name)
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 1);

      return *(this->mScalars.find(name)->second);
   }

   const Datatypes::VectorVariableType& EquationData::vector(PhysicalNames::Id name) const
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 1);

      return *(this->mVectors.find(name)->second);
   }

   Datatypes::VectorVariableType& EquationData::rVector(PhysicalNames::Id name)
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 1);

      return *(this->mVectors.find(name)->second);
   }

   const IEquationParameters& EquationData::eqParams() const
   {
      return *this->mspEqParams;
   }

   const DecoupledZSparse& EquationData::explicitLinear(FieldComponents::Spectral::Id id, const int j) const
   {
      // Safety assert
      assert(this->mLMatrices.count(id) > 0);
      
      return this->mLMatrices.find(id)->second.at(j);
   }

   const CouplingInformation& EquationData::couplingInfo(FieldComponents::Spectral::Id comp) const
   {
      // Safety assert
      assert(this->mCouplingInfos.count(comp) > 0);
      
      return this->mCouplingInfos.find(comp)->second;
   }

   void EquationData::setName(PhysicalNames::Id name)
   {
      this->mName = name;
   }

   PhysicalNames::Id   EquationData::name() const
   {
      return this->mName;
   }

   const VariableRequirement&  EquationData::requirements() const
   {
      return this->mRequirements;
   }

   const FieldRequirement&  EquationData::requirements(PhysicalNames::Id id) const
   {
      return this->mRequirements.field(id);
   }
}
}
