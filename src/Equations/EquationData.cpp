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

   EquationData::EquationData(SharedEquationParameters spEqParams)
      : mspEqParams(spEqParams)
   {
   }

   EquationData::~EquationData()
   {
   }

   void EquationData::setField(PhysicalNames::Id name, Datatypes::SharedScalarVariableType spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 0);

      this->mScalars.insert(std::make_pair(name, spField));
   }

   void EquationData::setField(PhysicalNames::Id name, Datatypes::SharedVectorVariableType spField)
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 0);

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

   const EquationParameters& EquationData::eqParams() const
   {
      return *this->mspEqParams;
   }

   const SparseMatrix& EquationData::quasiInverse(const FieldComponents::Spectral::Id compId, const int j) const
   {
      // Safety assert
      assert(this->mNLMatrices.count(compId) > 0);
      
      return this->mNLMatrices.find(compId)->second.at(j);
   }

   bool EquationData::hasExplicitDLinear(const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      return (this->mLDMatrices.count(key) > 0);
   }

   bool EquationData::hasExplicitZLinear(const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      return (this->mLZMatrices.count(key) > 0);
   }

   const SparseMatrix& EquationData::explicitDLinear(const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      // Safety assert
      assert(this->mLDMatrices.count(key) > 0);
      
      return this->mLDMatrices.find(key)->second.at(j);
   }

   const SparseMatrixZ& EquationData::explicitZLinear(const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      // Safety assert
      assert(this->mLZMatrices.count(key) > 0);
      
      return this->mLZMatrices.find(key)->second.at(j);
   }

   const CouplingInformation& EquationData::couplingInfo(const FieldComponents::Spectral::Id compId) const
   {
      // Safety assert
      assert(this->mCouplingInfos.count(compId) > 0);
      
      return this->mCouplingInfos.find(compId)->second;
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
