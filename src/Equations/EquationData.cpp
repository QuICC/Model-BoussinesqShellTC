/** 
 * @file EquationData.cpp
 * @brief Source of building block for the implementation of a time dependend evolution equation
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
#include "Equations/EquationData.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

   EquationData::EquationData(SharedEquationParameters spEqParams)
      : mspEqParams(spEqParams), mTime(-1.0)
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

   const SimulationBoundary& EquationData::bcIds() const
   {
      return *this->mspBcIds;
   }

   const SparseMatrix& EquationData::galerkinStencil(const FieldComponents::Spectral::Id compId, const int j) const
   {
      // Safety assert
      assert(this->mGStencils.count(compId) > 0);
      
      return this->mGStencils.find(compId)->second.at(j);
   }

   bool EquationData::hasQID(const FieldComponents::Spectral::Id compId) const
   {
      return (this->mQIDMatrices.count(compId) > 0);
   }

   bool EquationData::hasQIZ(const FieldComponents::Spectral::Id compId) const
   {
      return (this->mQIZMatrices.count(compId) > 0);
   }

   bool EquationData::hasExplicitDTerm(const ModelOperator::Id opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         return (this->mELDMatrices.count(key) > 0);

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR) 
      {
         return (this->mENLDMatrices.count(key) > 0);

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP) 
      {
         return (this->mENSDMatrices.count(key) > 0);

      } else
      {
         throw Exception("Requested inexistant explicit matrices");
      }
   }

   bool EquationData::hasExplicitZTerm(const ModelOperator::Id opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         return (this->mELZMatrices.count(key) > 0);

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR) 
      {
         return (this->mENLZMatrices.count(key) > 0);

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP) 
      {
         return (this->mENSZMatrices.count(key) > 0);

      } else
      {
         throw Exception("Requested inexistant explicit matrices");
      }
   }

   template <> const SparseMatrix& EquationData::quasiInverse<SparseMatrix>(const FieldComponents::Spectral::Id compId, const int j) const
   {
      // Safety assert
      assert(this->mQIDMatrices.count(compId) > 0);
      
      return this->mQIDMatrices.find(compId)->second.at(j);
   }

   template <> const SparseMatrixZ& EquationData::quasiInverse<SparseMatrixZ>(const FieldComponents::Spectral::Id compId, const int j) const
   {
      // Safety assert
      assert(this->mQIZMatrices.count(compId) > 0);
      
      return this->mQIZMatrices.find(compId)->second.at(j);
   }

   template <> const SparseMatrix& EquationData::explicitOperator<SparseMatrix>(const ModelOperator::Id opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         // Safety assert
         assert(this->mELDMatrices.count(key) > 0);
      
         return this->mELDMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR)
      {
         // Safety assert
         assert(this->mENLDMatrices.count(key) > 0);
      
         return this->mENLDMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP)
      {
         // Safety assert
         assert(this->mENSDMatrices.count(key) > 0);
      
         return this->mENSDMatrices.find(key)->second.at(j);

      } else
      {
         throw Exception("Requested inexistant explicit matrices");
      }
   }

   template <> const SparseMatrixZ& EquationData::explicitOperator<SparseMatrixZ>(const ModelOperator::Id opId, const FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const int j) const
   {
      // Make key
      std::pair<FieldComponents::Spectral::Id, SpectralFieldId> key = std::make_pair(compId, fieldId);

      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         // Safety assert
         assert(this->mELZMatrices.count(key) > 0);
      
         return this->mELZMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR)
      {
         // Safety assert
         assert(this->mENLZMatrices.count(key) > 0);
      
         return this->mENLZMatrices.find(key)->second.at(j);

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP)
      {
         // Safety assert
         assert(this->mENSZMatrices.count(key) > 0);
      
         return this->mENSZMatrices.find(key)->second.at(j);

      } else
      {
         throw Exception("Requested inexistant explicit matrices");
      }
   }

   std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrix> >& EquationData::rEDMatrices(const ModelOperator::Id opId)
   {
      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         return this->mELDMatrices;

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR)
      {
         return this->mENLDMatrices;

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP)
      {
         return this->mENSDMatrices;

      } else
      {
         throw Exception("Requested inexistant explicit matrices");
      }
   }

   std::map<std::pair<FieldComponents::Spectral::Id, SpectralFieldId>, std::vector<SparseMatrixZ> >& EquationData::rEZMatrices(const ModelOperator::Id opId)
   {
      if(opId == ModelOperator::EXPLICIT_LINEAR)
      {
         return this->mELZMatrices;

      } else if(opId == ModelOperator::EXPLICIT_NONLINEAR)
      {
         return this->mENLZMatrices;

      } else if(opId == ModelOperator::EXPLICIT_NEXTSTEP)
      {
         return this->mENSZMatrices;

      } else
      {
         throw Exception("Requested inexistant explicit matrices");
      }
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

   const VariableRequirement&  EquationData::imposedRequirements() const
   {
      return this->mImposedRequirements;
   }

   const FieldRequirement&  EquationData::requirements(PhysicalNames::Id id) const
   {
      return this->mRequirements.field(id);
   }

   const FieldRequirement&  EquationData::imposedRequirements(PhysicalNames::Id id) const
   {
      return this->mImposedRequirements.field(id);
   }

   FieldRequirement&  EquationData::updateFieldRequirements(PhysicalNames::Id id)
   {
      return this->mRequirements.rField(id);
   }

   void EquationData::setSolverIndex(const FieldComponents::Spectral::Id compId, const int idx)
   {
      // Safety assert
      assert(this->mCouplingInfos.count(compId) > 0);
      
      this->mCouplingInfos.find(compId)->second.setSolverIndex(idx);
   }

   void EquationData::setSolveTiming(const SolveTiming::Id time)
   {
      this->mSolveTiming = time;
   }

   SolveTiming::Id  EquationData::solveTiming() const
   {
      return this->mSolveTiming;
   }

   MHDFloat  EquationData::time() const
   {
      return this->mTime;
   }

   void  EquationData::setTime(const MHDFloat time, const bool finished)
   {
      this->mTime = time;
   }

   const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& EquationData::nlComponents() const
   {
      return this->mNLComponents;
   }

   void EquationData::addNLComponent(const FieldComponents::Spectral::Id compId, const int flag)
   {
      this->mNLComponents.push_back(std::make_pair(compId,flag));
   }

   MHDFloat incrementTimeAverage(const MHDComplex avg, const MHDFloat newData, const MHDFloat time, const MHDFloat timestep)
   {
      throw Exception("Setup is wrong, should not have been called!");

      return std::numeric_limits<MHDFloat>::quiet_NaN();
   }

   MHDFloat noupdateTimeAverage(const MHDComplex avg, const MHDFloat newData)
   {
      throw Exception("Setup is wrong, should not have been called!");

      return std::numeric_limits<MHDFloat>::quiet_NaN();
   }
}
}
