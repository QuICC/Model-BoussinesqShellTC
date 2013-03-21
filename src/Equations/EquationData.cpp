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
      : mspEqParams(spEqParams), mEqIsComplex(false), mZeroIdx(0)
   {
   }

   EquationData::~EquationData()
   {
   }

   bool EquationData::isComplex(FieldComponents::Spectral::Id id) const
   {
      return this->mEqIsComplex;
   }

   int EquationData::startIndex(FieldComponents::Spectral::Id id) const
   {
      return this->mZeroIdx;
   }

   int EquationData::rowShift(FieldComponents::Spectral::Id id, const int j) const
   {
      if(this->mCMatrices.count(id) == 0)
      {
         return 0;
      } else
      {
         return this->mCMatrices.find(id)->second.at(j).first.cols();
      }
   }

   void EquationData::finalizeMatrices()
   {
      // Matrix vector iterator for operators
      std::vector<DecoupledZSparse>::iterator opIt;
      // Matrix vector iterator for boundary conditions
      std::vector<DecoupledZSparse>::iterator bcIt;

      // Matrix map iterator
      std::map<FieldComponents::Spectral::Id, std::vector<DecoupledZSparse> >::iterator mapIt;

      // Loop over all components in time matrices map
      for(mapIt = this->mTMatrices.begin(); mapIt != this->mTMatrices.end(); ++mapIt)
      {
         // Make sure there are the right number of operators vs boundary matrices
         if(mapIt->second.size() != this->mBCMatrices.find(mapIt->first)->second.size())
         {
            throw Exception("Can't finalize matrices, incompatible operator and boundary matrices");
         }

         // Find the corresponding boundary condition matrix set
         bcIt = this->mBCMatrices.find(mapIt->first)->second.begin();
   
         // Loop over matrices components in map
         for(opIt = mapIt->second.begin(); opIt != mapIt->second.end(); ++opIt,++bcIt)
         {
            // If real matrix is empty resize it to right size
            if(opIt->first.size() == 0)
            {
               opIt->first.resize(opIt->second.rows(), opIt->second.cols());
            }
            opIt->first.makeCompressed();

            // If imaginary matrix is empty resize it to right size
            if(opIt->second.size() == 0)
            {
               opIt->second.resize(opIt->first.rows(), opIt->first.cols());
            }
            opIt->second.makeCompressed();

            // If real boundary matrix is empty resize it to right size
            if(bcIt->first.size() == 0)
            {
               bcIt->first.resize(opIt->first.rows(), opIt->first.cols());
            }
            bcIt->first.makeCompressed();

            // If imaginary boundary matrix is empty resize it to right size
            if(bcIt->second.size() == 0)
            {
               bcIt->second.resize(opIt->first.rows(), opIt->first.cols());
            }
            bcIt->second.makeCompressed();
         }
      }

      // Loop over all components in linear matrices map
      for(mapIt = this->mLMatrices.begin(); mapIt != this->mLMatrices.end(); ++mapIt)
      {
         // Loop over matrices components in map
         for(opIt = mapIt->second.begin(); opIt != mapIt->second.end(); ++opIt)
         {
            // If real matrix is empty resize it to right size
            if(opIt->first.size() == 0)
            {
               opIt->first.resize(opIt->second.rows(), opIt->second.cols());
            }
            opIt->first.makeCompressed();

            // If imaginary matrix is empty resize it to right size
            if(opIt->second.size() == 0)
            {
               opIt->second.resize(opIt->first.rows(), opIt->first.cols());
            }
            opIt->second.makeCompressed();
         }
      }

      // Loop over all components in coupling matrices map
      for(mapIt = this->mCMatrices.begin(); mapIt != this->mCMatrices.end(); ++mapIt)
      {
         // Make sure there are the right number of operators vs boundary matrices
         if(mapIt->second.size() != this->mCBCMatrices.find(mapIt->first)->second.size())
         {
            throw Exception("Can't finalize matrices, incompatible operator and boundary matrices");
         }

         // Find the corresponding boundary condition matrix set
         bcIt = this->mCBCMatrices.find(mapIt->first)->second.begin();
   
         // Loop over matrices components in map
         for(opIt = mapIt->second.begin(); opIt != mapIt->second.end(); ++opIt,++bcIt)
         {
            // If real matrix is empty resize it to right size
            if(opIt->first.size() == 0)
            {
               opIt->first.resize(opIt->second.rows(), opIt->second.cols());
            }
            opIt->first.makeCompressed();

            // If imaginary matrix is empty resize it to right size
            if(opIt->second.size() == 0)
            {
               opIt->second.resize(opIt->first.rows(), opIt->first.cols());
            }
            opIt->second.makeCompressed();

            // If real boundary matrix is empty resize it to right size
            if(bcIt->first.size() == 0)
            {
               bcIt->first.resize(opIt->first.rows(), opIt->first.cols());
            }
            bcIt->first.makeCompressed();

            // If imaginary boundary matrix is empty resize it to right size
            if(bcIt->second.size() == 0)
            {
               bcIt->second.resize(opIt->first.rows(), opIt->first.cols());
            }
            bcIt->second.makeCompressed();
         }
      }
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

   const DecoupledZSparse& EquationData::timeMatrix(FieldComponents::Spectral::Id id, const int j) const
   {
      // Safety assert
      assert(this->mTMatrices.count(id) > 0);
      
      return this->mTMatrices.find(id)->second.at(j);
   }

   const DecoupledZSparse& EquationData::linearMatrix(FieldComponents::Spectral::Id id, const int j) const
   {
      // Safety assert
      assert(this->mLMatrices.count(id) > 0);
      
      return this->mLMatrices.find(id)->second.at(j);
   }

   const DecoupledZSparse& EquationData::couplingMatrix(FieldComponents::Spectral::Id id, const int j) const
   {
      // Safety assert
      assert(this->mCMatrices.count(id) > 0);
      
      return this->mCMatrices.find(id)->second.at(j);
   }

   const DecoupledZSparse& EquationData::bcMatrix(FieldComponents::Spectral::Id id, const int j) const
   {
      // Safety assert
      assert(this->mBCMatrices.count(id) > 0);
      
      return this->mBCMatrices.find(id)->second.at(j);
   }

   const DecoupledZSparse& EquationData::cbcMatrix(FieldComponents::Spectral::Id id, const int j) const
   {
      // Safety assert
      assert(this->mCBCMatrices.count(id) > 0);
      
      return this->mCBCMatrices.find(id)->second.at(j);
   }

   const CouplingInformation& EquationData::couplingInfo() const
   {
      return this->mCouplingInfo;
   }

   void EquationData::setName(PhysicalNames::Id name)
   {
      this->mName = name;
   }

   void EquationData::setComplex(bool isComplex)
   {
      this->mEqIsComplex = isComplex;
   }

   void EquationData::setStartIndex(int start)
   {
      this->mZeroIdx = start;
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
