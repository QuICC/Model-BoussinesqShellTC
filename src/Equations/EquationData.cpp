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

   EquationData::EquationData(SharedIEquationParameters spEqParams)
      : mspEqParams(spEqParams), mBCCounter(-1), mEqIsComplex(false)
   {
   }

   EquationData::~EquationData()
   {
   }

   bool EquationData::isInitialized() const
   {
      bool status;

      // Check boundary conditions initialisation
      status = this->initializedBCs();

      return status;
   }

   bool EquationData::isComplex(FieldComponents::Spectral::Id id) const
   {
      return (this->mEqIsComplex || this->boundaryIsComplex(id));
   }

   bool EquationData::boundaryIsComplex(FieldComponents::Spectral::Id id) const
   {
      bool status = false;

      return status;
   }

   int EquationData::nBC(FieldComponents::Spectral::Id id, Dimensions::Transform::Id dim) const
   {
      int n = 0;

      if(this->mBCs.count(std::make_pair(id,dim)) > 0)
      {
         std::map<std::pair<FieldComponents::Spectral::Id,Dimensions::Transform::Id>, std::map<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position> >::const_iterator bcIt = this->mBCs.find(std::make_pair(id,dim));
         std::map<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position>::const_iterator mapIt;
         for(mapIt = bcIt->second.begin(); mapIt != bcIt->second.end(); mapIt++)
         {
            n++;
         }
      }

      return n;
   }

   void EquationData::addBC(FieldComponents::Spectral::Id id, Dimensions::Transform::Id dim, const std::pair<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position>& bc, ArrayZ val)
   {
      // Create boundary condition and boundary value maps if required
      std::pair<FieldComponents::Spectral::Id,Dimensions::Transform::Id> myId = std::make_pair(id,dim);
      if(this->mBCs.count(myId) == 0)
      {
         // Create boundary condition map
         this->mBCs.insert(std::make_pair(myId, std::map<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position>()));

         // Create boundary value map
         this->mBVals.insert(std::make_pair(myId, std::map<Spectral::BoundaryConditions::Id,ArrayZ>()));
      }

      // Insert boundary condition
      this->mBCs.find(myId)->second.insert(bc);

      // Insert boundary value
      this->mBVals.find(myId)->second.insert(std::make_pair(bc.first,val));
   }

   void EquationData::addCBC(FieldComponents::Spectral::Id id, Dimensions::Transform::Id dim, const std::pair<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position>& bc)
   {
      std::pair<FieldComponents::Spectral::Id,Dimensions::Transform::Id> myId = std::make_pair(id,dim);
      if(this->mCBCs.count(myId) == 0)
      {
         this->mCBCs.insert(std::make_pair(myId, std::map<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position>()));
      }

      this->mCBCs.find(myId)->second.insert(bc);
   }

   bool EquationData::initializedBCs() const
   {
      // Check that all the boundary conditions have been set
      if(this->mBCCounter == 0)
      {
         return true;
      } else
      {
         return false;
      }
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

   const std::map<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position>& EquationData::getBCs(FieldComponents::Spectral::Id id, Dimensions::Transform::Id dim) const
   {
      // Safety assert
      assert(this->mBCs.count(std::make_pair(id,dim)) > 0);

      return this->mBCs.find(std::make_pair(id,dim))->second;
   }

   const std::map<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position>& EquationData::getCBCs(FieldComponents::Spectral::Id id, Dimensions::Transform::Id dim) const
   {
      // Safety assert
      assert(this->mCBCs.count(std::make_pair(id,dim)) > 0);
      
      return this->mCBCs.find(std::make_pair(id,dim))->second;
   }

   const std::map<Spectral::BoundaryConditions::Id,ArrayZ>& EquationData::getBVals(FieldComponents::Spectral::Id id, Dimensions::Transform::Id dim) const
   {
      // Safety assert
      assert(this->mBVals.count(std::make_pair(id,dim)) > 0);

      return this->mBVals.find(std::make_pair(id,dim))->second;
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
