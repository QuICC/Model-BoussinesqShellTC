/** \file IEvolutionEquation.cpp
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
#include "Equations/IEvolutionEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   IEvolutionEquation::IEvolutionEquation(SharedEquationParameters spEqParams)
      : mspEqParams(spEqParams), mBCCounter(-1), mEqIsComplex(false)
   {
   }

   bool IEvolutionEquation::isComplex(FieldComponents::Spectral::Component id) const
   {
      return (this->mEqIsComplex || this->boundaryIsComplex(id));
   }

   bool IEvolutionEquation::boundaryIsComplex(FieldComponents::Spectral::Component id) const
   {
      bool status = false;

      return status;
   }

   int IEvolutionEquation::nBC(FieldComponents::Spectral::Component id, Dimensions::Type dim) const
   {
      int n = 0;

      if(this->mBCs.count(std::make_pair(id,dim)) > 0)
      {
         std::map<std::pair<FieldComponents::Spectral::Component,Dimensions::Type>, std::map<BoundaryConditions::Id,BoundaryConditions::Position> >::const_iterator bcIt = this->mBCs.find(std::make_pair(id,dim));
         std::map<BoundaryConditions::Id,BoundaryConditions::Position>::const_iterator mapIt;
         for(mapIt = bcIt->second.begin(); mapIt != bcIt->second.end(); mapIt++)
         {
            if(mapIt->second == BoundaryConditions::BOTH)
            {
               n += 2;
            } else
            {
               n++;
            }
         }
      }

      return n;
   }

   void IEvolutionEquation::addBC(FieldComponents::Spectral::Component id, Dimensions::Type dim, const std::pair<BoundaryConditions::Id,BoundaryConditions::Position>& bc, ArrayZ val)
   {
      // Create boundary condition and boundary value maps if required
      std::pair<FieldComponents::Spectral::Component,Dimensions::Type> myId = std::make_pair(id,dim);
      if(this->mBCs.count(myId) == 0)
      {
         // Create boundary condition map
         this->mBCs.insert(std::make_pair(myId, std::map<BoundaryConditions::Id,BoundaryConditions::Position>()));

         // Create boundary value map
         this->mBVals.insert(std::make_pair(myId, std::map<BoundaryConditions::Id,ArrayZ>()));
      }

      // Insert boundary condition
      this->mBCs.find(myId)->second.insert(bc);

      // Insert boundary value
      this->mBVals.find(myId)->second.insert(std::make_pair(bc.first,val));
   }

   void IEvolutionEquation::addCBC(FieldComponents::Spectral::Component id, Dimensions::Type dim, const std::pair<BoundaryConditions::Id,BoundaryConditions::Position>& bc)
   {
      std::pair<FieldComponents::Spectral::Component,Dimensions::Type> myId = std::make_pair(id,dim);
      if(this->mCBCs.count(myId) == 0)
      {
         this->mCBCs.insert(std::make_pair(myId, std::map<BoundaryConditions::Id,BoundaryConditions::Position>()));
      }

      this->mCBCs.find(myId)->second.insert(bc);
   }

   void IEvolutionEquation::init()
   {
      this->setCoupling();
   }

   bool IEvolutionEquation::isInitialised() const
   {
      bool status;

      // Check boundary conditions initialisation
      status = this->initialisedBCs();

      return status;
   }

   bool IEvolutionEquation::initialisedBCs() const
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

   int IEvolutionEquation::rowShift(FieldComponents::Spectral::Component id, const int j) const
   {
      if(this->mCMatrices.count(id) == 0)
      {
         return 0;
      } else
      {
         return this->mCMatrices.find(id)->second.at(j).first.cols();
      }
   }

   void IEvolutionEquation::finalizeMatrices()
   {
      // Matrix vector iterator for operators
      std::vector<DecoupledZSparse>::iterator opIt;
      // Matrix vector iterator for boundary conditions
      std::vector<DecoupledZSparse>::iterator bcIt;

      // Matrix map iterator
      std::map<FieldComponents::Spectral::Component, std::vector<DecoupledZSparse> >::iterator mapIt;

      // Loop over all components in time matrices map
      for(mapIt = this->mTMatrices.begin(); mapIt != this->mTMatrices.end(); ++mapIt)
      {
         // Make sure there are the right number of operators vs boundary matrices
         if(mapIt->second.size() != this->mBCMatrices.find(mapIt->first)->second.size())
         {
            throw Exception("IEvolutionEquation::finalizeMatrices", "Incompatible opertor and boundary matrices vectors");
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
            throw Exception("IEvolutionEquation::finalizeMatrices", "Incompatible opertor and boundary matrices vectors");
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
}
