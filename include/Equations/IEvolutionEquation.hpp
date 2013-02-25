/** \file IEvolutionEquation.hpp
 *  \brief Base building block for the implementation of a time dependend evolution equation
 */

#ifndef IEVOLUTIONEQUATION_HPP
#define IEVOLUTIONEQUATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"
#include "Simulation/PrepMacros/VariableTypedefsMacro.h"
#include "Simulation/PrepMacros/SpectralOperatorTypedefsMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "Enums/PhysicalNames.hpp"
#include "Enums/FieldComponents.h"
#include "Simulation/Enums/BoundaryConditions.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Base building block for the implementation of a time dependend evolution equation
    */
   class IEvolutionEquation
   {
      public:
         /**
          * \brief Simple constructor
          */
         IEvolutionEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IEvolutionEquation();

         /**
          * @brief Get name ID of the unknown
          */
         PhysicalNames::Id name() const;

         /**
          * @brief Set the smart pointer to the scalar field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the scalar field
          */
         void setField(PhysicalNames::Id name, Code::SharedScalarVariable spField);

         /**
          * @brief Set the smart pointer to the vector field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         void setField(PhysicalNames::Id name, Code::SharedVectorVariable spField);

         /**
          * @brief Initialise the equation
          */
         virtual void init();
   
         /**
          * @brief Get the number of implemented boundary conditions
          */
         int nBC(FieldComponents::Spectral::Component id, Dimensions::Type dim) const;

         /**
          * @brief Add boundary condition
          *
          * @param id   Spectral component ID
          * @param dim  ID of the dimension
          * @param bc   Boundary condition 
          * @param val  Boundary value (default to zero)
          */
         void addBC(FieldComponents::Spectral::Component id, Dimensions::Type dim, const std::pair<BoundaryConditions::Id,BoundaryConditions::Position>& bc, ArrayZ val = ArrayZ());

         /**
          * @brief Get boundary conditions
          */
         const std::map<BoundaryConditions::Id,BoundaryConditions::Position>& getBCs(FieldComponents::Spectral::Component id, Dimensions::Type dim) const;

         /**
          * @brief Add boundary condition
          */
         void addCBC(FieldComponents::Spectral::Component id, Dimensions::Type dim, const std::pair<BoundaryConditions::Id,BoundaryConditions::Position>& bc);

         /**
          * @brief Get coupled boundary conditions
          */
         const std::map<BoundaryConditions::Id,BoundaryConditions::Position>& getCBCs(FieldComponents::Spectral::Component id, Dimensions::Type dim) const;

         /**
          * @brief Get boundary values
          */
         const std::map<BoundaryConditions::Id,ArrayZ>& getBVals(FieldComponents::Spectral::Component id, Dimensions::Type dim) const;

         /**
          * @brief Check if equation is initialised correctly
          */
         virtual bool isInitialised() const;

         /**
          * @brief Are equation timestepping matrices complex?
          */
         virtual bool isComplex(FieldComponents::Spectral::Component id) const;

         /**
          * @brief Get the Time derivative matrices (sparse matrices)
          */
         const DecoupledZSparse& timeMatrix(FieldComponents::Spectral::Component id, const int j) const;

         /**
          * @brief Get the linear operator matrices (sparse matrices) 
          */
         const DecoupledZSparse& linearMatrix(FieldComponents::Spectral::Component id, const int j) const;

         /**
          * @brief Get the coupling matrices (sparse matrices) 
          */
         const DecoupledZSparse& couplingMatrix(FieldComponents::Spectral::Component id, const int j) const;

         /**
          * @brief Get the boundary condition matrices (sparse matrices) 
          */
         const DecoupledZSparse& bcMatrix(FieldComponents::Spectral::Component id, const int j) const;

         /**
          * @brief Get the coupling boundary condition matrices (sparse matrices) 
          */
         const DecoupledZSparse& cbcMatrix(FieldComponents::Spectral::Component id, const int j) const;

         /**
          * @brief Get the coupling information
          */
         const CouplingInformation&  couplingInfo() const;

         /**
          * @brief Get the row shift due to coupling matrices
          */
         int rowShift(FieldComponents::Spectral::Component id, const int j) const;

         /**
          * @brief Pure virtual method to transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Component id, DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Pure virtual method to transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Component id, MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Pure virtual method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Component id, const DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Pure virtual method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Component id, const MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Get map of field storage requirements information
          */
         const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& requirements() const;

         /**
          * @brief Get the field storage requirements information
          *
          * The TriBool vector contains the following information
          *    TriBool(0): Requires spectral coefficients?
          *    TriBool(1): Requires physical space values?
          *    TriBool(2): Requires physical space gradient/curl values?
          *
          * @param name Physical name of the field
          */
         TriBool requirements(PhysicalNames::Id name) const;

         /**
          * @brief Set the equation matrices
          */
         virtual void setSpectralMatrices(Code::SpectralOperator1DType& spec1D, Code::SpectralOperator2DType& spec2D, Code::SpectralOperator3DType& spec3D) = 0;

         /**
          * @brief Finalize the initialise equation matrices
          */
         void finalizeMatrices();

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Component id, DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Component id, MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Component id, const DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Component id, const MatrixZ& storage, const int matIdx, const int start) = 0;
         
      protected:
         /**
          * @brief Apply quasi-inverse to nonlinear values
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void applyQuasiInverse(FieldComponents::Spectral::Component id, DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Apply quasi-inverse to nonlinear values
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void applyQuasiInverse(FieldComponents::Spectral::Component id, MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Check if boundary conditions are complex
          */
         bool boundaryIsComplex(FieldComponents::Spectral::Component id) const;

         /**
          * @brief Set the equation variable requirements
          */
         virtual void setRequirements() = 0;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() = 0;

         /**
          * @brief Get the equation parameters
          */
         const EquationParameters& eqParams() const;

         /**
          * @brief Check if all BCs are set
          */
         bool initialisedBCs() const;

         /**
          * @brief Get scalar variable
          *
          * @param name Physical name of the field
          */
         const Code::ScalarVariable& scalar(PhysicalNames::Id name) const;

         /**
          * @brief Set scalar variable
          *
          * @param name Physical name of the field
          */
         Code::ScalarVariable& rScalar(PhysicalNames::Id name);

         /**
          * @brief Get vector variable
          *
          * @param name Physical name of the field
          */
         const Code::VectorVariable& vector(PhysicalNames::Id name) const;

         /**
          * @brief Set vector variable
          *
          * @param name Physical name of the field
          */
         Code::VectorVariable& rVector(PhysicalNames::Id name);

         /**
          * @brief Storage for smart equation parameters
          */
         SharedEquationParameters   mspEqParams;

         /**
          * @brief Boundary condition counter
          */
         int mBCCounter;

         /**
          * @brief Map of name and pointer for the scalar variables
          */
         std::map<PhysicalNames::Id, Code::SharedScalarVariable>  mScalars;

         /**
          * @brief Map of name and pointer for the vector variables
          */
         std::map<PhysicalNames::Id, Code::SharedVectorVariable>  mVectors;

         /**
          * @brief Map of component and time matrices 
          */
         std::map<FieldComponents::Spectral::Component, std::vector<DecoupledZSparse> > mTMatrices;

         /**
          * @brief Map of component and linear matrices
          */
         std::map<FieldComponents::Spectral::Component, std::vector<DecoupledZSparse> > mLMatrices;

         /**
          * @brief Map of component and coupling matrices
          */
         std::map<FieldComponents::Spectral::Component, std::vector<DecoupledZSparse> > mCMatrices;

         /**
          * @brief Map of component and boundary condition matrices
          */
         std::map<FieldComponents::Spectral::Component, std::vector<DecoupledZSparse> > mBCMatrices;

         /**
          * @brief Map of component and coupled boundary condition matrices
          */
         std::map<FieldComponents::Spectral::Component, std::vector<DecoupledZSparse> > mCBCMatrices;

         /**
          * @brief Map of component and nonlinear term multiplication matrices
          */
         std::map<FieldComponents::Spectral::Component, std::vector<SparseMatrix> > mNLMatrices;

         /**
          * @brief Storage for the boundary condition
          */
         std::map<std::pair<FieldComponents::Spectral::Component,Dimensions::Type>, std::map<BoundaryConditions::Id,BoundaryConditions::Position> > mBCs;

         /**
          * @brief Storage for the coupled boundary condition
          */
         std::map<std::pair<FieldComponents::Spectral::Component,Dimensions::Type>, std::map<BoundaryConditions::Id,BoundaryConditions::Position> > mCBCs;

         /**
          * @brief Storage for the boundary value
          */
         std::map<std::pair<FieldComponents::Spectral::Component,Dimensions::Type>, std::map<BoundaryConditions::Id,ArrayZ> > mBVals;

         /**
          * @brief Coupling information of the equation
          */
         CouplingInformation  mCouplingInfo;

         /**
          * @brief Name ID of the unknown
          */
         PhysicalNames::Id mName;

         /**
          * @brief Storage for the variable requirements
          */
         std::map<PhysicalNames::Id, std::pair<bool, TriBool> >   mRequirements;

         /**
          * @brief Are timestepping matrices for this equation complex?
          */
         bool mEqIsComplex;

      private:
   };

   inline void IEvolutionEquation::setField(PhysicalNames::Id name, Code::SharedScalarVariable spField)
   {
      this->mScalars.insert(std::make_pair(name, spField));
   }

   inline void IEvolutionEquation::setField(PhysicalNames::Id name, Code::SharedVectorVariable spField)
   {
      this->mVectors.insert(std::make_pair(name, spField));
   }

   inline const Code::ScalarVariable& IEvolutionEquation::scalar(PhysicalNames::Id name) const
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 1);

      return *(this->mScalars.find(name)->second);
   }

   inline Code::ScalarVariable& IEvolutionEquation::rScalar(PhysicalNames::Id name)
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 1);

      return *(this->mScalars.find(name)->second);
   }

   inline const Code::VectorVariable& IEvolutionEquation::vector(PhysicalNames::Id name) const
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 1);

      return *(this->mVectors.find(name)->second);
   }

   inline Code::VectorVariable& IEvolutionEquation::rVector(PhysicalNames::Id name)
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 1);

      return *(this->mVectors.find(name)->second);
   }

   inline const EquationParameters& IEvolutionEquation::eqParams() const
   {
      return *this->mspEqParams;
   }

   inline const std::map<BoundaryConditions::Id,BoundaryConditions::Position>& IEvolutionEquation::getBCs(FieldComponents::Spectral::Component id, Dimensions::Type dim) const
   {
      // Safety assert
      assert(this->mBCs.count(std::make_pair(id,dim)) > 0);

      return this->mBCs.find(std::make_pair(id,dim))->second;
   }

   inline const std::map<BoundaryConditions::Id,BoundaryConditions::Position>& IEvolutionEquation::getCBCs(FieldComponents::Spectral::Component id, Dimensions::Type dim) const
   {
      // Safety assert
      assert(this->mCBCs.count(std::make_pair(id,dim)) > 0);
      
      return this->mCBCs.find(std::make_pair(id,dim))->second;
   }

   inline const std::map<BoundaryConditions::Id,ArrayZ>& IEvolutionEquation::getBVals(FieldComponents::Spectral::Component id, Dimensions::Type dim) const
   {
      // Safety assert
      assert(this->mBVals.count(std::make_pair(id,dim)) > 0);

      return this->mBVals.find(std::make_pair(id,dim))->second;
   }

   inline const DecoupledZSparse& IEvolutionEquation::timeMatrix(FieldComponents::Spectral::Component id, const int j) const
   {
      // Safety assert
      assert(this->mTMatrices.count(id) > 0);
      
      return this->mTMatrices.find(id)->second.at(j);
   }

   inline const DecoupledZSparse& IEvolutionEquation::linearMatrix(FieldComponents::Spectral::Component id, const int j) const
   {
      // Safety assert
      assert(this->mLMatrices.count(id) > 0);
      
      return this->mLMatrices.find(id)->second.at(j);
   }

   inline const DecoupledZSparse& IEvolutionEquation::couplingMatrix(FieldComponents::Spectral::Component id, const int j) const
   {
      // Safety assert
      assert(this->mCMatrices.count(id) > 0);
      
      return this->mCMatrices.find(id)->second.at(j);
   }

   inline const DecoupledZSparse& IEvolutionEquation::bcMatrix(FieldComponents::Spectral::Component id, const int j) const
   {
      // Safety assert
      assert(this->mBCMatrices.count(id) > 0);
      
      return this->mBCMatrices.find(id)->second.at(j);
   }

   inline const DecoupledZSparse& IEvolutionEquation::cbcMatrix(FieldComponents::Spectral::Component id, const int j) const
   {
      // Safety assert
      assert(this->mCBCMatrices.count(id) > 0);
      
      return this->mCBCMatrices.find(id)->second.at(j);
   }

   inline const CouplingInformation& IEvolutionEquation::couplingInfo() const
   {
      return this->mCouplingInfo;
   }

   inline PhysicalNames::Id   IEvolutionEquation::name() const
   {
      return this->mName;
   }

   inline const std::map<PhysicalNames::Id, std::pair<bool, TriBool> >&  IEvolutionEquation::requirements() const
   {
      return this->mRequirements;
   }

   inline TriBool IEvolutionEquation::requirements(PhysicalNames::Id name) const
   {
      if(this->mRequirements.count(name) > 0)
      {
         return this->mRequirements.find(name)->second.second;
      } else
      {
         return TriBool(false,false,false);
      }
   }

   /// Typedef for a smart IEvolutionEquation
   typedef SharedPtrMacro<IEvolutionEquation> SharedIEvolutionEquation;
}

#endif // IEVOLUTIONEQUATION_HPP
