/**
 * @file IScalarEquation.hpp
 * @brief Base for the implementation of a scalar equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ISCALAREQUATION_HPP
#define ISCALAREQUATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Arithmetics.hpp"
#include "Enums/ModelOperatorBoundary.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IEquation.hpp"
#include "Base/DecoupledComplexInternal.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base for the implementation of a scalar equation
    */
   class IScalarEquation: public IEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IScalarEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarEquation();

         /**
          * @brief Set the shared pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Datatypes::SharedScalarVariableType spUnknown);

         /**
          * @brief Get the unknown variable
          */
         const Datatypes::ScalarVariableType& unknown() const;

         /**
          * @brief Access the shared resolution
          */
         virtual SharedResolution spRes() const;

         /**
          * @brief Update unknown from dealised data
          *
          * @param rhs     Dealised input
          * @param compId  ID of the component (allows for a more general implementation)
          * @param arithId ID of the arithmetic operation
          */
         void updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId, Arithmetics::Id arithId);

         /**
          * @brief Transfer solver solution to equation unknown
          *
          * @param compId  Component ID
          * @param storage Solver solution
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         template <typename TData> void storeSolution(FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start);

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds);

         /**
          * @brief Generic model operator dispatcher to python scripts
          */
         virtual void buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType) const;
         
      protected:
         /**
          * @brief Set the galerkin stencil
          */
         virtual void setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const;

         /**
          * @brief Set the quasi inverse matrix operator
          */
         virtual void setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx) const;

         /**
          * @brief Build coupling information from Python scripts
          */
         void defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource, const bool allowExplicit = true);

         /**
          * @brief Set the unknown variable
          */
         Datatypes::ScalarVariableType& rUnknown();

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedScalarVariableType mspUnknown;
   };

   /// Typedef for a shared IScalarEquation
   typedef SharedPtrMacro<IScalarEquation> SharedIScalarEquation;

   /**
    * @brief Solve for galerkin unknown using the stencil
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   template <typename TData> void solveStencilUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Copy unknown spectral values to solver
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    * @param useShift   Use galerkin shifts
    */
   template <typename TData> void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift);

   /**
    * @brief Transfer nonlinear spectral values from unknown to solver
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Add source term
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   template <typename TData> void IScalarEquation::storeSolution(FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      const TData * solution;
      TData tmp;
      int solStart;
      if(this->couplingInfo(compId).isGalerkin())
      {
         // Temporary storage is required
         tmp = TData(this->couplingInfo(compId).tauN(matIdx), this->couplingInfo(compId).rhsCols(matIdx));

         // Apply Galerkin stencil
         applyGalerkinStencil(*this, compId, tmp, start, matIdx, storage);

         solStart = 0;
         solution = &tmp;

      } else
      {
         solStart = start;
         solution = &storage;
      }

      if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
      {
         int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

         // Copy data
         int k = solStart;
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy timestep output into field
               this->rUnknown().rDom(0).rPerturbation().setPoint(Datatypes::internal::getScalar(*solution, k),i,j,matIdx);

               // increase linear storage counter
               k++;
            }
         }

      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

         // Copy data
         for(int j = 0; j < cols; j++)
         {
            for(int i = solStart; i < rows; i++)
            {
               // Copy timestep output into field
               this->rUnknown().rDom(0).rPerturbation().setPoint(Datatypes::internal::getScalar(*solution, i,j),i,j,matIdx);
            }
         }

      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = solStart;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            this->rUnknown().rDom(0).rPerturbation().setPoint(Datatypes::internal::getScalar(*solution,k),i,mode(1),mode(0));

            // increase linear storage counter
            k++;
         }

      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
      {
         assert(matIdx == 0);

         // Copy data
         int l;
         int k_;
         int j_;
         int dimK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
         int dimJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         for(int k = 0; k < this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
         {
            k_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
            for(int j = 0; j < this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               j_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
               for(int i = 0; i < this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
               {
                  // Compute correct position
                  l = solStart + k_ + j_ + i;

                  // Copy timestep output into field
                  this->rUnknown().rDom(0).rPerturbation().setPoint(Datatypes::internal::getScalar(*solution, l),i,j,k);
               }
            }
         }
      }
   }

   template <typename TData> void solveStencilUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Create temporary storage for tau data
      TData tmp(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      Equations::copyUnknown(eq, compId, tmp, matIdx, 0, false);
      TData rhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      internal::setTopBlock(rhs, 0, eq.couplingInfo(compId).galerkinN(matIdx), tmp);

      // Get a restricted stencil matrix
      SparseMatrix stencil = eq.galerkinStencil(compId, matIdx).topRows(eq.couplingInfo(compId).galerkinN(matIdx));
      stencil.makeCompressed();

      // Create solver and factorize stencil
      Solver::SparseSelector<SparseMatrix>::Type solver;
      solver.compute(stencil);
      // Safety assert for successful factorisation
      assert(solver.info() == Eigen::Success);

      // solve for galerkin expansion
      TData lhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      Solver::internal::solveWrapper(lhs, solver, rhs);
      internal::setTopBlock(storage, start, eq.couplingInfo(compId).galerkinN(matIdx), lhs);
   }

   template <typename TData> void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      int zeroRow = 0;
      int zeroCol = 0;
      if(useShift)
      {
         zeroRow = eq.couplingInfo(compId).galerkinShift(0);
         zeroCol = eq.couplingInfo(compId).galerkinShift(1);
      }

      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
      {
         int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);

         // Copy data
         int k = start;
         for(int j = zeroCol; j < cols; j++)
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Copy field value into storage
               Datatypes::internal::setScalar(storage, k, eq.unknown().dom(0).perturbation().point(i,j,matIdx));

               // increase storage counter
               k++;
            }
         }

      // matIdx is the index of the slowest varying direction with multiple RHS
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);

         // Copy data
         for(int j = zeroCol; j < cols; j++)
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Copy field value into storage
               Datatypes::internal::setScalar(storage, i, j, eq.unknown().dom(0).perturbation().point(i,j,matIdx));
            }
         }

      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         //Safety assertion
         assert(start >= 0);

         // Get mode indexes
         ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

         // Copy data
         int k = start;
         for(int i = zeroRow; i < rows; i++)
         {
            // Copy field value into storage
            Datatypes::internal::setScalar(storage, k, eq.unknown().dom(0).perturbation().point(i,mode(1),mode(0)));

            // increase storage counter
            k++;
         }

      // There is a single matrix
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
      {
         assert(matIdx == 0);

         //Safety assertion
         assert(start >= 0);

         // Copy data
         int l;
         int k_;
         int j_;
         int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
         int dimJ = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
         for(int k = 0; k < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
         {
            k_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
            for(int j = 0; j < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
               for(int i = 0; i < eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
               {
                  // Compute correct position
                  l = start + k_ + j_ + i;

                  // Copy field value into storage
                  Datatypes::internal::setScalar(storage, l, eq.unknown().dom(0).perturbation().point(i,j,k));
               }
            }
         }
      }
   }

   template <typename TData> void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);
      assert((!eq.couplingInfo(compId).isGalerkin() || eq.couplingInfo(compId).indexType() != CouplingInformation::SINGLE) && "Current version does not support galerkin basis");

      // Check if a nonlinear computation took place and a quasi-inverse has to be applied
      if(eq.couplingInfo(compId).hasNonlinear() && eq.couplingInfo(compId).hasQuasiInverse())
      {
         TData * rhs;
         bool useShift;
         int copyStart;
         TData tmp;
         if(eq.couplingInfo(compId).isGalerkin())
         {
            // Temporary storage is required
            tmp = TData(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
            rhs = &tmp;
            useShift = false;
            copyStart = 0;
         } else
         {
            // Storage can be used a RHS and LHS
            rhs = &storage;
            useShift = true;
            copyStart = start;
         }

         // simply copy values from unknown
         copyUnknown(eq, compId, *rhs, matIdx, copyStart, useShift);

         // Multiply nonlinear term by quasi-inverse
         applyQuasiInverse(eq, compId, storage, start, matIdx, copyStart, *rhs);

      /// Nonlinear computation took place but no quas-inverse is required
      } else if(eq.couplingInfo(compId).hasNonlinear())
      {
         // simply copy values from unknown
         copyUnknown(eq, compId, storage, matIdx, start, true);

      // Without nonlinear computation the values have to be initialised to zero   
      } else
      {
         // matIdx is the index of the slowest varying direction with a single RHS
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);
            int zeroCol = eq.couplingInfo(compId).galerkinShift(1);

            //Safety assertion
            assert(start >= 0);

            // Set data to zero
            int k = start;
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Set value to zero
                  Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

                  // increase storage counter
                  k++;
               }
            }

         // matIdx is the index of the slowest varying direction with multiple RHS
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);
            int zeroCol = eq.couplingInfo(compId).galerkinShift(1);

            //Safety assertion
            assert(start >= 0);

            // Set data to zero
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Set value to zero
                  Datatypes::internal::setScalar(storage, i, j, typename TData::Scalar(0.0));
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);

            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
            int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();

            // Set data to zero
            int k = start;
            for(int i = zeroRow; i < rows; i++)
            {
               // Set value to zero
               Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

               // increase storage counter
               k++;
            }

         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
         {
            //Safety assertion
            assert(matIdx == 0);
            assert(start >= 0);

            // Set data to zero
            int l;
            int k_;
            int j_;
            int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
            int dimJ = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            for(int k = 0; k < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
                  for(int i = 0; i < eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                  {
                     // Compute correct position
                     l = start + k_ + j_ + i;

                     // Set value to zero
                     Datatypes::internal::setScalar(storage, l, typename TData::Scalar(0.0));
                  }
               }
            }
         }
      }
   }

   template <typename TData> void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         // matIdx is the index of the slowest varying direction with a single RHS
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);
            int zeroCol = eq.couplingInfo(compId).galerkinShift(1);

            //Safety assertion
            assert(start >= 0);

            // Copy data
            int k = start;
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Add source value
                  Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, j, matIdx));

                  // increase storage counter
                  k++;
               }
            }

         // matIdx is the index of the slowest varying direction with multiple RHS
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);
            int zeroCol = eq.couplingInfo(compId).galerkinShift(1);

            //Safety assertion
            assert(start >= 0);

            // Copy data
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Add source value
                  Datatypes::internal::addScalar(storage, i, j, eq.sourceTerm(compId, i, j, matIdx));
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
            int rows = eq.unknown().dom(0).perturbation().slice(mode(0)).rows();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);

            // Copy data
            int k = start;
            for(int i = zeroRow; i < rows; i++)
            {
               // Get source value
               Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, mode(1), mode(0)));

               // increase storage counter
               k++;
            }

         // There is a single matrix
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
         {
            assert(matIdx == 0);

            //int zeroRow = eq.couplingInfo(compId).galerkinShift(0);
            //int zeroCol = eq.couplingInfo(compId).galerkinShift(1);
            //int zeroBlock = eq.couplingInfo(compId).galerkinShift(2);

            //Safety assertion
            assert(start >= 0);

            // Add source data
            int l;
            int k_;
            int j_;
            int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
            int dimJ = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            for(int k = 0; k < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
                  for(int i = 0; i < eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                  {
                     // Compute correct position
                     l = start + k_ + j_ + i;

                     // Add source value
                     Datatypes::internal::addScalar(storage, l, eq.sourceTerm(compId, i, j, k));
                  }
               }
            }
         }
      }
   }
}
}

#endif // ISCALAREQUATION_HPP
