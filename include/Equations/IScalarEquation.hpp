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
#include "TypeSelectors/EquationEigenSelector.hpp"
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
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents();

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
    * @param isSet      Arithmetic operation is set
    */
   template <typename TData> void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet);

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

   /**
    * @brief Set nonlinear spectral values to zero
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void setZeroNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

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

         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            // Add source data
            int l;
            int j_;
            int dimI = this->spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            #if defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
               int corrDim = this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
            for(int j = 0; j < cols; j++)
            {
               j_ = this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
               #if defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
                  j_ -= corrDim;
               #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
               for(int i = 0; i < rows; i++)
               {
                  // Compute correct position
                  l = start + j_ + i;

                  // Copy timestep output into field
                  this->rUnknown().rDom(0).rPerturbation().setPoint(Datatypes::internal::getScalar(*solution, l),i,j,matIdx);
               }
            }
         #else
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
         #endif //defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         int rows = this->unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().slice(matIdx).cols();

         // Copy data
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy timestep output into field
               this->rUnknown().rDom(0).rPerturbation().setPoint(Datatypes::internal::getScalar(*solution, i + solStart, j),i,j,matIdx);
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
      TData rhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      Equations::copyUnknown(eq, compId, rhs, matIdx, 0, true, true);

      // Get a restricted stencil matrix
      SparseMatrix stencil(eq.couplingInfo(compId).galerkinN(matIdx),eq.couplingInfo(compId).galerkinN(matIdx));
      eq.dispatchGalerkinStencil(compId, stencil, matIdx, eq.unknown().dom(0).spRes(), EigenSelector::getEigs(eq, matIdx), true);
      stencil.makeCompressed();

      // Create solver and factorize stencil
      Solver::SparseSelector<SparseMatrix>::Type solver;
      solver.compute(stencil);
      // Safety assert for successful factorisation
      if(solver.info() != Eigen::Success)
      {
         throw Exception("Stencil factorization for initial solution failed!");
      }

      // solve for galerkin expansion
      TData lhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      Solver::internal::solveWrapper(lhs, solver, rhs);
      internal::setTopBlock(storage, start, eq.couplingInfo(compId).galerkinN(matIdx), lhs);
   }

   template <typename TData> void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet)
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

         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            // Add source data
            int l;
            int j_;
            int dimI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            #if defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
               int corrDim = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
            if(isSet)
            {
               ///\mhdBug This is overkill
               // Set storage to zero
               setZeroNonlinear(eq, compId, storage, matIdx, start);

               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  #if defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
                     j_ -= corrDim;
                  #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::setScalar(storage, l, eq.unknown().dom(0).perturbation().point(i,j,matIdx));
                  }
               }
            } else
            {
               for(int j = zeroCol; j < eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(matIdx); j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::addScalar(storage, l, eq.unknown().dom(0).perturbation().point(i,j,matIdx));
                  }
               }
            }
         #else
            // Copy data
            int k = start;
            if(isSet)
            {
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
            } else
            {
               for(int j = zeroCol; j < cols; j++)
               {
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Copy field value into storage
                     Datatypes::internal::addScalar(storage, k, eq.unknown().dom(0).perturbation().point(i,j,matIdx));

                     // increase storage counter
                     k++;
                  }
               }
            }
         #endif //defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

      // matIdx is the index of the slowest varying direction with multiple RHS
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);

         // Copy data
         if(isSet)
         {
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Copy field value into storage
                  Datatypes::internal::setScalar(storage, i - zeroRow + start, j - zeroCol, eq.unknown().dom(0).perturbation().point(i,j,matIdx));
               }
            }
         } else
         {
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Copy field value into storage
                  Datatypes::internal::addScalar(storage, i - zeroRow + start, j - zeroCol, eq.unknown().dom(0).perturbation().point(i,j,matIdx));
               }
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
         if(isSet)
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Copy field value into storage
               Datatypes::internal::setScalar(storage, k, eq.unknown().dom(0).perturbation().point(i,mode(1),mode(0)));

               // increase storage counter
               k++;
            }
         } else
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Copy field value into storage
               Datatypes::internal::addScalar(storage, k, eq.unknown().dom(0).perturbation().point(i,mode(1),mode(0)));

               // increase storage counter
               k++;
            }
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
         if(isSet)
         {
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
         } else
         {
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
                     Datatypes::internal::addScalar(storage, l, eq.unknown().dom(0).perturbation().point(i,j,k));
                  }
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

      bool isInitialized = (eq.explicitTiming(compId) == ExplicitTiming::LINEAR && eq.couplingInfo(compId).nExplicit() > 0);
      bool copyIsSet = true;
      bool applyIsSet = true;

      // Check if a nonlinear computation took place and a quasi-inverse has to be applied
      if(eq.couplingInfo(compId).hasNonlinear() && eq.couplingInfo(compId).hasQuasiInverse())
      {
         TData * rhs;
         bool useShift;
         int copyStart;
         TData tmp;
         if(eq.couplingInfo(compId).isGalerkin() || isInitialized)
         {
            // Temporary storage is required
            tmp = TData(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
            rhs = &tmp;
            useShift = false;
            copyStart = 0;
            if(isInitialized)
            {
               applyIsSet = false;
            }
         } else
         {
            // Storage can be used a RHS and LHS
            rhs = &storage;
            useShift = true;
            copyStart = start;
         }

         // simply copy values from unknown
         copyUnknown(eq, compId, *rhs, matIdx, copyStart, useShift, copyIsSet);

         // Multiply nonlinear term by quasi-inverse
         applyQuasiInverse(eq, compId, storage, start, matIdx, copyStart, *rhs, applyIsSet);

      /// Nonlinear computation took place but no quas-inverse is required
      } else if(eq.couplingInfo(compId).hasNonlinear())
      {
         if(isInitialized)
         {
            copyIsSet = false;
         }

         // simply copy values from unknown
         copyUnknown(eq, compId, storage, matIdx, start, true, copyIsSet);

      // Without nonlinear computation the values have to be initialised to zero   
      } else if(! isInitialized)
      {
         setZeroNonlinear(eq, compId, storage, matIdx, start);
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

            #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
               // Add source data
               int l;
               int j_;
               int dimI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               #if defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
                  int corrDim = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
               #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  #if defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
                     j_ -= corrDim;
                  #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_BLFM
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Add source value
                     Datatypes::internal::addScalar(storage, l, eq.sourceTerm(compId, i, j, matIdx));
                  }
               }
            #else
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
            #endif //defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

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
                  Datatypes::internal::addScalar(storage, i - zeroRow + start, j - zeroCol, eq.sourceTerm(compId, i, j, matIdx));
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

   template <typename TData> void setZeroNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
      {

         //Safety assertion
         assert(start >= 0);

         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            for(int k = start; k < eq.couplingInfo(compId).galerkinN(matIdx); ++k)
            {
               // Set value to zero
               Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));
            }
         #else
            int rows = eq.unknown().dom(0).perturbation().slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);
            int zeroCol = eq.couplingInfo(compId).galerkinShift(1);

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
         #endif //defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

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
               Datatypes::internal::setScalar(storage, i - zeroRow + start, j - zeroCol, typename TData::Scalar(0.0));
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
}

#endif // ISCALAREQUATION_HPP
