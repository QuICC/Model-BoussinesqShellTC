/**
 * @file IVectorEquation.hpp
 * @brief Base for the implementation of a vector equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IVECTOREQUATION_HPP
#define IVECTOREQUATION_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/Arithmetics.hpp"
#include "Enums/ModelOperatorBoundary.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IEquation.hpp"
#include "Base/DecoupledComplexInternal.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SparseSolvers/SparseLinearSolverTools.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base for the implementation of a vector equation
    */
   class IVectorEquation: public IEquation
   {
      public:
         /// Typedef for the the spectral field component ID iterator
         typedef std::vector<FieldComponents::Spectral::Id>::const_iterator   SpectralComponent_iterator;

         /// Typedef for the the spectral field component ID iterator range
         typedef std::pair<SpectralComponent_iterator,SpectralComponent_iterator>  SpectralComponent_range;

         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IVectorEquation(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorEquation();

         /**
          * @brief Set the smart pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         void setUnknown(Datatypes::SharedVectorVariableType spUnknown);
         
         /**
          * @brief Get the unknown variable
          */
         const Datatypes::VectorVariableType& unknown() const;

         /**
          * @brief Access the shared resolution
          */
         virtual SharedResolution spRes() const;

         /**
          * @brief Get the number of spectral components
          */
         int nSpectral() const; 

         /**
          * @brief Get vector spectral component range
          */
         SpectralComponent_range spectralRange() const;

         /**
          * @brief Update unknown from dealised data
          *
          * @param rhs     Dealised input
          * @param compId  ID of the vector component
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
         virtual void initSpectralMatrices();

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
          * @brief Set the explicit matrix operator
          */
         virtual void setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const ModelOperator::Id opId, const SpectralFieldId fieldId, const int matIdx) const;

         /**
          * @brief Build coupling information from Python scripts
          */
         void defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasSource, const bool hasBoundaryValue = false, const bool allowExplicit = true);

         /**
          * @brief Set the unknown variable
          */
         Datatypes::VectorVariableType& rUnknown();

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedVectorVariableType mspUnknown;
   };

   /// Typedef for a shared IVectorEquation
   typedef SharedPtrMacro<IVectorEquation> SharedIVectorEquation;

   /**
    * @brief Solve for galerkin unknown using the stencil
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   template <typename TData> void solveStencilUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Copy unknown spectral values to solver
    *
    * @param eq         Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    * @param useShift   Use galerkin shifts
    * @param isSet      Arithmetic operation is set
    */
   template <typename TData> void copyUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet);

   /**
    * @brief Transfer nonlinear spectral values from unknown to solver
    *
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void copyNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Add source term
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void addSource(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Impose boundary value
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void setBoundaryValue(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Set nonlinear spectral values to zero
    *
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void setZeroNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   template <typename TData> void IVectorEquation::storeSolution(FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
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
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            // Add source data
            int l;
            int j_;
            int dimI = this->spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
               int corrDim = this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
            for(int j = 0; j < cols; j++)
            {
               j_ = this->spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
               #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                  j_ -= corrDim;
               #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
               for(int i = 0; i < rows; i++)
               {
                  // Compute correct position
                  l = start + j_ + i;

                  // Copy timestep output into field
                  typename Datatypes::internal::GetScalar<TData>::Scalar dataPoint = Datatypes::internal::getScalar(*solution, l);
                  dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, matIdx);
                  this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(dataPoint,i,j,matIdx);
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
                  typename Datatypes::internal::GetScalar<TData>::Scalar dataPoint = Datatypes::internal::getScalar(*solution, k);
                  dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, matIdx);
                  this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(dataPoint,i,j,matIdx);

                  // increase linear storage counter
                  k++;
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

         // Copy data
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy timestep output into field
               typename Datatypes::internal::GetScalar<TData>::Scalar dataPoint = Datatypes::internal::getScalar(*solution, i + solStart,j);
               dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, matIdx);
               this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(dataPoint,i,j,matIdx);
            }
         }

      } else if(this->couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         // Get mode indexes
         ArrayI mode = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

         // Copy data
         int k = solStart;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            typename Datatypes::internal::GetScalar<TData>::Scalar dataPoint = Datatypes::internal::getScalar(*solution, k);
            dataPoint = this->updateStoredSolution(dataPoint, compId, i, mode(1), mode(0));
            this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(dataPoint,i,mode(1),mode(0));

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
         #ifdef QUICC_SPATIALDIMENSION_3D
            int dimK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
         #else
            int dimK = 1;
         #endif //QUICC_SPATIALDIMENSION_3D
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
                  typename Datatypes::internal::GetScalar<TData>::Scalar dataPoint = Datatypes::internal::getScalar(*solution, l);
                  dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, k);
                  this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(dataPoint,i,j,k);
               }
            }
         }
      }
   }

   template <typename TData> void solveStencilUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Create temporary storage for tau data
      TData tmp(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      Equations::copyUnknown(eq, compId, tmp, matIdx, 0, false, true);
      TData rhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      #if defined QUICC_SPATIALSCHEME_WFT
         Datatypes::internal::setTopBlock(rhs, 0, eq.couplingInfo(compId).galerkinN(matIdx), eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL), eq.couplingInfo(compId).galerkinShift(matIdx, 0), tmp);
      #else 
         Datatypes::internal::setTopBlock(rhs, 0, eq.couplingInfo(compId).galerkinN(matIdx), tmp);
      #endif //defined QUICC_SPATIALSCHEME_WFT

      // Get a restricted stencil matrix
      SparseMatrix stencil(eq.couplingInfo(compId).galerkinN(matIdx),eq.couplingInfo(compId).galerkinN(matIdx));
      eq.dispatchGalerkinStencil(compId, stencil, matIdx, eq.unknown().dom(0).spRes(), eq.couplingInfo(compId).eigenTools().getEigs(eq.spRes(), matIdx), true);
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
      Datatypes::internal::setTopBlock(storage, start, eq.couplingInfo(compId).galerkinN(matIdx), lhs);
   }

   template <typename TData> void copyUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet)
   {
      int zeroRow = 0;
      int zeroCol = 0;
      if(useShift)
      {
         zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
         #if defined QUICC_SPATIALSCHEME_WFT
            zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
         #else
            zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
         #endif //defined QUICC_SPATIALSCHEME_WFT
      }

      // matIdx is the index of the slowest varying direction with single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
      {
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            // Add source data
            int l;
            int j_;
            int dimI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
               int corrDim = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
            if(isSet)
            {
               ///\mhdBug This is overkill
               // Set storage to zero
               setZeroNonlinear(eq, compId, storage, matIdx, start);

               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                     j_ -= corrDim;
                  #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::setScalar(storage, l, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx));
                  }
               }
            } else
            {
               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                     j_ -= corrDim;
                  #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::addScalar(storage, l, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx));
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
                     Datatypes::internal::setScalar(storage, k, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx));

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
                     Datatypes::internal::addScalar(storage, k, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx));

                     // increase storage counter
                     k++;
                  }
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

      // matIdx is the index of the slowest varying direction with single RHS
      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

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
                  Datatypes::internal::setScalar(storage, i - zeroRow + start, j - zeroCol, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx));
               }
            }
         } else
         {
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Copy field value into storage
                  Datatypes::internal::addScalar(storage, i - zeroRow + start, j - zeroCol, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,matIdx));
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
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

         #ifdef QUICC_SPATIALSCHEME_TFF
         // Filter out complex conjugate modes to be safe
         if(!(mode(3) == 0 && mode(2) > eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2))
         {
         #endif //QUICC_SPATIALSCHEME_TFF

         // Copy data
         int k = start;
         if(isSet)
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Copy field value into storage
               Datatypes::internal::setScalar(storage, k, eq.unknown().dom(0).perturbation().comp(compId).point(i,mode(1),mode(0)));

               // increase storage counter
               k++;
            }
         } else
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Copy field value into storage
               Datatypes::internal::addScalar(storage, k, eq.unknown().dom(0).perturbation().comp(compId).point(i,mode(1),mode(0)));

               // increase storage counter
               k++;
            }
         }

         #ifdef QUICC_SPATIALSCHEME_TFF
         } else
         {
            setZeroNonlinear(eq, compId, storage, matIdx, start);
         }
         #endif //QUICC_SPATIALSCHEME_TFF

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
         #ifdef QUICC_SPATIALDIMENSION_3D
            int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
         #else
            int dimK = 1;
         #endif //QUICC_SPATIALDIMENSION_3D
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
                     Datatypes::internal::setScalar(storage, l, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,k));
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
                     Datatypes::internal::addScalar(storage, l, eq.unknown().dom(0).perturbation().comp(compId).point(i,j,k));
                  }
               }
            }
         }
      }
   }

   template <typename TData> void copyNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Check if a nonlinear computation took place and a quasi-inverse has to be applied
      if(eq.couplingInfo(compId).hasNonlinear() && eq.couplingInfo(compId).hasQuasiInverse())
      {
         // Temporary storage is required
         TData tmp;
         tmp = TData(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));

         // simply copy values from unknown
         copyUnknown(eq, compId, tmp, matIdx, 0, false, true);

         // Multiply nonlinear term by quasi-inverse
         applyQuasiInverse(eq, compId, storage, start, matIdx, 0, tmp);

      /// Nonlinear computation took place but no quas-inverse is required
      } else if(eq.couplingInfo(compId).hasNonlinear())
      {
         // simply copy values from unknown
         copyUnknown(eq, compId, storage, matIdx, start, true, false);
      }
   }

   template <typename TData> void addSource(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         // matIdx is the index of the slowest varying direction with a single RHS
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            #if defined QUICC_SPATIALSCHEME_WFT
               int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            #else
               int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            #endif //defined QUICC_SPATIALSCHEME_WFT

            //Safety assertion
            assert(start >= 0);

            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               // Add source data
               int l;
               int j_;
               int dimI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                  int corrDim = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
               #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                     j_ -= corrDim;
                  #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
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
                     // Add source term
                     Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, j, matIdx));

                     // increase storage counter
                     k++;
                  }
               }
            #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

         // matIdx is the index of the slowest varying direction with multiple RHS
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            #if defined QUICC_SPATIALSCHEME_WFT
               int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            #else
               int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            #endif //defined QUICC_SPATIALSCHEME_WFT

            //Safety assertion
            assert(start >= 0);

            // Copy data
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Add source term
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
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);

            // Copy data
            int k = start;
            for(int i = zeroRow; i < rows; i++)
            {
               // Add source term
               Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, mode(1), mode(0)));

               // increase storage counter
               k++;
            }

         // There is a single matrix
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
         {
            assert(matIdx == 0);

            //int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            //int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            //int zeroBlock = eq.couplingInfo(compId).galerkinShift(matIdx,2);

            //Safety assertion
            assert(start >= 0);

            // Add source data
            int l;
            int k_;
            int j_;
            #ifdef QUICC_SPATIALDIMENSION_3D
               int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
            #else
               int dimK = 1;
            #endif //QUICC_SPATIALDIMENSION_3D
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

   template <typename TData> void setBoundaryValue(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Set boundary value if required
      if(eq.couplingInfo(compId).hasBoundaryValue())
      {
         if(eq.couplingInfo(compId).isGalerkin())
         {
            throw std::logic_error("Galerkin expansion cannot have a nonzero boundary value!");
         }

         std::vector<Eigen::Triplet<typename TData::Scalar> > triplets;
         typename TData::Scalar tripVal;

         // matIdx is the index of the slowest varying direction with a single RHS
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

            //Safety assertion
            assert(start >= 0);

            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               // Add source data
               int l;
               int j_;
               int dimI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                  int corrDim = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
               #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
               for(int j = 0; j < cols; j++)
               {
                  j_ = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  #if defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                     j_ -= corrDim;
                  #endif //defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
                  for(int i = 0; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Set boundary value
                     tripVal = eq.boundaryValue(compId, i, j, matIdx);
                     if(tripVal != 0.0)
                     {
                        triplets.push_back(Eigen::Triplet<typename TData::Scalar>(l, 0, tripVal));
                     }
                  }
               }
            #else
               // Copy data
               int k = start;
               for(int j = 0; j < cols; j++)
               {
                  for(int i = 0; i < rows; i++)
                  {
                     // Set boundary value
                     tripVal = eq.boundaryValue(compId, i, j, matIdx);
                     if(tripVal != 0.0)
                     {
                        triplets.push_back(Eigen::Triplet<typename TData::Scalar>(k, 0, tripVal));
                     }

                     // increase storage counter
                     k++;
                  }
               }
            #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

         // matIdx is the index of the slowest varying direction with multiple RHS
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

            //Safety assertion
            assert(start >= 0);

            // Copy data
            for(int j = 0; j < cols; j++)
            {
               for(int i = 0; i < rows; i++)
               {
                  // Set boundary value
                  tripVal = eq.boundaryValue(compId, i, j, matIdx);
                  if(tripVal != 0.0)
                  {
                     triplets.push_back(Eigen::Triplet<typename TData::Scalar>(i + start, j, tripVal));
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
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();

            // Copy data
            int k = start;
            for(int i = 0; i < rows; i++)
            {
               // Set boundary value
               tripVal = eq.boundaryValue(compId, i, mode(1), mode(0));
               if(tripVal != 0.0)
               {
                  triplets.push_back(Eigen::Triplet<typename TData::Scalar>(k, 0, tripVal));
               }

               // increase storage counter
               k++;
            }

         // There is a single matrix
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
         {
            assert(matIdx == 0);

            //Safety assertion
            assert(start >= 0);

            // Add source data
            int l;
            int k_;
            int j_;
            #ifdef QUICC_SPATIALDIMENSION_3D
               int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
            #else
               int dimK = 1;
            #endif //QUICC_SPATIALDIMENSION_3D
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

                     // Set boundary value
                     tripVal = eq.boundaryValue(compId, i, j, k);
                     if(tripVal != 0.0)
                     {
                        triplets.push_back(Eigen::Triplet<typename TData::Scalar>(l, 0, tripVal));
                     }
                  }
               }
            }
         }
         Eigen::SparseMatrix<typename TData::Scalar>  tmpMat(storage.rows(), storage.cols());
         tmpMat.setFromTriplets(triplets.begin(), triplets.end());
         storage = tmpMat; 
      }
   }

   template <typename TData> void setZeroNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)
      {
         //Safety assertion
         assert(start >= 0);

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            for(int k = 0; k < eq.couplingInfo(compId).galerkinN(matIdx); ++k)
            {
               // Set value to zero
               Datatypes::internal::setScalar(storage, k + start, typename TData::Scalar(0.0));
            }
         #else
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            #if defined QUICC_SPATIALSCHEME_WFT
               int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            #else
               int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            #endif //defined QUICC_SPATIALSCHEME_WFT

            // Copy data
            int k = start;
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Set field to zero
                  Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

                  // increase storage counter
                  k++;
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)
      {
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();
         int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
         #if defined QUICC_SPATIALSCHEME_WFT
            int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
         #else
            int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
         #endif //defined QUICC_SPATIALSCHEME_WFT

         //Safety assertion
         assert(start >= 0);

         // Copy data
         for(int j = zeroCol; j < cols; j++)
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Set field to zero
               Datatypes::internal::setScalar(storage, i - zeroRow + start, j - zeroCol, typename TData::Scalar(0.0));
            }
         }

      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
      {
         //Safety assertion
         assert(start >= 0);

         // Get mode indexes
         ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();
         int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);

         // Copy data
         int k = start;
         for(int i = zeroRow; i < rows; i++)
         {
            // Set field to zero
            Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

            // increase storage counter
            k++;
         }

      } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::SINGLE)
      {
         //Safety assertion
         assert(matIdx == 0);
         assert(start >= 0);

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            for(int k = 0; k < eq.couplingInfo(compId).galerkinN(matIdx); ++k)
            {
               // Set value to zero
               Datatypes::internal::setScalar(storage, k + start, typename TData::Scalar(0.0));
            }
         #else
            // Set data to zero
            int l;
            int k_;
            int j_;
            #ifdef QUICC_SPATIALDIMENSION_3D
               int dimK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
            #else
               int dimK = 1;
            #endif //QUICC_SPATIALDIMENSION_3D
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
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
      }
   }
}
}

#endif // IVECTOREQUATION_HPP
