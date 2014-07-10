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
#include "Enums/ModelOperatorBoundary.hpp"
#include "Equations/EquationParameters.hpp"
#include "Equations/IEquation.hpp"
#include "Base/DecoupledComplexInternal.hpp"

namespace GeoMHDiSCC {

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
          */
         void updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId);

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
         void defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource);

         /**
          * @brief Set the unknown variable
          */
         Datatypes::VectorVariableType& rUnknown();

         /**
          * @brief List of the avaiable vector components
          */
         std::vector<FieldComponents::Spectral::Id>   mSpectralIds;

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Datatypes::SharedVectorVariableType mspUnknown;
   };

   /// Typedef for a shared IVectorEquation
   typedef SharedPtrMacro<IVectorEquation> SharedIVectorEquation;

   /**
    * @brief Copy unknown spectral values to solver
    *
    * @param eq         Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    * @param useShift   Use galerkin shifts
    */
   template <typename TData> void copyUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift);

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

      if(this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = this->unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

         // Copy data
         int k = solStart;
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy timestep output into field
               this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(Datatypes::internal::getScalar(*solution, k),i,j,matIdx);

               // increase linear storage counter
               k++;
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
            this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(Datatypes::internal::getScalar(*solution,k),i,mode(1),mode(0));

            // increase linear storage counter
            k++;
         }
      }
   }

   template <typename TData> void copyUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift)
   {
      int zeroRow = 0;
      int zeroCol = 0;
      if(useShift)
      {
         zeroRow = eq.couplingInfo(compId).galerkinShift(0);
         zeroCol = eq.couplingInfo(compId).galerkinShift(1);
      }

      // matIdx is the index of the slowest varying direction
      if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
      {
         int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
         int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);

         // Copy data
         int k = start;
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
         for(int i = zeroRow; i < rows; i++)
         {
            // Copy field value into storage
            Datatypes::internal::setScalar(storage, k, eq.unknown().dom(0).perturbation().comp(compId).point(i,mode(1),mode(0)));

            // increase storage counter
            k++;
         }
      }
   }

   template <typename TData> void copyNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Check if a nonlinear computation took place and a quasi-inverse has to be applied
      if(eq.couplingInfo(compId).hasNonlinear())
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
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();
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
                  // Set field to zero
                  Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

                  // increase storage counter
                  k++;
               }
            }
         } else if(eq.couplingInfo(compId).indexType() == CouplingInformation::MODE)
         {
            //Safety assertion
            assert(start >= 0);

            // Get mode indexes
            ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(mode(0)).rows();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);

            // Copy data
            int k = start;
            for(int i = zeroRow; i < rows; i++)
            {
               // Set field to zero
               Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

               // increase storage counter
               k++;
            }
         }
      }
   }

   template <typename TData> void addSource(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         // matIdx is the index of the slowest varying direction
         if(eq.couplingInfo(compId).indexType() == CouplingInformation::SLOWEST)
         {
            int rows = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).rows();
            int cols = eq.unknown().dom(0).perturbation().comp(compId).slice(matIdx).cols();
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
                  // Add source term
                  Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, j, matIdx));

                  // increase storage counter
                  k++;
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
            int zeroRow = eq.couplingInfo(compId).galerkinShift(0);

            // Copy data
            int k = start;
            for(int i = zeroRow; i < rows; i++)
            {
               // Add source term
               Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, mode(1), mode(0)));

               // increase storage counter
               k++;
            }
         }
      }
   }
}
}

#endif // IVECTOREQUATION_HPP
