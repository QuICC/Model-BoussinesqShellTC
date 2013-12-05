/**
 * @file SphereExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector states in a sphere 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHEREEXACTVECTORSTATE_HPP
#define SPHEREEXACTVECTORSTATE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>
#include <tr1/tuple>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact vector state in a sphere
    */
   class SphereExactVectorState: public IVectorEquation
   {
      public:
         /**
          * @brief Enums for the avaialable exact states
          */
         enum StateTypeId {
            CONSTANT,
            HARMONIC,
         };

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @paarm name       Name of the field
          */
         SphereExactVectorState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~SphereExactVectorState();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const;

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDComplex sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Generic operator row dispatcher
          */
         virtual DecoupledZSparse operatorRow(const OperatorRowId opId, FieldComponents::Spectral::Id comp, const int matIdx, const bool hasBoundary) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set the state type id
          */
         void setStateType(const SphereExactVectorState::StateTypeId id);

         /**
          * @brief Set options for the harmonics states
          *
          * @param rModes   List of harmonics with amplitude to create
          * @param tModes   List of harmonics with amplitude to create
          * @param pModes   List of harmonics with amplitude to create
          */
         void setHarmonicOptions(const std::vector<std::tr1::tuple<int, int, MHDComplex> >& rModes, const std::vector<std::tr1::tuple<int, int, MHDComplex> >& tModes, const std::vector<std::tr1::tuple<int, int, MHDComplex> >& pModes);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the quasi inverse matrix operator
          */
         virtual void setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix &mat) const;

         /**
          * @brief Set the explicit linear matrix operator
          */
         virtual void setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const;

      private:
         /**
          * @brief Type of the state to generate
          */
         StateTypeId mTypeId;

         /**
          * @brief Storage for the list of spherical harmonic modes to generate
          */
         std::vector<std::tr1::tuple<int,int,MHDComplex> > mRSHModes;

         /**
          * @brief Storage for the list of spherical harmonic modes to generate
          */
         std::vector<std::tr1::tuple<int,int,MHDComplex> > mTSHModes;

         /**
          * @brief Storage for the list of spherical harmonic modes to generate
          */
         std::vector<std::tr1::tuple<int,int,MHDComplex> > mPSHModes;
   };

   /// Typedef for a shared SphereExactVectorState
   typedef SharedPtrMacro<SphereExactVectorState> SharedSphereExactVectorState;

   /**
    * @brief Get the quasi-inverse matrix operator
    *
    * @param eq      Equation to work on
    * @param mat     Storage for output matrix
    */
   void quasiInverseBlock(const SphereExactVectorState& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat);

   /**
    * @brief Get the linear matrix block for an equation on given field
    *
    * @param mat     Storage for output matrix
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void linearBlock(const SphereExactVectorState& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary);

   /**
    * @brief Get the boundary condition matrix block for an equation on given field
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param fieldId Physical ID of the field
    * @parm eigs     Wave number k
    */
   void boundaryBlock(SphereExactVectorState& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx);

}
}

#endif // SPHEREEXACTVECTORSTATE_HPP
