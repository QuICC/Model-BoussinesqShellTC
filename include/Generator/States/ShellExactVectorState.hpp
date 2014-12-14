/**
 * @file ShellExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector states in a shell 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHELLEXACTVECTORSTATE_HPP
#define SHELLEXACTVECTORSTATE_HPP

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
#include "Generator/States/ShellExactStateIds.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact vector state in a shell
    */
   class ShellExactVectorState: public IVectorEquation
   {
      public:
         /// Typedef to simplify notations for harmonic mode
         typedef std::tr1::tuple<int, int, MHDComplex> HarmonicModeType;

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         ShellExactVectorState(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ShellExactVectorState();

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
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set the state type id
          */
         void setStateType(const FieldComponents::Physical::Id compId, const ShellExactStateIds::Id id);

         /**
          * @brief Set options for the harmonics states
          *
          * @param modes   List of harmonics with amplitude to create
          */
         void setHarmonicOptions(const FieldComponents::Physical::Id compId, const std::vector<HarmonicModeType>& modes);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Type of the state to generate
          */
         std::map<FieldComponents::Physical::Id,ShellExactStateIds::Id> mTypeId;

         /**
          * @brief Storage for the list of spherical harmonic modes to generate
          */
         std::map<FieldComponents::Physical::Id,std::vector<HarmonicModeType> > mSHModes;
   };

   /// Typedef for a shared ShellExactVectorState
   typedef SharedPtrMacro<ShellExactVectorState> SharedShellExactVectorState;

}
}

#endif // SHELLEXACTVECTORSTATE_HPP
