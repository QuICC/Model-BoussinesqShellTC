/** 
 * @file RequirementTools.hpp
 * @brief Implementation of requirement tools for equations and variables
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef REQUIREMENTTOOLS_HPP
#define REQUIREMENTTOOLS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Variables/VariableRequirement.hpp"
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/TransformTree.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of requirement tools for equations and variables
    */
   class RequirementTools
   {
      public:
         /**
          * @brief Initialise variables and variable requirements from equations
          *
          * @param projectorTree Transform tree for backward projection
          * @param rScalarVars   Scalar variables
          * @param rVectorVars   Vector variables
          * @param scalarEqs     Scalar equations
          * @param vectorEqs     Vector equations
          * @param spRes      Shared resolution
          */
         static void initVariables(std::vector<Transform::TransformTree>& projectorTree, std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& rScalarVars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& rVectorVars, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs, SharedResolution spRes);

         /**
          * @brief Map variables to the corresponding equation 
          *
          * @param integratorTree Transform tree for forward integration
          * @param rScalarEqs          Scalar equations
          * @param rVectorEqs          Vector equations
          * @param scalarVars          Scalar variables
          * @param vectorVars          Vector variables
          * @param forwardIsNonlinear  Forward transform works on nonlinear terms
          */
         static void mapEquationVariables(std::vector<Transform::TransformTree>& integratorTree, std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalarVars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectorVars, const bool forwardIsNonlinear, const SharedSimulationBoundary spBcs);

         /**
          * @brief Initialise imposed variables and variable requirements from equations
          *
          * @param projectorTree Transform tree for backward projection
          * @param rScalarVars   Scalar variables
          * @param rVectorVars   Vector variables
          * @param scalarEqs     Scalar equations
          * @param vectorEqs     Vector equations
          * @param spRes      Shared resolution
          */
         static void initImposedVariables(std::vector<Transform::TransformTree>& projectorTree, std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& rScalarVars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& rVectorVars, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs, SharedResolution spRes);

         /**
          * @brief Map imposed variables to the corresponding equation 
          *
          * @param rScalarEqs          Scalar equations
          * @param rVectorEqs          Vector equations
          * @param scalarVars          Scalar variables
          * @param vectorVars          Vector variables
          */
         static void mapImposedVariables(std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalarVars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectorVars);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         RequirementTools();

         /**
          * @brief Destructor
          */
         ~RequirementTools();
   };

}

#endif // REQUIREMENTTOOLS_HPP
