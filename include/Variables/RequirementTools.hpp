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
          * @param varInfo       Variable requirements
          * @param rScalarVars   Scalar variables
          * @param rVectorVars   Vector variables
          * @param scalarEqs     Scalar equations
          * @param vectorEqs     Vector equations
          * @param spRes      Shared resolution
          */
         static void initVariables(VariableRequirement& varInfo, std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& rScalarVars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& rVectorVars, const std::vector<Equations::SharedIScalarEquation>& scalarEqs, const std::vector<Equations::SharedIVectorEquation>& vectorEqs, SharedResolution spRes);

         /**
          * @brief Map variables to the corresponding equation 
          *
          * @param nonInfo    Nonlinear computations requirements
          * @param rScalarEqs Scalar equations
          * @param rVectorEqs Vector equations
          * @param scalarVars Scalar variables
          * @param vectorVars Vector variables
          */
         static void mapEquationVariables(std::set<PhysicalNames::Id>& nonInfo, std::vector<Equations::SharedIScalarEquation>& rScalarEqs, std::vector<Equations::SharedIVectorEquation>& rVectorEqs, const std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalarVars, const std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectorVars);
         
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
