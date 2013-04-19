/** \file EquationData.hpp
 *  \brief Base building block for the implementation of a time dependend evolution equation
 */

#ifndef EQUATIONDATA_HPP
#define EQUATIONDATA_HPP

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
#include "Enums/Dimensions.hpp"
#include "Enums/PhysicalNames.hpp"
#include "Enums/FieldComponents.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"
#include "Equations/IEquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Variables/VariableRequirement.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of a time dependend evolution equation
    */
   class EquationData
   {
      public:
         /**
          * \brief Simple constructor
          */
         explicit EquationData(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~EquationData();

         /**
          * @brief Get name ID of the unknown
          */
         PhysicalNames::Id name() const;

         /**
          * @brief Get scalar variable
          *
          * @param name Physical name of the field
          */
         const Datatypes::ScalarVariableType& scalar(PhysicalNames::Id name) const;

         /**
          * @brief Set scalar variable
          *
          * @param name Physical name of the field
          */
         Datatypes::ScalarVariableType& rScalar(PhysicalNames::Id name);

         /**
          * @brief Get vector variable
          *
          * @param name Physical name of the field
          */
         const Datatypes::VectorVariableType& vector(PhysicalNames::Id name) const;

         /**
          * @brief Set vector variable
          *
          * @param name Physical name of the field
          */
         Datatypes::VectorVariableType& rVector(PhysicalNames::Id name);

         /**
          * @brief Set the smart pointer to the scalar field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the scalar field
          */
         void setField(PhysicalNames::Id name, Datatypes::SharedScalarVariableType spField);

         /**
          * @brief Set the smart pointer to the vector field
          *
          * \param name Name of the field
          * \param spField Shared pointer to the vector field
          */
         void setField(PhysicalNames::Id name, Datatypes::SharedVectorVariableType spField);
   
         /**
          * @brief Get the number of implemented boundary conditions
          */
         int nBC(FieldComponents::Spectral::Id id, Dimensions::Transform::Id dim) const;

         /**
          * @brief Get the explicit linear matrices (sparse matrices)
          */
         const DecoupledZSparse& explicitLinear(FieldComponents::Spectral::Id id, const int j) const;

         /**
          * @brief Get the coupling information
          */
         const CouplingInformation&  couplingInfo(FieldComponents::Spectral::Id id) const;

         /**
          * @brief Get map of field storage requirements information
          *
          * \mhdBug Ultimatively this should depend on component
          */
         const VariableRequirement& requirements() const;

         /**
          * @brief Get map of field storage requirements information
          */
         const FieldRequirement& requirements(PhysicalNames::Id id) const;

         /**
          * @brief Get the equation parameters
          */
         const IEquationParameters& eqParams() const;

      protected:
         /**
          * @brief Set the unknown name of equation
          */
         void setName(PhysicalNames::Id name);

         /**
          * @brief Storage for the variable requirements
          *
          * \mhdBug Ultimatively this should depend on component
          */
         VariableRequirement mRequirements;

         /**
          * @brief Coupling information of the equation
          */
         std::map<FieldComponents::Spectral::Id, CouplingInformation>  mCouplingInfos;

         /**
          * @brief Map of component and explicit linear matrices 
          */
         std::map<FieldComponents::Spectral::Id, std::vector<DecoupledZSparse> > mLMatrices;

         /**
          * @brief Map of component and nonlinear term multiplication matrices
          */
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> > mNLMatrices;

      private:
         /**
          * @brief Storage for smart equation parameters
          */
         SharedIEquationParameters   mspEqParams;

         /**
          * @brief Map of name and pointer for the scalar variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>  mScalars;

         /**
          * @brief Map of name and pointer for the vector variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>  mVectors;

         /**
          * @brief Name ID of the unknown
          */
         PhysicalNames::Id mName;
   };

   /// Typedef for a smart EquationData
   typedef SharedPtrMacro<EquationData> SharedEquationData;
}
}

#endif // EQUATIONDATA_HPP
