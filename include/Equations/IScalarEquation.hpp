/**
 * @file IScalarEquation.hpp
 * @brief Base for the implementation of a scalar equation 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "Equations/EquationParameters.hpp"
#include "Equations/IEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Base for the implementation of a scalar equation
    */
   class IScalarEquation: public IEquation
   {
      public:
         /**
          * \brief Simple constructor
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
         void storeSolution(FieldComponents::Spectral::Id compId, const DecoupledZMatrix& storage, const int matIdx, const int start);
         void storeSolution(FieldComponents::Spectral::Id compId, const MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Initialise the spectral equation matrices
          *
          * It has only a semi-dummy implementation
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds);
         
      protected:
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
    * @brief Copy unknown spectral values to solver
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start);
   void copyUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start);

   /**
    * @brief Transfer nonlinear spectral values from unknown to solver
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start);
   void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start);

   /**
    * @brief Add source term
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& storage, const int matIdx, const int start);
   void addSource(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, MatrixZ& storage, const int matIdx, const int start);

}
}

#endif // ISCALAREQUATION_HPP
