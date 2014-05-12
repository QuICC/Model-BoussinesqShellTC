/** 
 * @file IScalarEquation.cpp
 * @brief Source of the base implementation of a scalar equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <Python.h>
#include <iostream>

// External includes
//

// Class include
//
#include "Equations/IScalarEquation.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IScalarEquation::IScalarEquation(const std::string& pyName, SharedEquationParameters spEqParams)
      : IEquation(pyName, spEqParams)
   {
   }

   IScalarEquation::~IScalarEquation()
   {
   }

   void IScalarEquation::setUnknown(Datatypes::SharedScalarVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::ScalarVariableType& IScalarEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::ScalarVariableType& IScalarEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   SharedResolution IScalarEquation::spRes() const
   {
      return this->unknown().dom(0).spRes();
   }

   void IScalarEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().data().cols() == rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
   }

   void IScalarEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Make sure it is safe to do nothing
      bool needInit = this->couplingInfo(FieldComponents::Spectral::SCALAR).hasQuasiInverse();

      CouplingInformation::FieldId_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange();
      needInit = needInit && (fRange.first == fRange.second);

      // Initialise spectral matrices
      if(needInit)
      {
         this->initSpectralMatricesComponent(spBcIds, FieldComponents::Spectral::SCALAR);
      }
   }

   const Boundary::CoordinatorSelector& IScalarEquation::bcCoord(FieldComponents::Spectral::Id compId) const
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      return this->mBcCoord;
   }

   Boundary::CoordinatorSelector& IScalarEquation::rBcCoord(FieldComponents::Spectral::Id compId)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      return this->mBcCoord;
   }

   void  IScalarEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const bool hasBoundary) const
   {
      std::cerr << "CALLING PYTHON FOR SCALAR" << std::endl;

      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(4);

      // Get resolution
      pValue = PythonWrapper::makeTuple(this->unknown().dom(0).spRes()->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get the eigen direction values
      pValue = PythonWrapper::makeTuple(EigenSelector::getEigs(*this, matIdx));
      PyTuple_SetItem(pArgs, 1, pValue);

      // Get equation parameters
      std::vector<std::string> eq_names = this->eqParams().names();
      std::vector<NonDimensional::Id> eq_ids = this->eqParams().ids();
      std::vector<MHDFloat> eq_vals;
      for(unsigned int i = 0; i < eq_names.size(); i++)
      {
         eq_vals.push_back(this->eqParams().nd(eq_ids.at(i)));
      }
      pValue = PythonWrapper::makeDict(eq_names, eq_vals);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get boundary conditions
      pTmp = PyTuple_New(1);
      pValue = PyLong_FromLong(1);
      PyTuple_SetItem(pTmp, 0, pValue);
      PyTuple_SetItem(pArgs, 3, pTmp);

      // Call model operator Python routine
      PythonWrapper::setFunction(IoTools::IdToHuman::toString(opId));
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonWrapper::fillMatrix(rModelMatrix, pValue);

      // Finalise Python interpreter
      PythonWrapper::cleanup();
   }
}
}
