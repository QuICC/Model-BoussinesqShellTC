/** 
 * @file IVectorEquation.cpp
 * @brief Source of the base implementation of a vector equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <Python.h>

// External includes
//

// Class include
//
#include "Equations/IVectorEquation.hpp"

// Project includes
//
#include "TypeSelectors/EquationEigenSelector.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IVectorEquation::IVectorEquation(const std::string& pyName, SharedEquationParameters spEqParams)
      : IEquation(pyName, spEqParams)
   {
   }

   IVectorEquation::~IVectorEquation()
   {
   }

   void IVectorEquation::setUnknown(Datatypes::SharedVectorVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::VectorVariableType& IVectorEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::VectorVariableType& IVectorEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   SharedResolution IVectorEquation::spRes() const
   {
      return this->unknown().dom(0).spRes();
   }

   int IVectorEquation::nSpectral() const
   {
      return this->mSpectralIds.size();
   }

   IVectorEquation::SpectralComponent_range IVectorEquation::spectralRange() const
   {
      return std::make_pair(this->mSpectralIds.begin(), this->mSpectralIds.end());
   }

   void IVectorEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId)
   {
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().cols() < rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().rComp(compId).setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().rComp(compId).data().rows()));
   }

   void IVectorEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      for(SpectralComponent_iterator it = this->mSpectralIds.begin(); it != this->mSpectralIds.end(); ++it)
      {
         // Make sure it is safe to do nothing
         bool needInit = this->couplingInfo(*it).hasQuasiInverse();

         CouplingInformation::FieldId_range fRange = this->couplingInfo(*it).explicitRange();
         needInit = needInit && (fRange.first == fRange.second);

         // Initialise spectral matrices
         if(needInit)
         {
            this->initSpectralMatricesComponent(spBcIds, *it);
         }
      }
   }

   const Boundary::CoordinatorSelector& IVectorEquation::bcCoord(FieldComponents::Spectral::Id compId) const
   {
      return this->mBcCoord.at(compId);
   }

   Boundary::CoordinatorSelector& IVectorEquation::rBcCoord(FieldComponents::Spectral::Id compId)
   {
      return this->mBcCoord.at(compId);
   }

   void  IVectorEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const bool hasBoundary) const
   {
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

      // Get boundary condition
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
