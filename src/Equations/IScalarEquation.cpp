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

#include <iostream>
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

   void IScalarEquation::defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource)
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(2);

      // Get resolution
      pValue = PythonWrapper::makeTuple(this->unknown().dom(0).spRes()->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(this->name()).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(comp).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 1, pTmp);

      // Call model operator Python routine
      PythonWrapper::setFunction((char *)"equation_info");
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      // Get Complex solver flag
      pTmp = PyTuple_GetItem(pValue, 0);
      bool isComplex = PyObject_IsTrue(pTmp);

      // Get Implicit fields
      pTmp = PyTuple_GetItem(pValue, 1);
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >  imFields;
      PythonWrapper::getList(imFields, pTmp);

      // Get Explicit fields
      pTmp = PyTuple_GetItem(pValue, 2);
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >  exFields;
      PythonWrapper::getList(exFields, pTmp);

      // Get geometric coupling flag
      pTmp = PyTuple_GetItem(pValue, 3);
      bool hasGeometricCoupling = PyObject_IsTrue(pTmp);

      // Get block information
      pArgs = PyTuple_GetItem(pValue, 4);
      pTmp = PyTuple_GetItem(pArgs, 0);
      int blockSize = PyLong_AsLong(pTmp);
      pTmp = PyTuple_GetItem(pArgs, 1);
      int rhsSize = PyLong_AsLong(pTmp);

      // Finalise Python interpreter
      PythonWrapper::cleanup();

      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(comp,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), comp);

      // General setup: equation type? real/complex solver? start from m = ?
      infoIt.first->second.setGeneral(eqType, isComplex, iZero);

      // Set nonlinear flags: has nonlinear term? has quasi-inverse?
      infoIt.first->second.setNonlinear(hasNL, hasQI);

      // Set source flags: has source term?
      infoIt.first->second.setSource(hasSource);

      // Set index type: SLOWEST or MODE
      /// \mhdBug Always sets indexType to SLOWEST
      infoIt.first->second.setIndexType(CouplingInformation::SLOWEST);

      // Create implicit field coupling
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >::const_iterator fIt;
      for(fIt = imFields.begin(); fIt != imFields.end(); ++fIt)
      {
         infoIt.first->second.addImplicitField(fIt->first, fIt->second);
      }

      // Create explicit field
      for(fIt = exFields.begin(); fIt != exFields.end(); ++fIt)
      {
         infoIt.first->second.addExplicitField(fIt->first, fIt->second);
      }

      // Set geometric coupling information
      if(hasGeometricCoupling)
      {
         throw Exception("Geometric coupling is not yet implemented!");

      // Set field coupling information
      } else
      {
         int nMat = EigenSelector::fieldCouplingNMat(this->unknown().dom(0).spRes());
         ArrayI blockNs(nMat);
         blockNs.setConstant(blockSize);
         ArrayI rhsCols(nMat);
         rhsCols.setConstant(rhsSize);
         infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 
      }

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void  IScalarEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const bool hasBoundary) const
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(5);

      // Get resolution
      pValue = PythonWrapper::makeTuple(this->unknown().dom(0).spRes()->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get equation parameters
      std::vector<std::string> eq_names = this->eqParams().names();
      std::vector<NonDimensional::Id> eq_ids = this->eqParams().ids();
      std::vector<MHDFloat> eq_vals;
      for(unsigned int i = 0; i < eq_names.size(); i++)
      {
         eq_vals.push_back(this->eqParams().nd(eq_ids.at(i)));
      }
      pValue = PythonWrapper::makeDict(eq_names, eq_vals);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Get the eigen direction values
      pValue = PythonWrapper::makeTuple(EigenSelector::getEigs(*this, matIdx));
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      pValue = PythonWrapper::makeDict(bcMap);
      PyTuple_SetItem(pArgs, 3, pValue);

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(this->name()).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(comp).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 4, pTmp);

      // Call model operator Python routine
      PythonWrapper::setFunction(IoTools::IdToHuman::toString(opId));
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonWrapper::fillMatrix(rModelMatrix, pValue);

      // Finalise Python interpreter
      PythonWrapper::cleanup();
   }

   void IScalarEquation::setQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
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

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      pValue = PythonWrapper::makeDict(bcMap);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(this->name()).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(comp).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 3, pTmp);

      // Call model operator Python routine
      PythonWrapper::setFunction(IoTools::IdToHuman::toString(ModelOperator::QI));
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonWrapper::fillMatrix(mat, pValue);

      // Finalise Python interpreter
      PythonWrapper::cleanup();
   }

   void IScalarEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx) const
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(6);

      // Get resolution
      pValue = PythonWrapper::makeTuple(this->unknown().dom(0).spRes()->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get equation parameters
      std::vector<std::string> eq_names = this->eqParams().names();
      std::vector<NonDimensional::Id> eq_ids = this->eqParams().ids();
      std::vector<MHDFloat> eq_vals;
      for(unsigned int i = 0; i < eq_names.size(); i++)
      {
         eq_vals.push_back(this->eqParams().nd(eq_ids.at(i)));
      }
      pValue = PythonWrapper::makeDict(eq_names, eq_vals);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Get the eigen direction values
      pValue = PythonWrapper::makeTuple(EigenSelector::getEigs(*this, matIdx));
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      pValue = PythonWrapper::makeDict(bcMap);
      PyTuple_SetItem(pArgs, 3, pValue);

      // Get field row
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(this->name()).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(compId).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 4, pTmp);

      // Get field col
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(fieldId.first).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(fieldId.second).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 5, pTmp);

      // Call model operator Python routine
      PythonWrapper::setFunction(IoTools::IdToHuman::toString(ModelOperator::EXPLICIT_LINEAR));
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonWrapper::fillMatrix(mat, pValue);

      // Finalise Python interpreter
      PythonWrapper::cleanup();
   }
}
}
