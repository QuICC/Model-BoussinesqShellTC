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

   void IVectorEquation::defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource)
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

      // Set mininal matrix coupling
      int nMat = 0;
      ArrayI blockNs;
      ArrayI rhsCols;
      EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void  IVectorEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const bool hasBoundary) const
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

      // Get boundary condition
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

   void IVectorEquation::setQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
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

   void IVectorEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx) const
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
