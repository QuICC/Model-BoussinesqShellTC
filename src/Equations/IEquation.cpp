/** 
 * @file IEquation.cpp
 * @brief Source of building block for the implementation of a time dependend evolution equation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IEquation.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IEquation::IEquation(const std::string& pyName, SharedEquationParameters spEqParams)
      : EquationData(pyName, spEqParams)
   {
   }

   IEquation::~IEquation()
   {
   }

   void IEquation::init()
   {
      this->setCoupling();
   }

   void IEquation::initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      // Get the number of systems
      int nSystems = this->couplingInfo(compId).nSystems();

      //
      // Initialise the quasi-inverse operators for the nonlinear terms (if required)
      //
      if(this->couplingInfo(compId).hasQuasiInverse())
      {
         this->mNLMatrices.insert(std::make_pair(compId, std::vector<SparseMatrix>()));
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator qIt = this->mNLMatrices.find(compId);
         qIt->second.reserve(nSystems);
         for(int i = 0; i < nSystems; ++i)
         {
            qIt->second.push_back(SparseMatrix());

            this->setQuasiInverse(compId, qIt->second.back(), i);
         }
      }

      //
      // Initialise the explicit linear operators
      //
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = this->couplingInfo(compId).explicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         std::vector<DecoupledZSparse> tmpMat;
         tmpMat.reserve(nSystems);

         bool isComplex = false;

         // Create matrices
         for(int i = 0; i < nSystems; ++i)
         {
            // Get linear block
            tmpMat.push_back(DecoupledZSparse());
            this->setExplicitLinearBlock(compId, tmpMat.at(i), *fIt, i);

            // Explicit operator requires an additional minus sign
            tmpMat.at(i).real() = -tmpMat.at(i).real();
            tmpMat.at(i).imag() = -tmpMat.at(i).imag();

            isComplex = isComplex || (tmpMat.at(i).imag().nonZeros() > 0);
         }

         // Create key
         std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, *fIt);

         // Select real or complex operator
         if(isComplex)
         {
            this->mLZMatrices.insert(std::make_pair(key, std::vector<SparseMatrixZ>()));
            this->mLZMatrices.find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               SparseMatrixZ tmp = tmpMat.at(i).real().cast<MHDComplex>() + Math::cI*tmpMat.at(i).imag();
               this->mLZMatrices.find(key)->second.push_back(tmp);
            }
         } else
         {
            this->mLDMatrices.insert(std::make_pair(key, std::vector<SparseMatrix>()));
            this->mLDMatrices.find(key)->second.back().reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->mLDMatrices.find(key)->second.push_back(SparseMatrix());

               this->mLDMatrices.find(key)->second.back().swap(tmpMat.at(i).real());
            }
         }
      }
   }

   void IEquation::dispatchCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasQI, const bool hasSource, const SharedResolution spRes)
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(2);

      // Get resolution
      pValue = PythonWrapper::makeTuple(spRes->sim()->dimensions(Dimensions::Space::SPECTRAL));
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

      // Get geometric coupling flag
      pTmp = PyTuple_GetItem(pValue, 4);
      int indexMode = PyLong_AsLong(pTmp);

      // Get block information
      pArgs = PyTuple_GetItem(pValue, 5);
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
      infoIt.first->second.setIndexType(static_cast<CouplingInformation::IndexType>(indexMode));

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
         int nMat = EigenSelector::fieldCouplingNMat(spRes);
         ArrayI blockNs(nMat);
         blockNs.setConstant(blockSize);
         ArrayI rhsCols(nMat);
         rhsCols.setConstant(rhsSize);
         infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 
      }

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void  IEquation::dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue, *pList;
      pArgs = PyTuple_New(5);

      // Get resolution
      pValue = PythonWrapper::makeTuple(spRes->sim()->dimensions(Dimensions::Space::SPECTRAL));
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
      pValue = PythonWrapper::makeTuple(eigs);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      bcMap.insert(std::make_pair("bcType", bcType));
      pValue = PythonWrapper::makeDict(bcMap);
      PyTuple_SetItem(pArgs, 3, pValue);

      // Get list of implicit fields
      CouplingInformation::FieldId_range impRange = this->couplingInfo(comp).implicitRange();
      pList = PyList_New(0);
      for(CouplingInformation::FieldId_iterator fIt = impRange.first; fIt != impRange.second; ++fIt)
      {
         pTmp = PyTuple_New(2);
         pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(fIt->first).c_str());
         PyTuple_SetItem(pTmp, 0, pValue);
         pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(fIt->second).c_str());
         PyTuple_SetItem(pTmp, 1, pValue);

         PyList_Append(pList, pTmp);
      }
      PyTuple_SetItem(pArgs, 4, pList);

      // Call model operator Python routine
      PythonWrapper::setFunction(IoTools::IdToHuman::toString(opId));
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonWrapper::fillMatrix(rModelMatrix, pValue);

      // Finalise Python interpreter
      PythonWrapper::cleanup();
   }

   void IEquation::dispatchQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(4);

      // Get resolution
      pValue = PythonWrapper::makeTuple(spRes->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get the eigen direction values
      pValue = PythonWrapper::makeTuple(eigs);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      bcMap.insert(std::make_pair("bcType", ModelOperatorBoundary::NO_BC));
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

   void IEquation::dispatchExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(this->pyName());

      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(6);

      // Get resolution
      pValue = PythonWrapper::makeTuple(spRes->sim()->dimensions(Dimensions::Space::SPECTRAL));
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
      pValue = PythonWrapper::makeTuple(eigs);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      bcMap.insert(std::make_pair("bcType", ModelOperatorBoundary::NO_BC));
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

   void IEquation::setQuasiInverse(FieldComponents::Spectral::Id comp, SparseMatrix &mat, const int matIdx) const
   {
      // This implementation should never get called!
      throw Exception("Called dummy implementation of setQuasiInverse!");
   }

   void IEquation::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int matIdx) const
   {
      // This implementation should never get called!
      throw Exception("Called dummy implementation of setExplicitLinearBlock!");
   }

   void IEquation::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // This implementation should never get called!
      throw Exception("Activated nonlinear term without implementation!");
   }

   Datatypes::SpectralScalarType::PointType IEquation::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // This implementation should never get called!
      throw Exception("Activated source term without implementation!");

      return Datatypes::SpectralScalarType::PointType();
   }

   void  IEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id comp, const int matIdx, const ModelOperatorBoundary::Id bcType) const
   {
      // This implementation should never get called!
      throw Exception("Called dummy implementation of buildModelMatrix!");
   }
}
}
