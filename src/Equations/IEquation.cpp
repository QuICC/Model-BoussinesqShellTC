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
#include "Python/PythonTools.hpp"
#include "Python/PythonModelWrapper.hpp"
#include "Base/MathConstants.hpp"
#include "IoTools/IdToHuman.hpp"

namespace QuICC {

namespace Equations {

   IEquation::IEquation(SharedEquationParameters spEqParams)
      : EquationData(spEqParams)
   {
   }

   IEquation::~IEquation()
   {
   }

   void IEquation::init(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Set the coupling
      this->setCoupling();

      // Add the nonlinear integration components
      this->setNLComponents();
   }

   void IEquation::initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      //
      // Initialise the galerkin stencils (if activated and required)
      //
      #ifdef QUICC_BOUNDARYMETHOD_GALERKIN
      // Get the number of systems
      int nSystems = this->couplingInfo(compId).nSystems();

      if(this->couplingInfo(compId).isGalerkin())
      {
         this->mGStencils.insert(std::make_pair(compId, std::vector<SparseMatrix>()));
         std::map<FieldComponents::Spectral::Id, std::vector<SparseMatrix> >::iterator sIt = this->mGStencils.find(compId);
         sIt->second.reserve(nSystems);
         for(int i = 0; i < nSystems; ++i)
         {
            sIt->second.push_back(SparseMatrix());

            this->setGalerkinStencil(compId, sIt->second.back(), i);
         }
      }
      #endif //QUICC_BOUNDARYMETHOD_GALERKIN

      //
      // Initialise quasi inverse operator
      //
      this->initQIMatrices(spBcIds, compId);

      //
      // Initialise the explicit linear operators
      //
      this->initExplicitMatrices(spBcIds, compId, ModelOperator::EXPLICIT_LINEAR);

      //
      // Initialise the explicit nonlinear operators
      //
      this->initExplicitMatrices(spBcIds, compId, ModelOperator::EXPLICIT_NONLINEAR);

      //
      // Initialise the explicit nextstep operators
      //
      this->initExplicitMatrices(spBcIds, compId, ModelOperator::EXPLICIT_NEXTSTEP);
   }

   void IEquation::initQIMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId)
   {
      if(this->couplingInfo(compId).hasQuasiInverse())
      {
         // Get the number of systems
         int nSystems = this->couplingInfo(compId).nSystems();

         //
         // Initialise the quasi inverse operators
         //
         SpectralFieldId  fieldId = std::make_pair(this->name(), compId);

         std::vector<DecoupledZSparse> tmpMat;
         tmpMat.reserve(nSystems);

         bool isComplex = false;

         // Create matrices
         for(int i = 0; i < nSystems; ++i)
         {
            // Get block
            tmpMat.push_back(DecoupledZSparse());
            this->setExplicitBlock(compId, tmpMat.at(i), ModelOperator::EXPLICIT_NONLINEAR, fieldId, i);

            isComplex = isComplex || (tmpMat.at(i).imag().nonZeros() > 0);
         }

         // Select real or complex operator
         if(isComplex)
         {
            this->mQIZMatrices.insert(std::make_pair(compId, std::vector<SparseMatrixZ>()));
            this->mQIZMatrices.find(compId)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               SparseMatrixZ tmp = tmpMat.at(i).real().cast<MHDComplex>() + Math::cI*tmpMat.at(i).imag();
               this->mQIZMatrices.find(compId)->second.push_back(tmp);
            }
         } else
         {
            this->mQIDMatrices.insert(std::make_pair(compId, std::vector<SparseMatrix>()));
            this->mQIDMatrices.find(compId)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->mQIDMatrices.find(compId)->second.push_back(SparseMatrix());

               this->mQIDMatrices.find(compId)->second.back().swap(tmpMat.at(i).real());
            }
         }
      }
   }

   void IEquation::initExplicitMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId, const ModelOperator::Id opId)
   {
      // Get the number of systems
      int nSystems = this->couplingInfo(compId).nSystems();

      //
      // Initialise the explicit operators
      //
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = this->couplingInfo(compId).explicitRange(opId);
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
            this->setExplicitBlock(compId, tmpMat.at(i), opId, *fIt, i);

            isComplex = isComplex || (tmpMat.at(i).imag().nonZeros() > 0);
         }

         // Create key
         std::pair<FieldComponents::Spectral::Id, SpectralFieldId>   key = std::make_pair(compId, *fIt);

         // Select real or complex operator
         if(isComplex)
         {
            this->rEZMatrices(opId).insert(std::make_pair(key, std::vector<SparseMatrixZ>()));
            this->rEZMatrices(opId).find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               SparseMatrixZ tmp = tmpMat.at(i).real().cast<MHDComplex>() + Math::cI*tmpMat.at(i).imag();
               this->rEZMatrices(opId).find(key)->second.push_back(tmp);
            }
         } else
         {
            this->rEDMatrices(opId).insert(std::make_pair(key, std::vector<SparseMatrix>()));
            this->rEDMatrices(opId).find(key)->second.reserve(nSystems);

            for(int i = 0; i < nSystems; ++i)
            {
               this->rEDMatrices(opId).find(key)->second.push_back(SparseMatrix());

               this->rEDMatrices(opId).find(key)->second.back().swap(tmpMat.at(i).real());
            }
         }
      }
   }

   void IEquation::dispatchCoupling(FieldComponents::Spectral::Id compId, CouplingInformation::EquationTypeId eqType, const int iZero, const bool hasNL, const bool hasSource, const SharedResolution spRes, const bool allowExplicit)
   {
      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(2);

      // Get resolution
      pValue = PythonTools::makeTuple(spRes->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(this->name()).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(compId).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 1, pTmp);

      // Call model operator Python routine
      PythonModelWrapper::setMethod((char *)"equation_info");
      pValue = PythonModelWrapper::callMethod(pArgs);
      Py_DECREF(pArgs);

      // Get Complex solver flag
      pTmp = PyTuple_GetItem(pValue, 0);
      bool isComplex = PyObject_IsTrue(pTmp);

      // Get Implicit fields
      pTmp = PyTuple_GetItem(pValue, 1);
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >  imFields;
      PythonTools::getList(imFields, pTmp);

      // Get Explicit linear fields
      pTmp = PyTuple_GetItem(pValue, 2);
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >  exLFields;
      PythonTools::getList(exLFields, pTmp);

      // Get Explicit nonlinear fields
      pTmp = PyTuple_GetItem(pValue, 3);
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >  exNLFields;
      PythonTools::getList(exNLFields, pTmp);

      // Get Explicit nextstep fields
      pTmp = PyTuple_GetItem(pValue, 4);
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >  exNSFields;
      PythonTools::getList(exNSFields, pTmp);

      // Get index mode
      pTmp = PyTuple_GetItem(pValue, 5);
      int indexMode = PyLong_AsLong(pTmp);

      // Finalise Python interpreter
      Py_DECREF(pValue);

      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(compId,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), compId);

      // Compute effective starting index for local CPU
      int cpuIZero = iZero;
      if(iZero == 1)
      {
         if(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0) == 0 && spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0) == 0)
         {
            cpuIZero = 1;
         } else
         {
            cpuIZero = 0;
         }
      } else if(iZero > 1)
      {
         throw Exception("Matrix starting index > 1 is not implemented yet!");
      }

      // General setup: equation type? real/complex solver? start from m = ?
      infoIt.first->second.setGeneral(eqType, isComplex, cpuIZero);

      // Set source flags: has source term?
      infoIt.first->second.setSource(hasSource);

      // Set index type: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
      infoIt.first->second.setIndexType(static_cast<CouplingInformation::IndexType>(indexMode));

      // Create implicit field coupling
      std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> >::const_iterator fIt;
      for(fIt = imFields.begin(); fIt != imFields.end(); ++fIt)
      {
         infoIt.first->second.addImplicitField(fIt->first, fIt->second);
      }

      // Create explicit fields
      bool hasQI = false;
      if(allowExplicit)
      {
         // explicit linear
         for(fIt = exLFields.begin(); fIt != exLFields.end(); ++fIt)
         {
            infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::EXPLICIT_LINEAR);
         }

         // explicit nonlinear
         for(fIt = exNLFields.begin(); fIt != exNLFields.end(); ++fIt)
         {
            if(!(fIt->first == this->name() && fIt->second == compId))
            {
               infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::EXPLICIT_NONLINEAR);
            }
         }

         // explicit nextstep
         for(fIt = exNSFields.begin(); fIt != exNSFields.end(); ++fIt)
         {
            infoIt.first->second.addExplicitField(fIt->first, fIt->second, ModelOperator::EXPLICIT_NEXTSTEP);
         }

         // Extract quasi inverse
         fIt = std::find(exNLFields.begin(), exNLFields.end(), std::make_pair(this->name(), compId));
         if(fIt != exNLFields.end())
         {
            hasQI = true;
         }
      }

      // Set nonlinear flags: has nonlinear term? has quasi-inverse?
      infoIt.first->second.setNonlinear(hasNL, hasNL && hasQI);

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, eqId.second);

      // Get number of matrices
      int nMat = infoIt.first->second.eigenTools().nMat(spRes);

      // Prepare Python call
      pArgs = PyTuple_New(4);

      // Get resolution
      pValue = PythonTools::makeTuple(spRes->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      pValue = PythonTools::makeDict(bcMap);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(this->name()).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(compId).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 3, pTmp);

      // Call model operator Python routine
      PythonModelWrapper::setMethod((char *)"operator_info");

      // Set field coupling information
      ArrayI tauNs(nMat);
      ArrayI galerkinNs(nMat);
      MatrixI galerkinShifts(nMat, 3);
      ArrayI rhsCols(nMat);
      ArrayI systemNs(nMat);

      // Loop overall matrices/eigs
      for(int matIdx = 0; matIdx < nMat; ++matIdx)
      {
         // Get the eigen direction values
         pValue = PythonTools::makeTuple(infoIt.first->second.eigenTools().getEigs(spRes, matIdx));
         PyTuple_SetItem(pArgs, 1, pValue);

         pValue = PythonModelWrapper::callMethod(pArgs);

         // Get block information
         pTmp = PyTuple_GetItem(pValue, 0);
         tauNs(matIdx) = PyLong_AsLong(pTmp);
         pTmp = PyTuple_GetItem(pValue, 1);
         galerkinNs(matIdx) = PyLong_AsLong(pTmp);
         pTmp = PyTuple_GetItem(pValue, 2);
         ArrayI shifts(3);
         galerkinShifts(matIdx,0) = PyLong_AsLong(PyTuple_GetItem(pTmp, 0));
         galerkinShifts(matIdx,1) = PyLong_AsLong(PyTuple_GetItem(pTmp, 1));
         galerkinShifts(matIdx,2) = PyLong_AsLong(PyTuple_GetItem(pTmp, 2));
         pTmp = PyTuple_GetItem(pValue, 3);
         rhsCols(matIdx) = PyLong_AsLong(pTmp);
         pTmp = PyTuple_GetItem(pValue, 4);
         systemNs(matIdx) = PyLong_AsLong(pTmp);

         Py_DECREF(pValue);
      }
      Py_DECREF(pArgs);

      infoIt.first->second.eigenTools().setTauN(tauNs, spRes);
      infoIt.first->second.eigenTools().setGalerkinN(galerkinNs, spRes);
      infoIt.first->second.eigenTools().setRhsN(rhsCols, spRes);
      infoIt.first->second.eigenTools().setSystemN(systemNs, spRes);
      infoIt.first->second.setSizes(nMat, tauNs, galerkinNs, galerkinShifts, rhsCols, systemNs); 

      // Finalise Python interpreter
      PythonModelWrapper::cleanup();
   }

   PyObject*  IEquation::dispatchBaseArguments(const int tupleSize, const ModelOperatorBoundary::Id bcType, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const
   {
      // Prepare Python call arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(tupleSize);

      // Get resolution
      pValue = PythonTools::makeTuple(spRes->sim()->dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get equation parameters
      std::vector<std::string> eq_names = this->eqParams().names();
      std::vector<NonDimensional::Id> eq_ids = this->eqParams().ids();
      std::vector<MHDFloat> eq_vals;
      for(unsigned int i = 0; i < eq_names.size(); i++)
      {
         eq_vals.push_back(this->eqParams().nd(eq_ids.at(i)));
      }
      pValue = PythonTools::makeDict(eq_names, eq_vals);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Get the eigen direction values
      pValue = PythonTools::makeTuple(eigs);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get boundary conditions
      std::map<std::string,int> bcMap = this->bcIds().getTagMap();
      bcMap.insert(std::make_pair("bcType", bcType));
      pValue = PythonTools::makeDict(bcMap);
      PyTuple_SetItem(pArgs, 3, pValue);

      return pArgs;
   }

   void  IEquation::dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id compId, const int matIdx, const ModelOperatorBoundary::Id bcType, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const
   {
      // Get first four standard arguments in a tuple of size 6 (7 for inhomogeneous boundary condition)
      PyObject *pArgs;
      int argStart;
      if(opId == ModelOperator::INHOMOGENEOUS)
      {
         pArgs = this->dispatchBaseArguments(7, bcType, spRes, eigs);
         argStart = 5;
      } else
      {
         pArgs = this->dispatchBaseArguments(6, bcType, spRes, eigs);
         argStart = 4;
      }

      // Prepare Python call arguments
      PyObject *pTmp, *pValue, *pList;

      // Add list of modes
      if(opId == ModelOperator::INHOMOGENEOUS)
      {
         if((this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_MULTI_RHS)/* ||
        		 (this->couplingInfo(compId).indexType() == CouplingInformation::SLOWEST_SINGLE_RHS)*/)
         {
            pTmp = PyList_New(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(matIdx));
            for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(matIdx); ++i)
            {
               PyList_SetItem(pTmp, i, PyFloat_FromDouble(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(i,matIdx)));
            }
            PyTuple_SetItem(pArgs, 4, pTmp);

         } else
         {
            Py_INCREF(Py_None);
            PyTuple_SetItem(pArgs, 4, Py_None);
         }
      }

      // Get list of implicit fields
      CouplingInformation::FieldId_range impRange = this->couplingInfo(compId).implicitRange();
      pList = PyList_New(0);
      for(CouplingInformation::FieldId_iterator fIt = impRange.first; fIt != impRange.second; ++fIt)
      {
         pTmp = PyTuple_New(2);
         pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(fIt->first).c_str());
         PyTuple_SetItem(pTmp, 0, pValue);
         pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(fIt->second).c_str());
         PyTuple_SetItem(pTmp, 1, pValue);

         PyList_Append(pList, pTmp);
         Py_DECREF(pTmp);
      }
      PyTuple_SetItem(pArgs, argStart, pList);

      // Set the restriction option
      #ifdef QUICC_MPI
         #ifdef QUICC_MPISPSOLVE
            std::vector<int> slow;
            std::vector<std::vector<int> > middle;

            spRes->buildRestriction(slow, middle, matIdx);
            PyObject *pSlow, *pMiddle;

            if(middle.size() > 0)
            {
               pSlow = PythonTools::makeList(slow);
               pMiddle = PythonTools::makeList(middle);
               pTmp = PyTuple_New(2);
               PyTuple_SetItem(pTmp, 0, pSlow);
               PyTuple_SetItem(pTmp, 1, pMiddle);
            } else
            {
               pTmp = PythonTools::makeList(slow);
            }

            PyTuple_SetItem(pArgs, argStart+1, pTmp);
         #else
            // MPI code with serial sparse solver
            Py_INCREF(Py_None);
            PyTuple_SetItem(pArgs, argStart+1, Py_None);
         #endif //QUICC_MPISPSOLVE
      #else
         // Serial code can't have any restriction
         Py_INCREF(Py_None);
         PyTuple_SetItem(pArgs, argStart+1, Py_None);
      #endif //QUICC_MPI

      // Call model operator Python routine
      PythonModelWrapper::setMethod(IoTools::IdToHuman::toString(opId));
      pValue = PythonModelWrapper::callMethod(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonModelWrapper::fillMatrix(rModelMatrix, pValue);
      Py_DECREF(pValue);

      // Finalise Python interpreter
      PythonModelWrapper::cleanup();
   }

   void IEquation::dispatchGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs, const bool makeSquare) const
   {
      // Get first four standard arguments in a tuple of size 5
      PyObject *pArgs = this->dispatchBaseArguments(6, ModelOperatorBoundary::STENCIL, spRes, eigs);

      // Prepare Python call arguments
      PyObject *pTmp, *pValue;

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(this->name()).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(IoTools::IdToHuman::toTag(compId).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 4, pTmp);
      if(makeSquare)
      {
         Py_INCREF(Py_True);
         PyTuple_SetItem(pArgs, 5, Py_True);
      } else
      {
         Py_INCREF(Py_False);
         PyTuple_SetItem(pArgs, 5, Py_False);
      }

      // Call model operator Python routine
      PythonModelWrapper::setMethod(IoTools::IdToHuman::toString(ModelOperator::STENCIL));
      pValue = PythonModelWrapper::callMethod(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonModelWrapper::fillMatrix(mat, pValue);
      Py_DECREF(pValue);

      // Finalise Python interpreter
      PythonModelWrapper::cleanup();
   }

   void IEquation::dispatchExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const ModelOperator::Id opId,  const SpectralFieldId fieldId, const int matIdx, const SharedResolution spRes, const std::vector<MHDFloat>& eigs) const
   {
      // Get first four standard arguments in a tuple of size 7
      PyObject *pArgs = this->dispatchBaseArguments(7, ModelOperatorBoundary::FIELD_TO_RHS, spRes, eigs);

      // Prepare Python call arguments
      PyObject *pTmp, *pValue;

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

      // Set the restriction option
      #ifdef QUICC_MPI
         #ifdef QUICC_MPISPSOLVE
            std::vector<int> slow;
            std::vector<std::vector<int> > middle;

            spRes->buildRestriction(slow, middle, matIdx);
            PyObject *pSlow, *pMiddle;

            if(middle.size() > 0)
            {
               pSlow = PythonTools::makeList(slow);
               pMiddle = PythonTools::makeList(middle);
               pTmp = PyTuple_New(2);
               PyTuple_SetItem(pTmp, 0, pSlow);
               PyTuple_SetItem(pTmp, 1, pMiddle);
            } else
            {
               pTmp = PythonTools::makeList(slow);
            }

            PyTuple_SetItem(pArgs, 6, pTmp);
         #else
            // MPI code with serial sparse solver
            Py_INCREF(Py_None);
            PyTuple_SetItem(pArgs, 6, Py_None);
         #endif //QUICC_MPISPSOLVE
      #else
         // Serial code can't have any restriction
         Py_INCREF(Py_None);
         PyTuple_SetItem(pArgs, 6, Py_None);
      #endif //QUICC_MPI

      // Call model operator Python routine
      PythonModelWrapper::setMethod(IoTools::IdToHuman::toString(opId));
      pValue = PythonModelWrapper::callMethod(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      PythonModelWrapper::fillMatrix(mat, pValue);
      Py_DECREF(pValue);

      // Finalise Python interpreter
      PythonModelWrapper::cleanup();
   }

   void IEquation::setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const
   {
      // This implementation should never get called!
      throw Exception("Called dummy implementation of setGalerkinStencil!");
   }

   void IEquation::setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const ModelOperator::Id opId, const SpectralFieldId fieldId, const int matIdx) const
   {
      // This implementation should never get called!
      throw Exception("Called dummy implementation of setExplicitBlock!");
   }

   void IEquation::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // This implementation should never get called!
      throw Exception("Activated nonlinear term without implementation!");
   }

   void IEquation::useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id)
   {
      // DO NOTHING IN GENERAL!
   }

   Datatypes::SpectralScalarType::PointType IEquation::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
      // This implementation should never get called!
      throw Exception("Activated source term without implementation!");

      return Datatypes::SpectralScalarType::PointType();
   }

   void  IEquation::buildModelMatrix(DecoupledZSparse& rModelMatrix, const ModelOperator::Id opId, FieldComponents::Spectral::Id compId, const int matIdx, const ModelOperatorBoundary::Id bcType) const
   {
      // This implementation should never get called!
      throw Exception("Called dummy implementation of buildModelMatrix!");
   }
}
}
