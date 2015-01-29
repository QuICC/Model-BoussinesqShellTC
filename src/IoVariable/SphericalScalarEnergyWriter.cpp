/** 
 * @file SphericalScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "IoVariable/SphericalScalarEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   SphericalScalarEnergyWriter::SphericalScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mEnergy(-1.0)
   {
   }

   SphericalScalarEnergyWriter::~SphericalScalarEnergyWriter()
   {
   }

   void SphericalScalarEnergyWriter::init()
   {
      // Normalize by sphere volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
      MHDFloat ri = ro*this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;
      this->mVolume = (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.spherical.shell_radius");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(4);

      // ... compute a, b factors
      PyObject *pTmp = PyTuple_New(2);
      PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second));
      PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second));
      PythonWrapper::setFunction("linear_r2x");
      pValue = PythonWrapper::callFunction(pTmp);
      PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue, 0));
      PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue, 1));
      // ... create boundray condition (none)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyTuple_SetItem(pArgs, 3, pValue);

      // Get resolution
      int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>() + 2;
      pValue = PyLong_FromLong(cols);
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call x^2
      PythonWrapper::setFunction("x2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpR2;
      PythonWrapper::fillMatrix(tmpR2, pValue);
      Py_DECREF(pValue);

      pTmp = PyTuple_GetSlice(pArgs, 0, 3);
      // Call avg
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pTmp);
      // Fill matrix and cleanup
      SparseMatrix tmpAvg;
      PythonWrapper::fillMatrix(tmpAvg, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      this->mSphIntgOp = tmpAvg*tmpR2.leftCols(cols-2);

      IVariableAsciiEWriter::init();
   }

   void SphericalScalarEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      // Dealias variable data
      coord.communicator().dealiasSpectral(sRange.first->second->rDom(0).rTotal());
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVar.rData(), rInVar.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ, Arithmetics::SET);

      // Compute |f|^2
      rOutVar.rData() = rOutVar.rData().array()*rOutVar.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVar.rData(), rOutVar.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

      // Compute integral over Chebyshev expansion and sum harmonics
      this->mEnergy = 0.0;

      #ifdef GEOMHDISCC_SPATIALSCHEME_SLFM
         int start = 0;
         if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0) == 0)
         {
            this->mEnergy += (this->mSphIntgOp*rInVar.slice(0).real()).sum();
            start = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(0);
         }

         this->mEnergy += 2.0*(this->mSphIntgOp*rInVar.data().rightCols(rInVar.data().cols()-start).real()).sum();
      #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM
      #ifdef GEOMHDISCC_SPATIALSCHEME_SLFL
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               this->mEnergy += (this->mSphIntgOp*rInVar.slice(k).col(0).real())(0);
               start = 1;
            }
            this->mEnergy += 2.0*(this->mSphIntgOp*rInVar.slice(k).rightCols(rInVar.slice(k).cols()-start).real()).sum();
         }
      #endif //GEOMHDISCC_SPATIALSCHEME_SLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVar);

      // Normalize by the spherical volume
      this->mEnergy /= this->mVolume;
   }

   void SphericalScalarEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << this->mEnergy << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
