/** 
 * @file SphericalTorPolEnergyWriter.cpp
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
#include "IoVariable/SphericalTorPolEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   SphericalTorPolEnergyWriter::SphericalTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   SphericalTorPolEnergyWriter::~SphericalTorPolEnergyWriter()
   {
   }

   void SphericalTorPolEnergyWriter::init()
   {
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
      pValue = PyLong_FromLong(this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL) + 2);
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

      this->mIntgOp = tmpAvg*tmpR2.leftCols(tmpR2.rows()-2);

      IVariableAsciiEWriter::init();
   }

   void SphericalTorPolEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

      // Dealias variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rComp(FieldComponents::SpectralTor).rTotal());
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVar.rData(), rInVar.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ, Arithmetics::SET);

      // Compute |f|^2
      rOutVar.rData() = rOutVar.rData().array()*rOutVar.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate(rInVar.rData(), rOutVar.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

      // Compute integral over Chebyshev expansion and sum harmonics
      this->mEnergy = 0.0;

      #ifdef GEOMHDISCC_SPATIALSCHEME_SLFM
         int start = 0;
         if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0) == 0)
         {
            this->mEnergy += (this->mIntgOp*rInVar.slice(0).topRows(this->mIntgOp.cols()).real()).sum();
            start = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(0);
         }

         this->mEnergy += 2.0*(this->mIntgOp*rInVar.data().topRightCorner(this->mIntgOp.cols(),rInVar.data().cols()-start).real()).sum();
      #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFM
      #ifdef GEOMHDISCC_SPATIALSCHEME_SLFL
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               this->mEnergy += (this->mIntgOp*rInVar.slice(k).topLeftCorner(this->mIntgOp.cols(),1).real())(0);
               start = 1;
            }
            this->mEnergy += 2.0*(this->mIntgOp*rInVar.slice(k).topRightCorner(this->mIntgOp.cols(),rInVar.slice(k).cols()-start).real()).sum();
         }
      #endif //GEOMHDISCC_SPATIALSCHEME_SLFL

      // Normalize by sphere volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
      MHDFloat ri = ro*this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;
      this->mEnergy /= (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rInVar);
   }

   void SphericalTorPolEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef GEOMHDISCC_MPI
         Array energy(2);

         energy(0) = this->mTorEnergy;
         energy(1) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy = energy(0);
         this->mPolEnergy = energy(1);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << this->mTorEnergy + this->mPolEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
