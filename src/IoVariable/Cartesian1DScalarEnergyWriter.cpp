/** 
 * @file Cartesian1DScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) energy calculation for scalar field
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
#include "IoVariable/Cartesian1DScalarEnergyWriter.hpp"

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

   Cartesian1DScalarEnergyWriter::Cartesian1DScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mEnergy(-1.0)
   {
   }

   Cartesian1DScalarEnergyWriter::~Cartesian1DScalarEnergyWriter()
   {
   }

   void Cartesian1DScalarEnergyWriter::init()
   {
      // Normalize by Cartesian volume V = (2*pi)*(2*pi)*2 but FFT already includes 1/(2*pi)
      this->mVolume = 2.0;

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.cartesian.cartesian_1d");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(2);

      // Get resolution
      int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();
      pValue = PyLong_FromLong(cols);
      PyTuple_SetItem(pArgs, 0, pValue);

      // .. set scale factor
      pValue = PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::SCALE1D))->second);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Call integral
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      PythonWrapper::fillMatrix(this->mIntgOp, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void Cartesian1DScalarEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      // Initialize the energy
      this->mEnergy = 0.0;

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

      // Compute integration transform for first dimension 
      coord.transform1D().integrate_full(rInVar.rData(), rOutVar.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         int start = 0;
         if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
         {
            this->mEnergy += (this->mIntgOp*rInVar.slice(k).col(0).real())(0);
            start = 1;
         }
         this->mEnergy += 2.0*(this->mIntgOp*rInVar.slice(k).rightCols(rInVar.slice(k).cols()-start).real()).sum();
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVar);

      // Normalize by the Cartesian volume
      this->mEnergy /= this->mVolume;
   }

   void Cartesian1DScalarEnergyWriter::write()
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
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mEnergy))
      {
         #ifdef GEOMHDISCC_MPI
            MPI_Abort(MPI_COMM_WORLD, 99);
         #endif //GEOMHDISCC_MPI

         throw Exception("Kinetic energy is NaN!");
      }
   }

}
}
