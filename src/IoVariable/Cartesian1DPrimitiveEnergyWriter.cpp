/** 
 * @file Cartesian1DPrimitiveEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (primitive formulation)
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
#include "IoVariable/Cartesian1DPrimitiveEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   Cartesian1DPrimitiveEnergyWriter::Cartesian1DPrimitiveEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mXEnergy(-1.0), mYEnergy(-1.0), mZEnergy(-1.0)
   {
   }

   Cartesian1DPrimitiveEnergyWriter::~Cartesian1DPrimitiveEnergyWriter()
   {
   }

   void Cartesian1DPrimitiveEnergyWriter::init()
   {
      // Normalize by Cartesian volume: (2*pi)*(2*pi)*(2*
      this->mVolume = std::pow(2.0*Math::PI,2)*2;

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

      // Call avg
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpAvg(cols,cols);
      PythonWrapper::fillMatrix(this->mIntgOp, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void Cartesian1DPrimitiveEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::X);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::Y);
      assert(FieldComponents::Spectral::THREE == FieldComponents::Spectral::Z);

      // Initialize the energy
      this->mXEnergy = 0.0;
      this->mYEnergy = 0.0;
      this->mZEnergy = 0.0;

      std::vector<FieldComponents::Spectral::Id> comps;
      comps.push_back(FieldComponents::Spectral::X);
      comps.push_back(FieldComponents::Spectral::Y);
      comps.push_back(FieldComponents::Spectral::Z);

      for(std::vector<FieldComponents::Spectral::Id>::iterator cIt = comps.begin(); cIt != comps.end(); ++cIt)
      {
         // Dealias toroidal variable data
         coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(*cIt));

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension 
         coord.transform1D().project(rOutVarTor.rData(), rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ, Arithmetics::SET);

         // Compute |f|^2
         rOutVarTor.rData() = rOutVarTor.rData().array()*rOutVarTor.rData().conjugate().array();

         // Compute integration transform for first dimension 
         coord.transform1D().integrate_full(rInVarTor.rData(), rOutVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

         if(*cIt == FieldComponents::Spectral::X)
         {
            // Compute integral over Chebyshev expansion and sum Fourier coefficients
            this->mXEnergy += (this->mIntgOp*rInVar.data().rightCols(rInVar.data().cols()).real()).sum();
            // Normalize by box volume
            this->mXEnergy /= this->mVolume;

         } else if(*cIt == FieldComponents::Spectral::X)
         {
            // Compute integral over Chebyshev expansion and sum Fourier coefficients
            this->mYEnergy += (this->mIntgOp*rInVar.data().rightCols(rInVar.data().cols()).real()).sum();
            // Normalize by box volume
            this->mYEnergy /= this->mVolume;

         } else
         {
            // Compute integral over Chebyshev expansion and sum Fourier coefficients
            this->mZEnergy += (this->mIntgOp*rInVar.data().rightCols(rInVar.data().cols()).real()).sum();
            // Normalize by box volume
            this->mZEnergy /= this->mVolume;
         }

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTor);
      }
   }

   void Cartesian1DPrimitiveEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef GEOMHDISCC_MPI
         Array energy(3);

         energy(0) = this->mXEnergy;
         energy(1) = this->mYEnergy;
         energy(2) = this->mZEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mXEnergy = energy(0);
         this->mYEnergy = energy(1);
         this->mZEnergy = energy(2);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mXEnergy + this->mYEnergy + this->mZEnergy << "\t" << this->mXEnergy << "\t" << this->mYEnergy << "\t" << this->mZEnergy << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
