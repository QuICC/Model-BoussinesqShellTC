/** 
 * @file Cartesian1DStreamEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (streamfunction formulation)
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
#include "IoVariable/Cartesian1DStreamEnergyWriter.hpp"

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

   Cartesian1DStreamEnergyWriter::Cartesian1DStreamEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mXEnergy(-1.0), mYEnergy(-1.0), mZEnergy(-1.0)
   {
   }

   Cartesian1DStreamEnergyWriter::~Cartesian1DStreamEnergyWriter()
   {
   }

   void Cartesian1DStreamEnergyWriter::init()
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

      // Call avg
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      PythonWrapper::fillMatrix(this->mIntgOp, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void Cartesian1DStreamEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator sIt;
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 2);

      // Initialize the energy
      this->mXEnergy = 0.0;
      this->mYEnergy = 0.0;
      this->mZEnergy = 0.0;

      std::map<FieldComponents::Spectral::Id, scalar_iterator> comps;
      sIt = sRange.first;
      sIt++;
      if(sRange.first->first == PhysicalNames::VELOCITYZ)
      {
         comps.insert(std::make_pair(FieldComponents::Spectral::X, sIt));
         comps.insert(std::make_pair(FieldComponents::Spectral::Y, sIt));
         comps.insert(std::make_pair(FieldComponents::Spectral::Z, sRange.first));
      } else
      {
         comps.insert(std::make_pair(FieldComponents::Spectral::X, sRange.first));
         comps.insert(std::make_pair(FieldComponents::Spectral::Y, sRange.first));
         comps.insert(std::make_pair(FieldComponents::Spectral::Z, sIt));
      }

      for(std::map<FieldComponents::Spectral::Id,scalar_iterator>::iterator cIt = comps.begin(); cIt != comps.end(); ++cIt)
      {
         // Dealias variable data
         coord.communicator().dealiasSpectral(cIt->second->second->rDom(0).rTotal());

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension 
         coord.transform1D().project(rOutVar.rData(), rInVar.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ, Arithmetics::SET);

         if(cIt->first == FieldComponents::Spectral::Z)
         {
            // Compute |f|^2
            rOutVar.rData() = rOutVar.rData().array()*rOutVar.rData().conjugate().array();
         } else if(cIt->first == FieldComponents::Spectral::X)
         {
            for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
               {
                  int j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);

                  rOutVar.setProfile(std::pow(j_,2)*(rOutVar.profile(j,k).array()*rOutVar.profile(j,k).conjugate().array()).matrix(),j,k);
               }
            }
         } else
         {
            for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               int k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

               rOutVar.setSlice(std::pow(k_,2)*(rOutVar.slice(k).array()*rOutVar.slice(k).conjugate().array()).matrix(),k);
            }
         }

         // Compute integration transform for first dimension 
         coord.transform1D().integrate_full(rInVar.rData(), rOutVar.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

         MHDFloat *pEnergy;
         if(cIt->first == FieldComponents::Spectral::X)
         {
            pEnergy = &this->mXEnergy;
         } else if(cIt->first == FieldComponents::Spectral::Y)
         {
            pEnergy = &this->mYEnergy;
         } else
         {
            pEnergy = &this->mZEnergy;
         }

         // Compute integral over Chebyshev expansion and sum Fourier coefficients
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               *pEnergy += (this->mIntgOp*rInVar.slice(k).col(0).real())(0);
               start = 1;
            }
            *pEnergy += 2.0*(this->mIntgOp*rInVar.slice(k).rightCols(rInVar.slice(k).cols()-start).real()).sum();
         }

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVar);

         // Normalize by the Cartesian volume
         *pEnergy /= this->mVolume;
      }
   }

   void Cartesian1DStreamEnergyWriter::write()
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
