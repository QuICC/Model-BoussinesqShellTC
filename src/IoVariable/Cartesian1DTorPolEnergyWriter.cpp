/** 
 * @file Cartesian1DTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (toroidal/poloidal formulation)
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
#include "IoVariable/Cartesian1DTorPolEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   Cartesian1DTorPolEnergyWriter::Cartesian1DTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mXEnergy(-1.0), mYEnergy(-1.0), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   Cartesian1DTorPolEnergyWriter::~Cartesian1DTorPolEnergyWriter()
   {
   }

   void Cartesian1DTorPolEnergyWriter::init()
   {
      // Normalize by Cartesian volume V = (2*pi*Box1D/k1D)*(2*pi*Box2D/k2D)*2 but FFT already includes 1/(2*pi)
      this->mVolume = 2.0/(this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM2D)*this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM3D));

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("quicc.geometry.cartesian.cartesian_1d");

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

   void Cartesian1DTorPolEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);

      // Initialize the energy
      this->mXEnergy = 0.0;
      this->mYEnergy = 0.0;
      this->mTorEnergy = 0.0;
      this->mPolEnergy = 0.0;

      // Dealias toroidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVarTor.rData(), rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Compute |f|^2
      rOutVarTor.rData() = rOutVarTor.rData().array()*rOutVarTor.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarTor.rData(), rOutVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      int nK = this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         MHDFloat k_ = static_cast<MHDFloat>(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k));
         // Convert to double Fourier mode
         if(k_ >= nK/2 + (nK % 2))
         {
            k_ = k_ - nK;
         }
         // Include boxscale
         k_ *= this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM3D);
         

         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            MHDFloat j_ = static_cast<MHDFloat>(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k));
            // Include boxscale
            j_ *= this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM2D);
            MHDFloat factor = (k_*k_ + j_*j_);

            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0)
            {
               // Include ignored complex conjugate
               if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
               {
                  this->mXEnergy += (this->mIntgOp*rInVarTor.slice(k).col(0).real())(0);
               } else if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  this->mTorEnergy += 2.0*factor*(this->mIntgOp*rInVarTor.slice(k).col(0).real())(0);
               }
            } else
            {
               this->mTorEnergy += 2.0*factor*(this->mIntgOp*rInVarTor.slice(k).col(j).real())(0);
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTor);

      // Dealias poloidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVarPol.rData(), rInVarPol.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Compute |f|^2
      rOutVarPol.rData() = rOutVarPol.rData().array()*rOutVarPol.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarPol.rData(), rOutVarPol.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         MHDFloat k_ = static_cast<MHDFloat>(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k));
         // Convert to double Fourier mode
         if(k_ >= nK/2 + (nK % 2))
         {
            k_ = k_ - nK;
         }
         // Include boxscale
         k_ *= this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM3D);

         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            MHDFloat j_ = static_cast<MHDFloat>(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k));
            // Include boxscale
            j_ *= this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM2D);
            MHDFloat factor = std::pow((k_*k_ + j_*j_),2);

            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0)
            {
               // Include ignored complex conjugate
               if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
               {
                  this->mYEnergy += (this->mIntgOp*rInVarPol.slice(k).col(0).real())(0);
               } else if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  this->mPolEnergy += 2.0*factor*(this->mIntgOp*rInVarPol.slice(k).col(0).real())(0);
               }
            } else
            {
               this->mPolEnergy += 2.0*factor*(this->mIntgOp*rInVarPol.slice(k).col(j).real())(0);
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPol);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPol);

      // Dealias poloidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolDz = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolDz = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVarPolDz.rData(), rInVarPolDz.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::DIFF);

      // Compute |f|^2
      rOutVarPolDz.rData() = rOutVarPolDz.rData().array()*rOutVarPolDz.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarPolDz.rData(), rOutVarPolDz.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Compute integral over Chebyshev expansion and sum Fourier coefficients
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         MHDFloat k_ = static_cast<MHDFloat>(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k));
         // Convert to double Fourier mode
         if(k_ >= nK/2 + (nK % 2))
         {
            k_ = k_ - nK;
         }
         // Include boxscale
         k_ *= this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM3D);

         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            MHDFloat j_ = static_cast<MHDFloat>(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k));
            // Include boxscale
            j_ *= this->mspRes->sim()->boxScale(Dimensions::Simulation::SIM2D);
            MHDFloat factor = (k_*k_ + j_*j_);

            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) == 0)
            {
               // Include ignored complex conjugate
               if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
               {
                  // Do nothing
               } else if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  this->mPolEnergy += 2.0*factor*(this->mIntgOp*rInVarPolDz.slice(k).col(0).real())(0);
               }
            } else
            {
               this->mPolEnergy += 2.0*factor*(this->mIntgOp*rInVarPolDz.slice(k).col(j).real())(0);
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolDz);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolDz);

      // Normalize by the Cartesian volume
      this->mXEnergy /= 2.0*this->mVolume;
      this->mYEnergy /= 2.0*this->mVolume;
      this->mTorEnergy /= 2.0*this->mVolume;
      this->mPolEnergy /= 2.0*this->mVolume;
   }

   void Cartesian1DTorPolEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(4);

         energy(0) = this->mXEnergy;
         energy(1) = this->mYEnergy;
         energy(2) = this->mTorEnergy;
         energy(3) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), energy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mXEnergy = energy(0);
         this->mYEnergy = energy(1);
         this->mTorEnergy = energy(2);
         this->mPolEnergy = energy(3);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mXEnergy + this->mYEnergy + this->mTorEnergy + this->mPolEnergy << "\t" << this->mXEnergy << "\t" << this->mYEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mXEnergy) || std::isnan(this->mYEnergy) || std::isnan(this->mTorEnergy) || std::isnan(this->mPolEnergy))
      {
         FrameworkMacro::abort(99);

         throw Exception("Kinetic energy is NaN!");
      }
   }

}
}
