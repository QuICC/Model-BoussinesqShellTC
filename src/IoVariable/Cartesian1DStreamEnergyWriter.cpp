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
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   Cartesian1DStreamEnergyWriter::Cartesian1DStreamEnergyWriter(const std::string& prefix, const std::string& type, const bool hasZonalX, const bool hasZonalY)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mHasZonalX(hasZonalX), mHasZonalY(hasZonalY), mXEnergy(-1.0), mYEnergy(-1.0), mZEnergy(-1.0), mXZonalXEnergy(-1.0), mZZonalXEnergy(-1.0), mYZonalYEnergy(-1.0), mZZonalYEnergy(-1.0)
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
      // Compute total energy
      this->computeEnergy(coord, TOTAL);

      // Compute zonal energy with respect to X axis
      if(this->mHasZonalX)
      {
         this->computeEnergy(coord, ZONAL_X);
      }

      // Compute zonal energy with respect to Y axis
      if(this->mHasZonalY)
      {
         this->computeEnergy(coord, ZONAL_Y);
      }
   }

   void Cartesian1DStreamEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef GEOMHDISCC_MPI
         int n = 3;
         if(this->mHasZonalX)
         {
            n += 2;
         }
         if(this->mHasZonalY)
         {
            n += 2;
         }
         Array energy(n);

         energy(0) = this->mXEnergy;
         energy(1) = this->mYEnergy;
         energy(2) = this->mZEnergy;
         int i = 3;
         if(this->mHasZonalX)
         {
            energy(i) = this->mXZonalXEnergy;
            energy(i+1) = this->mZZonalXEnergy;
            i += 2;
         }
         if(this->mHasZonalY)
         {
            energy(i) = this->mYZonalYEnergy;
            energy(i+1) = this->mZZonalYEnergy;
         }

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mXEnergy = energy(0);
         this->mYEnergy = energy(1);
         this->mZEnergy = energy(2);
         i = 3;
         if(this->mHasZonalX)
         {
            this->mXZonalXEnergy = energy(i);
            this->mZZonalXEnergy = energy(i+1);
            i += 2;
         }
         if(this->mHasZonalY)
         {
            this->mYZonalYEnergy = energy(i);
            this->mZZonalYEnergy = energy(i+1);
         }
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mXEnergy + this->mYEnergy + this->mZEnergy << "\t" << this->mXEnergy << "\t" << this->mYEnergy << "\t" << this->mZEnergy << std::endl;

         if(this->mHasZonalX)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mXZonalXEnergy + this->mZZonalXEnergy << "\t" << this->mXZonalXEnergy << "\t" << this->mZZonalXEnergy;
         }

         if(this->mHasZonalY)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mYZonalYEnergy + this->mZZonalYEnergy << "\t" << this->mYZonalYEnergy << "\t" << this->mZZonalYEnergy;
         }

         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mXEnergy) || std::isnan(this->mYEnergy) || std::isnan(this->mZEnergy))
      {
         #ifdef GEOMHDISCC_MPI
            MPI_Abort(MPI_COMM_WORLD, 99);
         #endif //GEOMHDISCC_MPI

         throw Exception("Kinetic energy is NaN!");
      }
   }

   void Cartesian1DStreamEnergyWriter::computeEnergy(Transform::TransformCoordinatorType& coord, const Cartesian1DStreamEnergyWriter::EnergyTypeId flag)
   {
      scalar_iterator sIt;
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 2);

      MHDFloat *pX;
      MHDFloat *pY;
      MHDFloat *pZ;
      MHDFloat tmp;
      if(flag == TOTAL)
      {
         pX = &this->mXEnergy;
         pY = &this->mYEnergy;
         pZ = &this->mZEnergy;

      } else if(flag == ZONAL_X)
      {
         pX = &this->mXZonalXEnergy;
         pY = &tmp;
         pZ = &this->mZZonalXEnergy;
      } else if(flag == ZONAL_Y)
      {
         pX = &tmp;
         pY = &this->mYZonalYEnergy;
         pZ = &this->mZZonalYEnergy;
      }

      // Initialize the energy
      *pX = 0.0;
      *pY = 0.0;
      *pZ = 0.0;

      std::map<FieldComponents::Spectral::Id, scalar_iterator> comps;
      sIt = sRange.first;
      sIt++;
      if(sRange.first->first == PhysicalNames::VELOCITYZ)
      {
         if(flag == TOTAL || flag == ZONAL_X)
         {
            comps.insert(std::make_pair(FieldComponents::Spectral::X, sIt));
         }
         if(flag == TOTAL || flag == ZONAL_Y)
         {
            comps.insert(std::make_pair(FieldComponents::Spectral::Y, sIt));
         }
         comps.insert(std::make_pair(FieldComponents::Spectral::Z, sRange.first));
      } else
      {
         if(flag == TOTAL || flag == ZONAL_X)
         {
            comps.insert(std::make_pair(FieldComponents::Spectral::X, sRange.first));
         }
         if(flag == TOTAL || flag == ZONAL_Y)
         {
            comps.insert(std::make_pair(FieldComponents::Spectral::Y, sRange.first));
         }
         comps.insert(std::make_pair(FieldComponents::Spectral::Z, sIt));
      }

      for(std::map<FieldComponents::Spectral::Id,scalar_iterator>::iterator cIt = comps.begin(); cIt != comps.end(); ++cIt)
      {
         // Dealias variable data
         coord.communicator().dealiasSpectral(cIt->second->second->rDom(0).rTotal());

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Restrict field if necessary
         if(flag == ZONAL_X)
         {
            int nK = this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
            for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) != 0)
               {
                  rInVar.setSlice(MatrixZ::Zero(rInVar.slice(k).rows(), rInVar.slice(k).cols()),k);
               }
            }
         } else if(flag == ZONAL_Y)
         {
            for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
               {
                  if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k) != 0)
                  {
                     rInVar.setProfile(ArrayZ::Zero(rInVar.profile(j,k).rows()),j,k);
                  }
               }
            }
         }

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
            int nK = this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
            for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               int k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
               // Convert to double Fourier mode
               if(k_ >= nK/2 + (nK % 2))
               {
                  k_ = k_ - nK;
               }

               rOutVar.setSlice(std::pow(k_,2)*(rOutVar.slice(k).array()*rOutVar.slice(k).conjugate().array()).matrix(),k);
            }
         }

         // Compute integration transform for first dimension 
         coord.transform1D().integrate_full(rInVar.rData(), rOutVar.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

         MHDFloat *pEnergy;
         if(cIt->first == FieldComponents::Spectral::X)
         {
            pEnergy = pX;
         } else if(cIt->first == FieldComponents::Spectral::Y)
         {
            pEnergy = pY;
         } else
         {
            pEnergy = pZ;
         }

         // Compute integral over Chebyshev expansion and sum Fourier coefficients
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               // Include ignored complex conjugate
               if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
               {
                  *pEnergy += (this->mIntgOp*rInVar.slice(k).col(0).real())(0);
               } else if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
               {
                  *pEnergy += 2.0*(this->mIntgOp*rInVar.slice(k).col(0).real())(0);
               }
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

}
}
