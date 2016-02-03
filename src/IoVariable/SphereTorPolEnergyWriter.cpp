/** 
 * @file SphereTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a sphere
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
#include "IoVariable/SphereTorPolEnergyWriter.hpp"

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

   SphereTorPolEnergyWriter::SphereTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   SphereTorPolEnergyWriter::~SphereTorPolEnergyWriter()
   {
   }

   void SphereTorPolEnergyWriter::init()
   {
      // Spherical shell volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.spherical.sphere_radius");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(3);

      // ... create boundray condition (none)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get resolution
      int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>() + 2;
      pValue = PyLong_FromLong(cols);
      PyTuple_SetItem(pArgs, 0, pValue);

      // Set even basis
      pValue = PyLong_FromLong(0);
      PyTuple_SetItem(pArgs, 1, pValue);
      // Call r^2
      PythonWrapper::setFunction("r2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpR2;
      PythonWrapper::fillMatrix(tmpR2, pValue);
      Py_DECREF(pValue);

      PyObject *pTmp = PyTuple_New(2);
      pTmp = PyTuple_GetSlice(pArgs, 0, 2);
      // Call avg
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pTmp);
      // Fill matrix and cleanup
      SparseMatrix tmpAvg;
      PythonWrapper::fillMatrix(tmpAvg, pValue);
      Py_DECREF(pValue);

      // Store integral
      this->mIntgOpEven = tmpAvg.leftCols(cols-2);

      // Store spherical integral (include r^2 factor)
      this->mSphIntgOpEven = tmpAvg*tmpR2.leftCols(cols-2);

      // Set odd basis
      pValue = PyLong_FromLong(1);
      PyTuple_SetItem(pArgs, 1, pValue);
      // Call r^2
      PythonWrapper::setFunction("r2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      PythonWrapper::fillMatrix(tmpR2, pValue);
      Py_DECREF(pValue);

      pTmp = PyTuple_GetSlice(pArgs, 0, 2);
      // Call avg
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pTmp);
      // Fill matrix and cleanup
      PythonWrapper::fillMatrix(tmpAvg, pValue);
      Py_DECREF(pValue);

      // Store integral
      this->mIntgOpOdd = tmpAvg.leftCols(cols-2);

      // Store spherical integral (include r^2 factor)
      this->mSphIntgOpOdd = tmpAvg*tmpR2.leftCols(cols-2);
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void SphereTorPolEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

      // Dealias toroidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVarTor.rData(), rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ, Arithmetics::SET);

      // Compute |f|^2
      rOutVarTor.rData() = rOutVarTor.rData().array()*rOutVarTor.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarTor.rData(), rOutVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

      // Compute integral over Chebyshev expansion and sum harmonics
      this->mTorEnergy = 0.0;
      this->mPolEnergy = 0.0;

      SparseMatrix *op;

      #ifdef GEOMHDISCC_SPATIALSCHEME_BLFM
         MHDFloat lfactor = 0.0;
         int start = 0;
         int m0 = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         if(m0 == 0)
         {
            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(0); ++j)
            {
               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, 0);
               lfactor = lfactor*(lfactor+1.0);

               this->mTorEnergy += lfactor*(this->mSphIntgOp*rInVarTor.slice(0).col(j).real()).sum();
            }
            start = 1;
         }
         for(int k = start; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = lfactor*(lfactor+1.0);

               this->mTorEnergy += 2.0*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(j).real()).sum();
            }
         }
      #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFM
      #ifdef GEOMHDISCC_SPATIALSCHEME_BLFL
         MHDFloat lfactor = 0.0;
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) % 2 == 0)
            {
               op = &this->mSphIntgOpEven;
            } else
            {
               op = &this->mSphIntgOpOdd;
            }
            lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = lfactor*(lfactor+1.0);
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               this->mTorEnergy += lfactor*((*op)*rInVarTor.slice(k).col(0).real())(0);
               start = 1;
            }
            this->mTorEnergy += 2.0*lfactor*((*op)*rInVarTor.slice(k).rightCols(rInVarTor.slice(k).cols()-start).real()).sum();
         }
      #endif //GEOMHDISCC_SPATIALSCHEME_BLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTor);

      // Normalize by sphere volume: 4/3*pi
      this->mTorEnergy /= 2*this->mVolume;

      // Dealias poloidal variable data for Q component
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVarPolQ.rData(), rInVarPolQ.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ, Arithmetics::SET);

      // Compute |f|^2
      rOutVarPolQ.rData() = rOutVarPolQ.rData().array()*rOutVarPolQ.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarPolQ.rData(), rOutVarPolQ.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

      // Compute energy in Q component of QST decomposition
      #ifdef GEOMHDISCC_SPATIALSCHEME_BLFM
         start = 0;
         m0 = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         if(m0 == 0)
         {
            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(0); ++j)
            {
               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, 0);
               lfactor = std::pow(lfactor*(lfactor+1.0),2);

               this->mPolEnergy += lfactor*(this->mIntgOp*rInVarPolQ.slice(0).col(j).real()).sum();
            }
            start = 1;
         }
         for(int k = start; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = std::pow(lfactor*(lfactor+1.0),2);

               this->mPolEnergy += 2.0*lfactor*(this->mIntgOp*rInVarPolQ.slice(k).col(j).real()).sum();
            }
         }
      #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFM
      #ifdef GEOMHDISCC_SPATIALSCHEME_BLFL
         lfactor = 0.0;
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) % 2 == 0)
            {
               op = &this->mIntgOpEven;
            } else
            {
               op = &this->mIntgOpOdd;
            }
            lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = std::pow(lfactor*(lfactor+1.0),2);
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               this->mPolEnergy += lfactor*((*op)*rInVarPolQ.slice(k).col(0).real())(0);
               start = 1;
            }
            this->mPolEnergy += 2.0*lfactor*((*op)*rInVarPolQ.slice(k).rightCols(rInVarPolQ.slice(k).cols()-start).real()).sum();
         }
      #endif //GEOMHDISCC_SPATIALSCHEME_BLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolQ);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolQ);

      // Dealias poloidal variable data for S component
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolS = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolS = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVarPolS.rData(), rInVarPolS.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::DIFFR, Arithmetics::SET);

      // Compute |f|^2
      rOutVarPolS.rData() = rOutVarPolS.rData().array()*rOutVarPolS.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarPolS.rData(), rOutVarPolS.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG, Arithmetics::SET);

      // Compute energy in Q component of QST decomposition
      #ifdef GEOMHDISCC_SPATIALSCHEME_BLFM
         start = 0;
         m0 = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         if(m0 == 0)
         {
            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(0); ++j)
            {
               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, 0);
               lfactor = lfactor*(lfactor+1.0);

               this->mPolEnergy += lfactor*(this->mIntgOp*rInVarPolS.slice(0).col(j).real()).sum();
            }
            start = 1;
         }
         for(int k = start; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
            {
               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = lfactor*(lfactor+1.0);

               this->mPolEnergy += 2.0*lfactor*(this->mIntgOp*rInVarPolS.slice(k).col(j).real()).sum();
            }
         }
      #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFM
      #ifdef GEOMHDISCC_SPATIALSCHEME_BLFL
         lfactor = 0.0;
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) % 2 == 0)
            {
               op = &this->mIntgOpEven;
            } else
            {
               op = &this->mIntgOpOdd;
            }
            lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = lfactor*(lfactor+1.0);
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               this->mPolEnergy += lfactor*((*op)*rInVarPolS.slice(k).col(0).real())(0);
               start = 1;
            }
            this->mPolEnergy += 2.0*lfactor*((*op)*rInVarPolS.slice(k).rightCols(rInVarPolS.slice(k).cols()-start).real()).sum();
         }
      #endif //GEOMHDISCC_SPATIALSCHEME_BLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolS);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolS);

      // Normalize by the volume
      this->mPolEnergy /= 2*this->mVolume;
   }

   void SphereTorPolEnergyWriter::write()
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
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorEnergy + this->mPolEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mTorEnergy) || std::isnan(this->mPolEnergy))
      {
         #ifdef GEOMHDISCC_MPI
            MPI_Abort(MPI_COMM_WORLD, 99);
         #endif //GEOMHDISCC_MPI

         throw Exception("Toroidal/Poloidal energy is NaN!");
      }
   }

}
}
