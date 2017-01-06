/** 
 * @file CylinderTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII energy calculation for a toroidal/poloidal field in a cylinder (Worland expansion)
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
#include "IoVariable/CylinderTorPolEnergyWriter.hpp"

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

   CylinderTorPolEnergyWriter::CylinderTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   CylinderTorPolEnergyWriter::~CylinderTorPolEnergyWriter()
   {
   }

   void CylinderTorPolEnergyWriter::init()
   {
      // Spherical cylindrical volume: 
      this->mVolume = Math::PI;

      // Python variables
      PyObject *pArgs, *pValue;
      SparseMatrix tmpMat;

      // Initialise python wrapper
      PythonWrapper::init();

      // Load module for R direction (Worland)
      PythonWrapper::import("quicc.geometry.cylindrical.cylinder_radius_worland");

      // Prepare arguments
      pArgs = PyTuple_New(2);

      // Get resolution
      int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();
      pValue = PyLong_FromLong(cols);
      PyTuple_SetItem(pArgs, 0, pValue);

      // Loop over m
      this->mRIntgOp.resize(cols,this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
      for(int i = 0; i < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); i++)
      {
         // Set m
         pValue = PyFloat_FromDouble(static_cast<MHDFloat>(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i)));
         PyTuple_SetItem(pArgs, 1, pValue);
         PythonWrapper::setFunction("integral");
         pValue = PythonWrapper::callFunction(pArgs);

         // Fill matrix and cleanup
         PythonWrapper::fillMatrix(tmpMat, pValue);
         this->mRIntgOp.col(i) = tmpMat.row(0).transpose();
         Py_DECREF(pValue);
      }

      // Cleanup wrapper and load module for Z direction (Chebyshev)
      PythonWrapper::cleanup();
      PythonWrapper::import("quicc.geometry.cartesian.cartesian_1d");

      // Prepare arguments
      pArgs = PyTuple_New(2);

      // Get resolution
      cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>();
      pValue = PyLong_FromLong(cols);
      PyTuple_SetItem(pArgs, 0, pValue);

      // Set scaling factor
      pValue = PyFloat_FromDouble(1.0);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Call integral
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      PythonWrapper::fillMatrix(tmpMat, pValue);
      this->mZIntgOp.resize(cols, 1);
      //this->mZIntgOp = tmpMat;
      Py_DECREF(pValue);

      // Finalize wrapper
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void CylinderTorPolEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

      this->mTorEnergy = vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).point(1,1,1).real();
      this->mPolEnergy = vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).point(1,1,1).real();

//      // Dealias toroidal variable data
//      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));
//      
//      // Recover dealiased BWD data
//      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();
//
//      // Get FWD storage
//      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().project(rOutVarTor.rData(), rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);
//
//      // Compute |f|^2
//      rOutVarTor.rData() = rOutVarTor.rData().array()*rOutVarTor.rData().conjugate().array();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().integrate_energy(rInVarTor.rData(), rOutVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);
//
//      // Compute integral over Chebyshev expansion and sum harmonics
//      this->mTorEnergy = 0.0;
//      this->mPolEnergy = 0.0;
//
//      MHDFloat lfactor = 0.0;
//      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//         double factor = 1.0;
//         // Loop over harmonic order m
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
//            // m = 0, no factor of two
//            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
//            {
//               factor = 1.0;
//            } else
//            { 
//               factor = 2.0;
//            }
//
//            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
//            {
//               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
//               lfactor = lfactor*(lfactor+1.0);
//
//               this->mTorEnergy += factor*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(j).real()).sum();
//            }
//         }
//      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
//         // Loop over harmonic degree l
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
//            lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
//            lfactor = lfactor*(lfactor+1.0);
//            int start = 0;
//            // m = 0, no factor of two
//            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
//            {
//               this->mTorEnergy += lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(0).real())(0);
//               start = 1;
//            }
//            this->mTorEnergy += 2.0*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).rightCols(rInVarTor.slice(k).cols()-start).real()).sum();
//         }
//      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
//
//      // Free BWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);
//
//      // Free FWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTor);
//
//      // Normalize by sphere volume: 4/3*pi
//      this->mTorEnergy /= 2*this->mVolume;
//
//      // Dealias poloidal variable data for Q component
//      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));
//
//      // Recover dealiased BWD data
//      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();
//
//      // Get FWD storage
//      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().project(rOutVarPolQ.rData(), rInVarPolQ.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);
//
//      // Compute |f|^2
//      rOutVarPolQ.rData() = rOutVarPolQ.rData().array()*rOutVarPolQ.rData().conjugate().array();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().integrate_energy(rInVarPolQ.rData(), rOutVarPolQ.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);
//
//      // Compute energy in Q component of QST decomposition
//      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//         // Loop over harmonic order m
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
//            // m = 0, no factor of two
//            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
//            {
//               factor = 1.0;
//            } else
//            { 
//               factor = 2.0;
//            }
//
//            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
//            {
//               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
//               lfactor = std::pow(lfactor*(lfactor+1.0),2);
//
//               this->mPolEnergy += factor*lfactor*(this->mIntgOp*rInVarPolQ.slice(k).col(j).real()).sum();
//            }
//         }
//      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
//         lfactor = 0.0;
//         // Loop over harmonic degree l
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
//            lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
//            lfactor = std::pow(lfactor*(lfactor+1.0),2);
//            int start = 0;
//            // m = 0, no factor of two
//            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
//            {
//               this->mPolEnergy += lfactor*(this->mIntgOp*rInVarPolQ.slice(k).col(0).real())(0);
//               start = 1;
//            }
//            this->mPolEnergy += 2.0*lfactor*(this->mIntgOp*rInVarPolQ.slice(k).rightCols(rInVarPolQ.slice(k).cols()-start).real()).sum();
//         }
//      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
//
//      // Free BWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolQ);
//
//      // Free FWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolQ);
//
//      // Dealias poloidal variable data for S component
//      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));
//
//      // Recover dealiased BWD data
//      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolS = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();
//
//      // Get FWD storage
//      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolS = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().project(rOutVarPolS.rData(), rInVarPolS.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::DIFFR);
//
//      // Compute |f|^2
//      rOutVarPolS.rData() = rOutVarPolS.rData().array()*rOutVarPolS.rData().conjugate().array();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().integrate_energy(rInVarPolS.rData(), rOutVarPolS.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);
//
//      // Compute energy in S component of QST decomposition
//      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//         // Loop over harmonic order m
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
//            // m = 0, no factor of two
//            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) == 0)
//            {
//               factor = 1.0;
//            } else
//            { 
//               factor = 2.0;
//            }
//
//            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
//            {
//               lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
//               lfactor = lfactor*(lfactor+1.0);
//
//               this->mPolEnergy += factor*lfactor*(this->mIntgOp*rInVarPolS.slice(k).col(j).real()).sum();
//            }
//         }
//      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFM
//         lfactor = 0.0;
//         // Loop over harmonic degree l
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
//            lfactor = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
//            lfactor = lfactor*(lfactor+1.0);
//            int start = 0;
//            // m = 0, no factor of two
//            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
//            {
//               this->mPolEnergy += lfactor*(this->mIntgOp*rInVarPolS.slice(k).col(0).real())(0);
//               start = 1;
//            }
//            this->mPolEnergy += 2.0*lfactor*(this->mIntgOp*rInVarPolS.slice(k).rightCols(rInVarPolS.slice(k).cols()-start).real()).sum();
//         }
//      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
//
//      // Free BWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolS);
//
//      // Free FWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolS);
//
//      // Normalize by the volume
//      this->mPolEnergy /= 2*this->mVolume;
   }

   void CylinderTorPolEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(2);

         energy(0) = this->mTorEnergy;
         energy(1) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy = energy(0);
         this->mPolEnergy = energy(1);
      #endif //QUICC_MPI

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
         FrameworkMacro::abort(99);

         throw Exception("Toroidal/Poloidal energy is NaN!");
      }
   }

}
}
