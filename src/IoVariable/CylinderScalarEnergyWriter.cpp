/** 
 * @file CylinderScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII energy calculation for scalar field in a cylinder (Worland expansion)
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
#include "IoVariable/CylinderScalarEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace QuICC {

namespace IoVariable {

   CylinderScalarEnergyWriter::CylinderScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL, IVariableAsciiWriter::EXTEND), mEnergy(-1.0)
   {
   }

   CylinderScalarEnergyWriter::~CylinderScalarEnergyWriter()
   {
   }

   void CylinderScalarEnergyWriter::init()
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

      IVariableAsciiWriter::init();
   }

   void CylinderScalarEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      this->mEnergy = sRange.first->second->dom(0).total().point(1,1,1).real();

//      // Dealias variable data
//      coord.communicator().dealiasSpectral(sRange.first->second->rDom(0).rTotal());
//      
//      // Recover dealiased BWD data
//      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();
//
//      // Get FWD storage
//      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().project(rOutVar.rData(), rInVar.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);
//
//      // Compute |f|^2
//      rOutVar.rData() = rOutVar.rData().array()*rOutVar.rData().conjugate().array();
//
//      // Compute projection transform for first dimension 
//      coord.transform1D().integrate_energy(rInVar.rData(), rOutVar.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);
//
//      // Compute integral over Chebyshev expansion and sum harmonics
//      this->mEnergy = 0.0;
//
//      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//         double factor = 1.0;
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
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
//               this->mEnergy += factor*(this->mSphIntgOpEven*rInVar.slice(k).col(j).real()).sum();
//            }
//         }
//      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
//      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
//         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//         {
//            int start = 0;
//            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
//            {
//               this->mEnergy += (this->mSphIntgOpEven*rInVar.slice(k).col(0).real())(0);
//               start = 1;
//            }
//            this->mEnergy += 2.0*(this->mSphIntgOpEven*rInVar.slice(k).rightCols(rInVar.slice(k).cols()-start).real()).sum();
//         }
//      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
//
//      // Free BWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);
//
//      // Free FWD storage
//      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVar);

      // Normalize by the spherical volume
      this->mEnergy /= this->mVolume;
   }

   void CylinderScalarEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

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
         FrameworkMacro::abort(99);

         throw Exception("Scalar energy is NaN!");
      }
   }

}
}
