/** 
 * @file NusseltBeta3DQGPerWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer for the Beta 3DQG periodic model
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
#include "IoVariable/NusseltBeta3DQGPerWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/NusseltTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   NusseltBeta3DQGPerWriter::NusseltBeta3DQGPerWriter(std::string type)
      : IVariableAsciiEWriter(NusseltTags::BASENAME, NusseltTags::EXTENSION, NusseltTags::HEADER, type, NusseltTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   NusseltBeta3DQGPerWriter::~NusseltBeta3DQGPerWriter()
   {
   }

   void NusseltBeta3DQGPerWriter::init()
   {
      NusseltBeta3DQGPerWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 2);
      NusseltBeta3DQGPerWriter::scalar_iterator sit = sRange.first;

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.cartesian.cartesian_1d");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(1);

      // Get resolution
      pValue = PyLong_FromLong(sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call avgFlux_z
      PythonWrapper::setFunction("avg");
      pValue = PythonWrapper::callFunction(pArgs);

      // Fill matrix and clenup
      PythonWrapper::fillMatrix(this->mNusseltOp, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void NusseltBeta3DQGPerWriter::write()
   {
      // Create file
      this->preWrite();

      NusseltBeta3DQGPerWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 2);
      NusseltBeta3DQGPerWriter::scalar_iterator nlIt;
      NusseltBeta3DQGPerWriter::scalar_iterator tIt;
      NusseltBeta3DQGPerWriter::scalar_iterator it = sRange.first;
      if(it->first == PhysicalNames::TEMPERATURE)
      {
         tIt = it;
         ++it;
         nlIt = it;
      } else
      {
         nlIt = it;
         ++it;
         tIt = it;
      }

      Array nusselt = Array::Zero(2);
      ArrayZ tmp;
      MHDComplex zk_;
      for(int k = 0; k < tIt->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         if(tIt->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
         {
            nusselt(0) += nlIt->second->dom(0).perturbation().point(0,0,k).real();

            if(tIt->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= tIt->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
            {
               zk_ = MHDComplex(0.0, static_cast<MHDFloat>(tIt->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)));
            } else
            {
               zk_ = MHDComplex(0.0, static_cast<MHDFloat>(tIt->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)) - tIt->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL));
            }
            tmp = this->mNusseltOp*tIt->second->dom(0).perturbation().profile(0,k);
            nusselt(1) += (zk_*tmp(0)).real();
         }
      }

      // Get the "global" Nusselt number from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, nusselt.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << 1.0 - nusselt(0)/(1 + nusselt(1)) << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
