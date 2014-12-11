/** 
 * @file NusseltCubicWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer for a 3D finite box
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
#include "IoVariable/NusseltCubicWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/NusseltTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"
#include "Base/DecoupledComplexInternal.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   NusseltCubicWriter::NusseltCubicWriter(std::string type)
      : IVariableAsciiEWriter(NusseltTags::BASENAME, NusseltTags::EXTENSION, NusseltTags::HEADER, type, NusseltTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   NusseltCubicWriter::~NusseltCubicWriter()
   {
   }

   void NusseltCubicWriter::init()
   {
      NusseltCubicWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      NusseltCubicWriter::scalar_iterator sit = sRange.first;

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.cartesian.cartesian_3d");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(4);

      // Get resolution
      pValue = PyLong_FromLong(sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);
      pValue = PyLong_FromLong(sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 1, pValue);
      pValue = PyLong_FromLong(sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 2, pValue);
      // Set scale
      pValue = PyLong_FromLong(2.0);
      PyTuple_SetItem(pArgs, 3, pValue);

      // Call avgFlux_z
      PythonWrapper::setFunction("avgFlux_z");
      pValue = PythonWrapper::callFunction(pArgs);

      // Fill matrix and clenup
      PythonWrapper::fillMatrix(this->mNusseltOp, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void NusseltCubicWriter::write()
   {
      // Create file
      this->preWrite();

      NusseltCubicWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      NusseltCubicWriter::scalar_iterator sit = sRange.first;

      // Copy data
      int l;
      int k_;
      int j_;
      int dimK = sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
      int dimJ = sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      Array field(dimK*sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL));
      for(int k = 0; k < sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
      {
         k_ = sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
         for(int j = 0; j < sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
         {
            j_ = sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
            for(int i = 0; i < sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
            {
               // Compute correct position
               l = k_ + j_ + i;

               // Copy field value into field
               Datatypes::internal::setScalar(field, l, sit->second->dom(0).perturbation().point(i,j,k));
            }
         }
      }

      // Get the "global" field
      #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE 
         MPI_Reduce(MPI_IN_PLACE, field.data(), field.rows(), MPI_DOUBLE, MPI_SUM, 0, GeoMHDiSCC::FramworkMacro::spectralComm());
      #endif //defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE

      Matrix nusselt = -this->mNusseltOp*field;
      assert(nusselt.rows() == nusselt.cols() && nusselt.rows() == 1);
      nusselt(0,0) += 1;

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << nusselt(0,0) << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
