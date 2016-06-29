/** 
 * @file Cartesian1DNusseltDZWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer through the Z boundary extracted from temperature field
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
#include "IoVariable/Cartesian1DNusseltDZWriter.hpp"

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

   Cartesian1DNusseltDZWriter::Cartesian1DNusseltDZWriter(std::string type)
      : IVariableAsciiEWriter(NusseltTags::BASENAME, NusseltTags::EXTENSION, NusseltTags::HEADER, type, NusseltTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   Cartesian1DNusseltDZWriter::~Cartesian1DNusseltDZWriter()
   {
   }

   void Cartesian1DNusseltDZWriter::init()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      scalar_iterator sit = sRange.first;

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("geomhdiscc.geometry.cartesian.cartesian_1d");

      // Prepare arguments
      PyObject *pArgs, *pValue, *pTmp;
      pArgs = PyTuple_New(2);

      // Get resolution
      pValue = PyLong_FromLong(sit->second->dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Set boundary condition
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(21));
      pTmp = PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::SCALE1D))->second);
      PyDict_SetItem(pValue, PyUnicode_FromString("c"), pTmp);
      Py_INCREF(Py_False);
      PyDict_SetItem(pValue, PyUnicode_FromString("use_parity"), Py_False);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Call zblk and use derivative boundary condition
      PythonWrapper::setFunction("zblk");
      pValue = PythonWrapper::callFunction(pArgs);

      // Fill matrix and cleanup
      PythonWrapper::fillMatrix(this->mNusseltOp, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      IVariableAsciiEWriter::init();
   }

   void Cartesian1DNusseltDZWriter::write()
   {
      // Create file
      this->preWrite();

      Cartesian1DNusseltDZWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      Cartesian1DNusseltDZWriter::scalar_iterator sit = sRange.first;

      ArrayI mode = sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(0);
      MHDFloat nusselt = 0.0;
      if(mode(2) == 0 && mode(3) == 0)
      {
         // Compute Nusselt number
         nusselt = -(this->mNusseltOp*sit->second->dom(0).perturbation().profile(0,0).real()).row(0)(0);
      }

      // Get the "global" Nusselt number from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &nusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << 1.0 + nusselt << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if Nusselt number is NaN
      if(std::isnan(nusselt))
      {
         FrameworkMacro::abort(99);

         throw Exception("Nusselt number is NaN!");
      }
   }

}
}
