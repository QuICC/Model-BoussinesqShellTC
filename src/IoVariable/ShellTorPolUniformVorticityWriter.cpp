/** 
 * @file ShellTorPolUniformVorticityWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics uniform vorticity calculation for scalar field in a spherical shell
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>
#include <assert.h>
// External includes
//

// Class include
//
#include "IoVariable/ShellTorPolUniformVorticityWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/AverageTags.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace QuICC {

namespace IoVariable {

   ShellTorPolUniformVorticityWriter::ShellTorPolUniformVorticityWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + AverageTags::BASENAME, AverageTags::EXTENSION, prefix + AverageTags::HEADER, type, AverageTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   ShellTorPolUniformVorticityWriter::~ShellTorPolUniformVorticityWriter()
   {
   }

   void ShellTorPolUniformVorticityWriter::init()
   {
      // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
      MHDFloat rratio = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;
      MHDFloat ri = ro*rratio;

      this->mDelta = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::EKMAN))->second;
      this->mDelta = std::pow(this->mDelta,0.5)*10.;
      //this->mVolume = (std::pow(ro-this->mDelta,5)-std::pow(ri+this->mDelta,5))/(5.* std::sqrt(3.));
      this->mVolume = (std::pow(ro,5)-std::pow(ri,5))/(5.* std::sqrt(3./(4*Math::PI)));

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("quicc.geometry.spherical.shell_radius");

      // Prepare arguments
      PyObject *pArgs, *pValue, *pBC;
      pArgs = PyTuple_New(4);

      // ... compute a, b factors
      PyObject *pTmp = PyTuple_New(2);
      PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(ro));
      PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(rratio));
      PythonWrapper::setFunction("linear_r2x");
      pValue = PythonWrapper::callFunction(pTmp);
      PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue, 0));
      PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue, 1));
      // ... create boundray condition (none)
      pBC = PyDict_New();
      PyDict_SetItem(pBC, PyLong_FromLong(0), PyLong_FromLong(0));
      PyDict_SetItem(pBC, PyUnicode_FromString("cr"), PyLong_FromLong(2));
      PyTuple_SetItem(pArgs, 3, pBC);

      // Get resolution
      int nR = this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      pValue = PyLong_FromLong(nR+2);
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call r^2
      PythonWrapper::setFunction("r2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpR2(nR+2, nR);
      PythonWrapper::fillMatrix(tmpR2, pValue);
      Py_DECREF(pValue);

      // change a couple of arguments
      pValue = PyLong_FromLong(nR+3);
      PyTuple_SetItem(pArgs, 0, pValue);
      pBC = PyDict_New();
      PyDict_SetItem(pBC, PyLong_FromLong(0), PyLong_FromLong(0));
      PyDict_SetItem(pBC, PyUnicode_FromString("cr"), PyLong_FromLong(1));
      PyTuple_SetItem(pArgs, 3, pBC);
      // Call r^1
      PythonWrapper::setFunction("r1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpR1(nR+3, nR+2);
      PythonWrapper::fillMatrix(tmpR1, pValue);
      Py_DECREF(pValue);

      // change a couple of arguments
      pValue = PyLong_FromLong(nR+4);
      PyTuple_SetItem(pArgs, 0, pValue);
      // Call i1
      PythonWrapper::setFunction("i1");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpI1(nR+4, nR+3);
      PythonWrapper::fillMatrix(tmpI1, pValue);
      Py_DECREF(pValue);

      // prepare the points of cut offs
      PyObject * pList;
      pList = PyList_New(2);
      PyList_SetItem(pList, 0, PyFloat_FromDouble(ro-this->mDelta));
      PyList_SetItem(pList, 1, PyFloat_FromDouble(ri+this->mDelta));
      PyTuple_SetItem(pArgs, 3, pList);

      //Call proj
      PythonWrapper::import("quicc.projection.shell");
      PythonWrapper::setFunction("proj");
      pValue = PythonWrapper::callFunction(pArgs);
      SparseMatrix tmpProj(2, nR+4);
      PythonWrapper::fillMatrix(tmpProj, pValue);
      Py_DECREF(pValue);

      // Finalize the Python wrapper
      PythonWrapper::finalize();

      // Store integral projector
      SparseMatrix temp = tmpProj*tmpI1*tmpR1*tmpR2;
      this->mIntgOp  = (temp.row(0)-temp.row(1));
      assert( this->mIntgOp.rows() == 1);
      IVariableAsciiEWriter::init();
   }

   void ShellTorPolUniformVorticityWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

      MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
      MHDFloat ri = ro*this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;

      // Compute integral over Chebyshev expansion and sum harmonics
      this->mUVx = 0.0;
      this->mUVy = 0.0;
      this->mUVz = 0.0;

      #ifdef QUICC_SPATIALSCHEME_SLFM
         double factor = 1.0;
         // Loop over harmonic order m
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
        	int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)

			if(m > 1){
				continue;
			}

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

               if( l!=1){
            	   continue;
               }

				if(m==0){
					this->mUVz = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real());
				} else {
					this->mUVx = - this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real())*std::sqrt(2);
					this->mUVy = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).imag())*std::sqrt(2);
				}
            }
         }
      #endif //defined QUICC_SPATIALSCHEME_SLFM
      #ifdef QUICC_SPATIALSCHEME_SLFL
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if( l!=1){
         	   continue;
            }

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
            	int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);


            	if(m==0){
            		this->mUVz = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real());
            	} else {
            		this->mUVx = - this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real())*std::sqrt(2);
            		this->mUVy = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).imag())*std::sqrt(2);
            	}
            }
         }
      #endif //QUICC_SPATIALSCHEME_SLFL

      // Normalize the integral with int_ri^ror^4dr minus the boundaries
      this->mUVz /= this->mVolume;
      this->mUVx /= this->mVolume;
      this->mUVy /= this->mVolume;
   }

   void ShellTorPolUniformVorticityWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array vorticity(3);

         vorticity(0) = this->mUVx;
         vorticity(1) = this->mUVy;
         vorticity(2) = this->mUVz;

         MPI_Allreduce(MPI_IN_PLACE, vorticity.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mUVx = vorticity(0);
         this->mUVy = vorticity(1);
         this->mUVz = vorticity(2);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mUVx << '\t' << this->mUVy << '\t' << this->mUVz << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mUVz) || std::isnan(this->mUVx) || std::isnan(this->mUVy))
      {
         FrameworkMacro::abort(99);

         throw Exception("Uniform Vorticity is NaN!");
      }
   }

}
}
