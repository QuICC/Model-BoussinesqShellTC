/** 
 * @file ShellTorPolTorqueWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical shell
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
#include "IoVariable/ShellTorPolTorqueWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace QuICC {

namespace IoVariable {

   ShellTorPolTorqueWriter::ShellTorPolTorqueWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   ShellTorPolTorqueWriter::~ShellTorPolTorqueWriter()
   {
   }

   void ShellTorPolTorqueWriter::init()
   {
	   this->mComputeFlag = false;
#ifdef QUICC_SPATIALSCHEME_SLFM
	 double factor = 1.0;
	 // Loop over harmonic order m
	 for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
	 {
		// determine current m
		int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
		// m = 0, no factor of two


		for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
		{
		   int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

		   if( m==0 && l==1){
			   this->mComputeFlag = true;
		   }
		}
	 }
#endif //defined QUICC_SPATIALSCHEME_SLFM
#ifdef QUICC_SPATIALSCHEME_SLFL
	 double factor=1.;
	 // Loop over harmonic degree l
	 for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
	 {
		int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);


		for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

			if( m ==0 && l==1){
				this->mComputeFlag = true;
			}
		}

	 }
#endif //QUICC_SPATIALSCHEME_SLFL

	 if(this->mComputeFlag){
		  // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
		  MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
		  MHDFloat ri = ro*this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;
		  this->mVolume = (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

		  // Initialise python wrapper
		  PythonWrapper::init();
		  PythonWrapper::import("quicc.geometry.spherical.shell_radius");

		  // Prepare arguments
		  PyObject *pArgs, *pValue, *pABBoundary;
		  pArgs = PyTuple_New(4);

		  // ... compute a, b factors
		  PyObject *pTmp = PyTuple_New(2);
		  PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(ro);
		  PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second));
		  PythonWrapper::setFunction("linear_r2x");
		  pValue = PythonWrapper::callFunction(pTmp);
		  MHDFloat a = PyTuple_GetItem(pValue,0);
		  MHDFloat b = PyTuple_GetItem(pValue,1);

		  // prepare the next function call
		  pArgs = PyTuple_New(2);

		  // create boundray conditions (none)
		  pValue = PyDict_New();

		  PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(20));
		  PyTuple_SetItem(pArgs, 1, pValue);

		  // Get resolution
		  //int nR = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>() + 2;
		  int nR = this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);


		  PyTuple_SetItem(pArgs, 0, PyLong_FromLong(nR));

		  // Call zblk
		  PythonWrapper::setFunction("zblk");
		  pValue = PythonWrapper::callFunction(pArgs);
		  // Fill matrix and cleanup
		  this->mProj = SparseMatrix(nR,nR);
		  PythonWrapper::fillMatrix(this->mProj, pValue);


		  pTmp = PyTuple_GetSlice(pArgs, 0, 3);

		  // add specifics to boundary conditions
		   = PyDict_New();
		  PyDict_SetItem(pABBoundary, PyString_FromString("a"), PyFloat_FromDouble(a));
		  PyDict_SetItem(pABBoundary, PyString_FromString("b"), PyFloat_FromDouble(b));
		  PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(21));
		  PyTuple_SetItem(pArgs, 1, pValue);
		  PyDict_SetItem(pValue, PyString_FromString("c"), pABBoundary);

		  // Call avg
		  PythonWrapper::setFunction("zblk");
		  pValue = PythonWrapper::callFunction(pTmp);
		  // Fill matrix and cleanup
		  this->mProjDr = SparseMatrix(nR,nR);
		  PythonWrapper::fillMatrix(this->mProjDr, pValue);
		  Py_DECREF(pValue);
		  Py_DECREF(pABBoundary);
		  PythonWrapper::finalize();
	 }


	  // call init in the base class
	  IVariableAsciiEWriter::init();
   }

   void ShellTorPolTorqueWriter::compute(Transform::TransformCoordinatorType& coord)
   {

	   // for the detection of the core carrying the Toroidal 1/0 function


	   /*
	    * compute the "bad way"
	    */

	   if(this->mComputeFlag){
	  // get iterator to field
	  vector_iterator vIt;
	  vector_iterator_range vRange = this->vectorRange();
	  assert(std::distance(vRange.first, vRange.second) == 1);
	  assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
	  assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

	  // retrieve ri
	  MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
	  MHDFloat ri = ro*this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;


	  #ifdef QUICC_SPATIALSCHEME_SLFM
		 double factor = 1.0;
		 // Loop over harmonic order m
		 for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
		 {
			// determine current m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
			// m = 0, no factor of two
			if( m == 0)
			{
			   factor = 1.0;
			} else
			{
			   factor = 2.0;
			}

			for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
			{
			   int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
			   lfactor = l*(l+1.0);

			   Array temp = -this->mProjDr*vRange.first.slice(k).col(j).real()*ri+2*this->mProj*vRange.first.slice(k).col(j).real();
			   this->mTorque = temp[1];

			}
		 }
	  #endif //defined QUICC_SPATIALSCHEME_SLFM
	  #ifdef QUICC_SPATIALSCHEME_SLFL
		 double factor=1.;
		 // Loop over harmonic degree l
		 for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
		 {
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
			lfactor = l*(l+1.0);

			for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				Array temp = -this->mProjDr*vRange.first.slice(k).col(j).real()*ri+2*this->mProj*vRange.first.slice(k).col(j).real();
				this->mTorque = temp[1];
			}

		 }
	  #endif //QUICC_SPATIALSCHEME_SLFL

	   }



   }

   void ShellTorPolTorqueWriter::write()
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

         /*
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
      }*/
   }

}
}
