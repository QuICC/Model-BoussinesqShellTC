/** 
 * @file ShellTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>
#include <cmath>
#include <iostream>

// External includes
//

// Class include
//
#include "IoVariable/ShellTorPolEnergyWriter.hpp"

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

   ShellTorPolEnergyWriter::ShellTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   ShellTorPolEnergyWriter::~ShellTorPolEnergyWriter()
   {
   }

   void ShellTorPolEnergyWriter::init()
   {
      // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
      MHDFloat ri = ro*this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;
      this->mVolume = (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

      // Initialise python wrapper
      PythonWrapper::init();
      PythonWrapper::import("quicc.geometry.spherical.shell_radius");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(4);

      // ... compute a, b factors
      PyObject *pTmp = PyTuple_New(2);
      PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second));
      PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second));
      PythonWrapper::setFunction("linear_r2x");
      pValue = PythonWrapper::callFunction(pTmp);
      PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue, 0));
      PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue, 1));
      // ... create boundray condition (none)
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
      PyTuple_SetItem(pArgs, 3, pValue);

      // Get resolution
      int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>() + 2;
      pValue = PyLong_FromLong(cols);
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call r^2
      PythonWrapper::setFunction("r2");
      pValue = PythonWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpR2(cols,cols);
      PythonWrapper::fillMatrix(tmpR2, pValue);
      Py_DECREF(pValue);

      pTmp = PyTuple_GetSlice(pArgs, 0, 3);
      // Call avg
      PythonWrapper::setFunction("integral");
      pValue = PythonWrapper::callFunction(pTmp);
      // Fill matrix and cleanup
      SparseMatrix tmpAvg(cols,cols);
      PythonWrapper::fillMatrix(tmpAvg, pValue);
      Py_DECREF(pValue);
      PythonWrapper::finalize();

      // Store integral
      this->mIntgOp = tmpAvg.leftCols(cols-2);

      // Store spherical integral (include r^2 factor)
      this->mSphIntgOp = tmpAvg*tmpR2.leftCols(cols-2);

      IVariableAsciiEWriter::init();
   }

   void ShellTorPolEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
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
      coord.transform1D().project(rOutVarTor.rData(), rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Compute |f|^2
      rOutVarTor.rData() = rOutVarTor.rData().array()*rOutVarTor.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarTor.rData(), rOutVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Compute integral over Chebyshev expansion and sum harmonics
      this->mTorEnergy = 0.0;
      this->mPolEnergy = 0.0;
		this->mCentroAntysymEnergy = 0.0;
		this->mCentroSymEnergy = 0.0;
		this->mEquaAntysymEnergy = 0.0;
		this->mEquaSymEnergy = 0.0;

      MHDFloat lfactor = 0.0;

      #ifdef QUICC_SPATIALSCHEME_SLFM
         double factor = 1.0;
         // Loop over harmonic order m
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if( m == 0)
            {
               factor = 1.0;
            } else
            { 
               factor = 2.0;
            }

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = l*(l+1.0);

					MHDFloat ModeEnergy = factor*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(j).real()).sum();
               this->mTorEnergy += ModeEnergy;

					// assign the centro symmetry energies
					if( (l % 2) == 1 )
					{
						this->mCentroSymEnergy += ModeEnergy;
					} else {

						this->mCentroAntysymEnergy += ModeEnergy;
					}

					// assign the equatorial symmetry energies
					if( ((l+m) % 2) == 1 )
					{
						this->mEquaSymEnergy += ModeEnergy;
					} else {

						this->mEquaAntysymEnergy += ModeEnergy;
					}
            }
         }
      #endif //defined QUICC_SPATIALSCHEME_SLFM
      #ifdef QUICC_SPATIALSCHEME_SLFL
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = l*(l+1.0);
            int start = 0;

            int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
            // m = 0, no factor of two
            if( firstM == 0)
            {
               MHDFloat ModeEnergy = lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(0).real())(0);
               this->mTorEnergy += ModeEnergy;
               start = 1;
               firstM = 1;

               if ( (l % 2) == 1){

                  this->mCentroSymEnergy += ModeEnergy;
                  this->mEquaSymEnergy += ModeEnergy;
               } else {

                  this->mCentroAntysymEnergy += ModeEnergy;
                  this->mEquaAntysymEnergy += ModeEnergy;
               }
            }

            Matrix MatrixModes = 2.0*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).rightCols(rInVarTor.slice(k).cols()-start).real());
            this->mTorEnergy += MatrixModes.sum();
            if ( (l % 2) == 1){

					this->mCentroSymEnergy += MatrixModes.sum();
				} else {

					this->mCentroAntysymEnergy += MatrixModes.sum();
				}

				for(int mm = 0; mm < MatrixModes.cols(); ++mm){

				   MHDFloat temp = MatrixModes.col(mm).sum();

				   if( ((l+mm+firstM) % 2) == 1)
				   {
                  this->mEquaSymEnergy += temp;
				   } else {
                  this->mEquaAntysymEnergy += temp;
				   }
				}


         }
      #endif //QUICC_SPATIALSCHEME_SLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTor);

      // Normalize by sphere volume: 4/3*pi*(r_o^3 - r_i^3)
      this->mTorEnergy /= 2*this->mVolume;

      // Dealias poloidal variable data for Q component
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVarPolQ.rData(), rInVarPolQ.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Compute |f|^2
      rOutVarPolQ.rData() = rOutVarPolQ.rData().array()*rOutVarPolQ.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarPolQ.rData(), rOutVarPolQ.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Compute energy in Q component of QST decomposition
      #ifdef QUICC_SPATIALSCHEME_SLFM
         // Loop over harmonic order m
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
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
               lfactor = std::pow(l*(l+1.0),2);

					MHDFloat ModeEnergy = factor*lfactor*(this->mIntgOp*rInVarPolQ.slice(k).col(j).real()).sum();
               this->mPolEnergy += ModeEnergy;

					// assign the centro symmetry energies
					if( (l % 2) == 0 )
					{

						this->mCentroSymEnergy += ModeEnergy;
					} else {

						this->mCentroAntysymEnergy += ModeEnergy;
					}

					// assign the equatorial symmetry energies
					if( ((l+m) % 2) == 0 )
					{
						this->mEquaSymEnergy += ModeEnergy;
					} else {

						this->mEquaAntysymEnergy += ModeEnergy;
					}
            }
         }
      #endif //defined QUICC_SPATIALSCHEME_SLFM
      #ifdef QUICC_SPATIALSCHEME_SLFL
         lfactor = 0.0;
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = std::pow(l*(l+1.0),2);
            int start = 0;
            // m = 0, no factor of two

            int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
            if( firstM == 0)
            {
               MHDFloat ModeEnergy = lfactor*(this->mIntgOp*rInVarPolQ.slice(k).col(0).real())(0);
               this->mPolEnergy += ModeEnergy;
               start = 1;
               firstM = 1;

					if ( (l % 2) == 0){

                  this->mCentroSymEnergy += ModeEnergy;
                  this->mEquaSymEnergy += ModeEnergy;
               } else {

                  this->mCentroAntysymEnergy += ModeEnergy;
                  this->mEquaAntysymEnergy += ModeEnergy;
               }
            }
            Matrix MatrixModes = 2.0*lfactor*(this->mIntgOp*rInVarPolQ.slice(k).rightCols(rInVarPolQ.slice(k).cols()-start).real());
            this->mPolEnergy += MatrixModes.sum();

				if ( (l % 2) == 0 ){

					this->mCentroSymEnergy += MatrixModes.sum();
				} else {

					this->mCentroAntysymEnergy += MatrixModes.sum();
				}

				for(int mm = 0; mm < MatrixModes.cols(); ++mm){

				   MHDFloat temp = MatrixModes.col(mm).sum();

				   if( ((l+mm+firstM) % 2) == 0){

                  this->mEquaSymEnergy += temp;
				   } else {

                  this->mEquaAntysymEnergy += temp;
				   }
				}


         }
      #endif //QUICC_SPATIALSCHEME_SLFL

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
      coord.transform1D().project(rOutVarPolS.rData(), rInVarPolS.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::DIFFR);

      // Compute |f|^2
      rOutVarPolS.rData() = rOutVarPolS.rData().array()*rOutVarPolS.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVarPolS.rData(), rOutVarPolS.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Compute energy in S component of QST decomposition
      #ifdef QUICC_SPATIALSCHEME_SLFM
         // Loop over harmonic order m
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if(m == 0)
            {
               factor = 1.0;
            } else
            { 
               factor = 2.0;
            }

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = l*(l+1.0);

					MHDFloat ModeEnergy = factor*lfactor*(this->mIntgOp*rInVarPolS.slice(k).col(j).real()).sum();
               this->mPolEnergy += ModeEnergy;

					// assign the centro symmetry energies
					if( (l % 2) == 0 )
					{
						this->mCentroSymEnergy += ModeEnergy;
					} else {

						this->mCentroAntysymEnergy += ModeEnergy;
					}

					// assign the equatorial symmetry energies
					if( ((l+m) % 2) == 0 )
					{

						this->mEquaSymEnergy += ModeEnergy;
					} else {

						this->mEquaAntysymEnergy += ModeEnergy;
					}
            }
         }
      #endif //defined QUICC_SPATIALSCHEME_SLFM
      #ifdef QUICC_SPATIALSCHEME_SLFL
         lfactor = 0.0;
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = l*(l+1.0);
            int start = 0;
            // m = 0, no factor of two
            int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
            if(firstM == 0)
            {
               MHDFloat ModeEnergy =  lfactor*(this->mIntgOp*rInVarPolS.slice(k).col(0).real())(0);
               this->mPolEnergy += ModeEnergy;
               start = 1;
               firstM = 1;
					if ( (l % 2) == 0){

                  this->mCentroSymEnergy += ModeEnergy;
                  this->mEquaSymEnergy += ModeEnergy;
               } else {

                  this->mCentroAntysymEnergy += ModeEnergy;
                  this->mEquaAntysymEnergy += ModeEnergy;
               }
            }

            Matrix MatrixModes = 2.0*lfactor*(this->mIntgOp*rInVarPolS.slice(k).rightCols(rInVarPolS.slice(k).cols()-start).real());
            this->mPolEnergy += MatrixModes.sum();

				if ( (l % 2) == 0 ){

					this->mCentroSymEnergy += MatrixModes.sum();
				} else {

					this->mCentroAntysymEnergy += MatrixModes.sum();
				}

				for(int mm = 0; mm < MatrixModes.cols(); ++mm){

				   MHDFloat temp = MatrixModes.col(mm).sum();

				   if( ((l+mm+firstM) % 2) == 0){

                  this->mEquaSymEnergy += temp;
				   } else {

                  this->mEquaAntysymEnergy += temp;
				   }
				}
         }
      #endif //QUICC_SPATIALSCHEME_SLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolS);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolS);

      // Normalize by the volume
      this->mPolEnergy /= 2*this->mVolume;

		// Normalize the remainder of the energies by volume
		this->mEquaSymEnergy /= 2*this->mVolume;
		this->mEquaAntysymEnergy /= 2*this->mVolume;
		this->mCentroSymEnergy /= 2*this->mVolume;
		this->mCentroAntysymEnergy /= 2*this->mVolume;
   }

   void ShellTorPolEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(6);

         energy(0) = this->mTorEnergy;
         energy(1) = this->mPolEnergy;
         energy(2) = this->mCentroSymEnergy;
         energy(3) = this->mCentroAntysymEnergy;
         energy(4) = this->mEquaSymEnergy;
         energy(5) = this->mCentroAntysymEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy = energy(0);
         this->mPolEnergy = energy(1);
		   this->mCentroSymEnergy = energy(2);
         this->mCentroAntysymEnergy = energy(3);
         this->mEquaSymEnergy = energy(4);
         this->mEquaAntysymEnergy = energy(5);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorEnergy + this->mPolEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy <<
							"\t" << this->mCentroSymEnergy <<  "\t" <<  this->mCentroAntysymEnergy << "\t" << this->mEquaSymEnergy << "\t" << this->mEquaAntysymEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(isnan(this->mTorEnergy) || isnan(this->mPolEnergy))
      {
         FrameworkMacro::abort(99);

         throw Exception("Toroidal/Poloidal energy is NaN!");
      }
   }

}
}
;