/** 
 * @file ShellTorPolEnergySpectraWriter.cpp
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

// External includes
//

// Class include
//
#include "IoVariable/ShellTorPolEnergySpectraWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

#include<iostream>

namespace QuICC {

namespace IoVariable {

   ShellTorPolEnergySpectraWriter::ShellTorPolEnergySpectraWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(), mPolEnergy(), mTorRadial(), mPolRadial()
   {
   }

   ShellTorPolEnergySpectraWriter::~ShellTorPolEnergySpectraWriter()
   {
   }

   void ShellTorPolEnergySpectraWriter::init()
   {
      // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
      MHDFloat ri = ro*this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;
      this->mVolume = (4.0/3.0)*Math::PI*(std::pow(ro,3) - std::pow(ri,3));

      // obtain Nmax
      int Nmax = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATB1D>();
      mTorRadial = Array(Nmax);
      mPolRadial = Array(Nmax);

      // obtain Lmax and Mmax
      int Lmax = this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL);
      int Mmax = this->mspRes->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::SPECTRAL);

      // resize the mTorEnergy and mPolEnergy matrices
      mTorEnergy = Matrix(Lmax, Mmax);
      mPolEnergy = Matrix(Lmax, Mmax);

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

      // function call linear_r2x
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

   void ShellTorPolEnergySpectraWriter::compute(Transform::TransformCoordinatorType& coord)
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
      this->mTorEnergy.setZero();
      this->mPolEnergy.setZero();

      this->mTorRadial.setZero();
      this->mPolRadial.setZero();

      MHDFloat lfactor = 0.0;
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

               this->mTorEnergy(l,m) += factor*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(j).real()).sum();
               this->mTorRadial += factor*lfactor*(rOutVarTor.slice(k).col(j).real());
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
            	if(m==0){
            		factor = 1.0;
            	} else {
            		factor = 2.0;
            	}
            	this->mTorEnergy(l,m) += factor*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(j).real()).sum();
            	this->mTorRadial += factor*lfactor*(rInVarTor.slice(k).col(j).real());
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

               this->mPolEnergy(l,m) += factor*lfactor*(this->mIntgOp*rInVarPolQ.slice(k).col(j).real()).sum();
               this->mPolRadial += factor*lfactor*(rInVarPolQ.slice(k).col(j).real());
            }
         }
      #endif //defined QUICC_SPATIALSCHEME_SLFM
      #ifdef QUICC_SPATIALSCHEME_SLFL

         lfactor = 0.0;
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = std::pow(l*(l+1.0),2);

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
				if(m==0){
					factor = 1.0;
				} else {
					factor = 2.0;
				}
				this->mPolEnergy(l,m) += factor*lfactor*(this->mIntgOp*rInVarPolQ.slice(k).col(j).real()).sum();
				this->mPolRadial += factor*lfactor*(rInVarPolQ.slice(k).col(j).real());
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

               this->mPolEnergy(l,m) += factor*lfactor*(this->mIntgOp*rInVarPolS.slice(k).col(j).real()).sum();
               this->mPolRadial += factor*lfactor*(rInVarPolS.slice(k).col(j).real());

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

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
				if(m==0){
					factor = 1.0;
				} else {
					factor = 2.0;
				}
				this->mPolEnergy(l,m) += factor*lfactor*(this->mIntgOp*rInVarPolS.slice(k).col(j).real()).sum();
				this->mPolRadial += factor*lfactor*(rInVarPolS.slice(k).col(j).real());

			}
         }
      #endif //QUICC_SPATIALSCHEME_SLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolS);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolS);

      // Normalize by the volume
      this->mPolEnergy /= 2*this->mVolume;
   }

   void ShellTorPolEnergySpectraWriter::write()
   {
      // Create file
      this->preWrite();

      // prepare the vector by either summing on the 1st or 2nd axis

      Array LTorSpectrum = this->mTorEnergy.rowwise().sum();
      Array MTorSpectrum = this->mTorEnergy.colwise().sum();
      Array LPolSpectrum = this->mPolEnergy.rowwise().sum();
      Array MPolSpectrum = this->mPolEnergy.colwise().sum();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
      	  //Array energy(2);
         MPI_Allreduce(MPI_IN_PLACE, LTorSpectrum.data(), LTorSpectrum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, MTorSpectrum.data(), MTorSpectrum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, LPolSpectrum.data(), LPolSpectrum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, MPolSpectrum.data(), MPolSpectrum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         //MPI_Allreduce(MPI_IN_PLACE, mTorRadial.data(), mTorRadial.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         //MPI_Allreduce(MPI_IN_PLACE, mPolRadial.data(), mPolRadial.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         //this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorEnergy + this->mPolEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy << std::endl;
    	 this->mFile << std::setprecision(14) << this->mTime << "\t" << LTorSpectrum.transpose() << '\t';
    	 this->mFile << std::setprecision(14) << '\t' << MTorSpectrum.transpose() << '\t';
    	 this->mFile << std::setprecision(14) << '\t' << LPolSpectrum.transpose() << '\t';
    	 this->mFile << std::setprecision(14) << '\t' << MPolSpectrum.transpose() << std::endl;
    	 //this->mFile << std::setprecision(14) << '\t' << mTorRadial.transpose() << '\t';
    	 //this->mFile << std::setprecision(14) << '\t' << mPolRadial.transpose() << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(isnan(this->mTorEnergy.norm()) || isnan(this->mPolEnergy.norm()))
      {
         FrameworkMacro::abort(99);

         throw Exception("Toroidal/Poloidal energy is NaN!");
      }
   }

}
}
