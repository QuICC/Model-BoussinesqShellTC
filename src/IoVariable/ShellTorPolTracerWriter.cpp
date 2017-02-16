/*
 * ShellTorPolTracerWriter.cpp
 *
 *  Created on: Nov 23, 2016
 *      Author: Nicolò Lardelli
 */

// Configuration includes
//
#include"Framework/FrameworkMacro.h"

// System includes
//
#include<iomanip>

// Class include
//
#include"IoVariable/ShellTorPolTracerWriter.hpp"

// Project includes
//
#include"Enums/Dimensions.hpp"
#include"Enums/FieldIds.hpp"
#include"IoTools/IdToHuman.hpp"
#include"IoVariable/EnergyTags.hpp"
#include"TypeSelectors/ScalarSelector.hpp"
#include"Python/PythonWrapper.hpp"
#include <iostream>

namespace QuICC{

namespace IoVariable{

	/*
	 * @brief Constructor
	 */
	ShellTorPolTracerWriter::ShellTorPolTracerWriter(const std::string& prefix, const std::string& type, const Matrix& Points)
		:IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix+EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mRpart(Points.rows()), mThetapart(Points.rows()), mPhipart(Points.rows()), mPoints(Points)
	{

	}

	/*
	 * @brief Destructor
	 */

	ShellTorPolTracerWriter::~ShellTorPolTracerWriter(){}


	/*
	 * @brief Subroutine init: generates the projections matrices
	 */

	void ShellTorPolTracerWriter::init(){

		// Retrieve size of the outer and inner shell
		MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
		MHDFloat rratio = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;

		// initialize the Python wrapper
		PythonWrapper::init();


		// prepare arguments for the linear_2x call
		PyObject *pArgs, * pValue, *pTmp;
		pArgs = PyTuple_New(4);
		pTmp = PyTuple_New(2);

		// set input arguments
		PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second));
		PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second));

		// import python module shell_radius and function call linear_r2x
		PythonWrapper::import("quicc.geometry.spherical.shell_radius");
		PythonWrapper::setFunction("linear_r2x");
		pValue = PythonWrapper::callFunction(pTmp);
		//MHDFloat a =PyFloat_AsDouble(PyTuple_GetItem(pValue,0));
		//MHDFloat b =PyFloat_AsDouble(PyTuple_GetItem(pValue,1));

		// store arguments and prepare for the proj_radial function call
		PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue,0));
		PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue,1));
		int nR = this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
		PyTuple_SetItem(pArgs,0,PyLong_FromLong(nR));

		// prepare the mPoints[:,0] argument for python
		Array xRadial = mPoints.col(0);
		int n = xRadial.size();

		// prepare the list to pass
		pTmp = PyList_New(n);
		for(int i=0; i<n; ++i){
			PyList_SetItem(pTmp, i, PyFloat_FromDouble(xRadial(i)));
		}
		PyTuple_SetItem(pArgs,3,pTmp);

		// load module quicc.projection.shell
		PythonWrapper::import("quicc.projection.shell");


		// function call proj_radial
		PythonWrapper::setFunction("proj_radial");
		pValue = PythonWrapper::callFunction(pArgs);
		pTmp = PyList_New(n);
		for(int i=0; i<n; ++i){
			PyList_SetItem(pTmp, i, PyFloat_FromDouble(xRadial(i)));
		}
		PyTuple_SetItem(pArgs,3,pTmp);
		//Py_INCREF(pTmp);

		// Fill matrix mProjMat and cleanup
		this->mProjMat = Matrix(n,nR);
		PythonWrapper::getMatrix(mProjMat, pValue);
		//Py_DECREF(pValue);

		// function call for the dT_dr part (proj_dradial_dr)
		PythonWrapper::setFunction("proj_dradial_dr");
		pValue = PythonWrapper::callFunction(pArgs);
		//Py_INCREF(pTmp);

		// Fill mProjDrMat and cleanup
		this->mProjDrMat = Matrix(n,nR);
		PythonWrapper::getMatrix(mProjDrMat,pValue);
		Py_DECREF(pValue);
		//Py_INCREF(pTmp);

		// function call for the T/r part (proj_radial_r)
		PythonWrapper::setFunction("proj_radial_r");
		pValue = PythonWrapper::callFunction(pArgs);


		// fill matrix mProjInvrMat and cleanup
		this->mProjInvrMat = Matrix(n, nR);
		PythonWrapper::getMatrix(mProjInvrMat, pValue);
		Py_DECREF(pValue);

		//Py_DECREF(pArgs);
		// create PyObjects for the 2 vectors theta and phi
		PyObject *vPhi, *vTheta;
		vPhi = PyList_New(n);
		vTheta = PyList_New(n);
		for(int i = 0; i<n; ++i){
			PyList_SetItem(vTheta,i,  PyFloat_FromDouble(mPoints.col(1)(i)));
			PyList_SetItem(vPhi,i,  PyFloat_FromDouble(mPoints.col(2)(i)));

		}

		// load quicc.projection spherical and prepare the eipm function call
		PythonWrapper::import("quicc.projection.spherical");
		PythonWrapper::setFunction("eipm");

		// prepare the temporary container Mparts
		ArrayZMap Mparts = ArrayZMap();

		// precompute the exponential aximuthal part of Y_l^m and friends
		#ifdef QUICC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				// prepare the call of eimp
				pArgs = PyTuple_New(3);
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2, vPhi);

				// set the function to eimp
				pValue = PythonWrapper::callFunction(pArgs);

				// retrieve the result
				ArrayZ eimp(n);
				PythonWrapper::getVector(eimp,pValue);
				Mparts[std::make_pair(l,m)] = eimp;
			}
		}
		#endif // QUICC_SPATIALSCHEME_SLFM
		#ifdef QUICC_SPATIALSCHEME_SLFL
		for( int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++j){
			// j  is the index to the order l
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(j);

			for(int k=0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(j); ++k){

				// k index is the m degree
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(k, j);

				// prepare the call of eimp
				pArgs = PyTuple_New(3);
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2, vPhi);

				// set the function to eimp
				pValue = PythonWrapper::callFunction(pArgs);

				// retrieve the result
				ArrayZ eimp(n);
				PythonWrapper::getVector(eimp,pValue);
				Mparts[std::make_pair(l,m)] = eimp;
			}
		}
		#endif //QUICC_SPATIALSCHEME_SLFL



		PythonWrapper::setFunction("lplm");

		// set up the containers for the computed vectors
		YlmParts = ArrayZMap();

		// precompute the Ylm of the projector
		#ifdef QUICC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				// prepare the call of lplm
				pArgs = PyTuple_New(3);
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2,vTheta);

				// set the function to lplm
				pValue = PythonWrapper::callFunction(pArgs);

				// retrieve the result
				Array lplm(n);
				PythonWrapper::getVector(lplm,pValue);

				// set the Ylm
				YlmParts[std::make_pair(l,m)] = lplm.array()*Mparts[std::make_pair(l,m)].array();

			}
		}
		#endif // QUICC_SPATIALSCHEME_SLFM
		#ifdef QUICC_SPATIALSCHEME_SLFL
		for( int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++j){
			// j  is the index to the order l
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(j);

			for(int k=0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(j); ++k){

				// k index is the m degree
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(k, j);

				// prepare the call of lplm
				pArgs = PyTuple_New(3);
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2,vTheta);

				// set the function to lplm
				pValue = PythonWrapper::callFunction(pArgs);

				// retrieve the result
				Array lplm(n);
				PythonWrapper::getVector(lplm,pValue);

				// set the Ylm
				YlmParts[std::make_pair(l,m)] = lplm.array()*Mparts[std::make_pair(l,m)].array();

			}
		}
		#endif //QUICC_SPATIALSCHEME_SLFL

		PythonWrapper::setFunction("dplm");

		// set up the containers for the computed vectors
		dYlmdthParts = ArrayZMap();

		// precompute the d Ylm d theta part of the projector
		#ifdef QUICC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				// prepare the call of lplm
				pArgs = PyTuple_New(3);
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2,vTheta);

				// set the function to lplm
				pValue = PythonWrapper::callFunction(pArgs);

				// retrieve the result
				Array dplm(n);
				PythonWrapper::getVector(dplm,pValue);

				// set the dYlm /d theta
				dYlmdthParts[std::make_pair(l,m)] = dplm.array()*Mparts[std::make_pair(l,m)].array();
				std::cout << l << '\t' << m << '\t' << "dYlm" << dYlmdthParts[std::make_pair(l,m)] << std::endl;

			}
		}
		#endif // QUICC_SPATIALSCHEME_SLFM
		#ifdef QUICC_SPATIALSCHEME_SLFL
		for( int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++j){
			// j  is the index to the order l
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(j);

			for(int k=0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(j); ++k){

				// k index is the m degree
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(k, j);

					// prepare the call of lplm
					pArgs = PyTuple_New(3);
					PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
					PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
					PyTuple_SetItem(pArgs,2,vTheta);

					// set the function to lplm
					pValue = PythonWrapper::callFunction(pArgs);

					// retrieve the result
					Array dplm(n);
					PythonWrapper::getVector(dplm,pValue);

					// set the dYlm /d theta
					dYlmdthParts[std::make_pair(l,m)] = dplm.array()*Mparts[std::make_pair(l,m)].array();
					std::cout << l << '\t' << m << '\t' << "dYlm" << dYlmdthParts[std::make_pair(l,m)] << std::endl;

			}
		}
		#endif //QUICC_SPATIALSCHEME_SLFL

		PythonWrapper::setFunction("lplm_sin");

		// set up the containers for the computed vectors
		YlmIsinParts = ArrayZMap();

		// precompute the d Ylm d theta part of the projector
		#ifdef QUICC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				// prepare the call of lplm
				pArgs = PyTuple_New(3);
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2,vTheta);

				// set the function to lplm
				pValue = PythonWrapper::callFunction(pArgs);

				// retrieve the result
				Array lplm_sin(n);
				PythonWrapper::getVector(lplm_sin,pValue);

				// set the dYlm /d theta
				YlmIsinParts[std::make_pair(l,m)] = lplm_sin.array()*Mparts[std::make_pair(l,m)].array()*m*1i;
				std::cout << l << '\t' << m << '\t' << "YlmIs" << YlmIsinParts[std::make_pair(l,m)].transpose() << std::endl;

			}
		}
		#endif // QUICC_SPATIALSCHEME_SLFM
		#ifdef QUICC_SPATIALSCHEME_SLFL
		for( int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++j){
			// j  is the index to the order l
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(j);

			for(int k=0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(j); ++k){

				// k index is the m degree
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(k, j);

				// prepare the call of lplm
				pArgs = PyTuple_New(3);
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2,vTheta);

				// set the function to lplm
				pValue = PythonWrapper::callFunction(pArgs);

				// retrieve the result
				Array lplm_sin(n);
				PythonWrapper::getVector(lplm_sin,pValue);

				// set the dYlm /d theta
				YlmIsinParts[std::make_pair(l,m)] = lplm_sin.array()*Mparts[std::make_pair(l,m)].array()*m*1i;
				std::cout << l << '\t' << m << '\t' << "YlmIs" << YlmIsinParts[std::make_pair(l,m)].transpose() << std::endl;
			}
		}
		#endif //QUICC_SPATIALSCHEME_SLFL





		IVariableAsciiEWriter::init();
	}

	/*
	 * @brief
	 */

	void ShellTorPolTracerWriter::compute(Transform::TransformCoordinatorType& coord){

		// get iterator to field
		vector_iterator vIt;
		vector_iterator_range vRange = this->vectorRange();
		assert(std::distance(vRange.first,vRange.second)==1);
		assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
		assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

		// Initialize the probes to zero
		this->mRpart.setZero();
		this->mThetapart.setZero();
		this->mPhipart.setZero();

		Array tempPolTraces(this->mPoints.rows());
		#ifdef  QUICC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			double factor = (m==0) ? 1.0 : 2.0;

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				// compute R component from poloidal part
				ArrayZ tempR = (this->mProjMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).slice(j).col(k)).array()*YlmParts[std::make_pair(l,m)].array();
				this->mRpart += tempR.real()*factor*l*(l+1));

				ArrayZ tempPol = this->mProjDrMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).slice(j).col(k);
				ArrayZ tempTor = this->mProjInvrMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(j).col(k);

				Array temp1 = (tempPol.array()*dYlmdthParts[std::make_pair(l,m)].array()+tempTor.array()*YlmIsinParts[std::make_pair(l,m)].array()).real();
				this->mThetapart += temp1*factor;

				Array temp2 = (tempPol.array()*YlmIsinParts[std::make_pair(l,m)].array() - tempTor.array()*dYlmdthParts[std::make_pair(l,m)].array()).real();
				this->mPhipart += temp2*factor;

			}
		}
		#endif //GEOMHDISCC_SPATIALSCHEME_SLFM
		#ifdef  QUICC_SPATIALSCHEME_SLFL
		for( int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++j){
			// j  is the index to the order l
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(j);

			for(int k=0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(j); ++k){

				// k index is the m degree
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(k, j);

				double factor = (m==0) ? 1.0 : 2.0;

				// compute R component from poloidal part
				ArrayZ tempR = (mProjMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).slice(j).col(k)).array()*YlmParts[std::make_pair(l,m)].array();
				this->mRpart += tempR.real()*factor*l*(l+1);
				std::cout << l << '\t' << m << '\t' << (mProjMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).slice(j).col(k)).transpose() << std::endl;

				ArrayZ tempPol = this->mProjDrMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).slice(j).col(k);
				ArrayZ tempTor = this->mProjInvrMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(j).col(k);

				Array temp1 = (tempPol.array()*dYlmdthParts[std::make_pair(l,m)].array()+tempTor.array()*YlmIsinParts[std::make_pair(l,m)].array()).real();
				this->mThetapart += temp1*factor;

				Array temp2 = (tempPol.array()*YlmIsinParts[std::make_pair(l,m)].array() - tempTor.array()*dYlmdthParts[std::make_pair(l,m)].array()).real();
				this->mPhipart += temp2*factor;

			}
		}
		#endif //QUICC_SPATIALSCHEME_SLFL




	}

	void ShellTorPolTracerWriter::write(){

		// Create file
		this->preWrite();

		// get global values from MPI code
		#ifdef QUICC_MPI

		// Allreduce the R components
		Array R = this ->mRpart;
		MPI_Allreduce(MPI_IN_PLACE, R.data(), R.rows(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		this->mPolTracer = R;

		// Allreduce the Theta component
		Array Theta = this->mThetapart;
		MPI_Allreduce(MPI_IN_PLACE, Theta.data(), Theta.rows(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		this->mThetapart = Theta;

		// Allreduce the Phi component
		Array Phi = this->mPhipart;
		MPI_Allreduce(MPI_IN_PLACE, Phi.data(), Phi.rows(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		this->mPhipart = Phi;
		#endif // QUICC_MPI

		// Check is the workflow allows IO to be performed
		if(FrameworkMacro::allowsIO()){
			this ->mFile << std::setprecision(14) << this->mTime << '\t' << this->mRpart.transpose();
			this ->mFile << std::setprecision(14) << this->mThetapart.transpose() << '\t' << this->mPhipart.transpose()<< std::endl;
		}

		// Close file
		this->postWrite();

	}
}
}



