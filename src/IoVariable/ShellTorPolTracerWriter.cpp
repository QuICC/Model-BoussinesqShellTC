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
		:IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix+EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorTracer(Points.rows()), mPolTracer(Points.rows()), mPoints(Points)
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
		PythonWrapper::import("quicc.geometry.spherical.shell_radius");


		// prepare arguments for the linear_2x call
		PyObject *pArgs, * pValue, *pTmp;
		pArgs = PyTuple_New(4);
		pTmp = PyTuple_New(2);

		// set input arguments
		PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second));
		PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second));

		// function call linear_2x
		PythonWrapper::setFunction("linear_r2x");
		pValue = PythonWrapper::callFunction(pTmp);
		MHDFloat a =PyFloat_AsDouble(PyTuple_GetItem(pValue,0));
		MHDFloat b =PyFloat_AsDouble(PyTuple_GetItem(pValue,1));
		Py_DECREF(pTmp);

		// store arguments and prepare for the proj_radial function call
		PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue,0));
		PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue,1));
		int nR = this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
		PyTuple_SetItem(pArgs,0,PyLong_FromLong(nR));

		// TODO: decide how the argument r positions is passed
		Array xRadial = mPoints.col(0);
		int m = xRadial.size();

		//pTmp = PyArray_SimpleNewFromData(1, &m, NPY_DOUBLE, ((void*) xRadial.data()));
		pTmp = PyList_New(m);
		for(int i=0; i<m; ++i){
			PyList_SetItem(pTmp, i, PyFloat_FromDouble(xRadial(i)));
		}
		PyTuple_SetItem(pArgs,3,pTmp);

		// function call proj_radial

		PythonWrapper::import("quicc.projection.shell");
		PythonWrapper::setFunction("proj_radial");
		pValue = PythonWrapper::callFunction(pArgs);

		// Fill matrix and cleanup
		mProjMat = Matrix(m,nR);
		PythonWrapper::getMatrix(mProjMat, pValue);
		Py_DECREF(pValue);

		std::cout << mProjMat << std::endl;

		// function call for the dT_dr and cleanup
		PythonWrapper::setFunction("proj_dradial_dr");
		pValue = PythonWrapper::callFunction(pArgs);
		mProjDrMat = Matrix(m,nR);
		PythonWrapper::getMatrix(mProjDrMat,pValue);
		Py_DECREF(pValue);
		std::cout << mProjDrMat << std::endl;

		// create PyObjects for the 2 vectors theta and phi
		PyObject *vPhi, *vTheta;
		//vTheta = PyArray_SimpleNewFromData(1,&m,NPY_FLOAT64,mPoints.col(1).data());
		//vPhi = PyArray_SimpleNewFromData(1,&m,NPY_FLOAT64,mPoints.col(2).data());
		vPhi = PyList_New(m);
		vTheta = PyList_New(m);
		for(int i = 0; i<m; ++i){
			PyList_SetItem(vTheta,i,  PyFloat_FromDouble(mPoints.col(1)(i)));
			PyList_SetItem(vPhi,i,  PyFloat_FromDouble(mPoints.col(2)(i)));

		}

		//
		PythonWrapper::import("quicc.projection.spherical");
		PythonWrapper::setFunction("lplm");

		// set up the containers for the computed vectors
		Lparts = ArrayMap();

		// precompute the assoc_legendre part of the projector
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
				Array lplm;
				PythonWrapper::getVector(lplm,pValue);
				Lparts[std::make_pair(l,m)] = lplm;

			}
		}
		#endif // GEOMHDISCC_SPATIALSCHEME_SLFM
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
				Array lplm;
				PythonWrapper::getVector(lplm,pValue);
				Lparts[std::make_pair(l,m)] = lplm;

			}
		}
		#endif //GEOMHDISCC_SPATIALSCHEME_SLFL


		PythonWrapper::setFunction("eipm");

		// set up the containers for the computed vector
		Mparts = ArrayZMap();

		// precompute the assoc_legendre part of the projector
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
					ArrayZ eimp;
					PythonWrapper::getVector(eimp,pValue);
					Mparts[std::make_pair(l,m)] = eimp;
				}
			}
			#endif // GEOMHDISCC_SPATIALSCHEME_SLFM
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
					ArrayZ eimp;
					PythonWrapper::getVector(eimp,pValue);
					Mparts[std::make_pair(l,m)] = eimp;
				}
			}
			#endif //GEOMHDISCC_SPATIALSCHEME_SLFL




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
		this->mPolTracer.setZero();

		Array tempPolTraces(this->mPoints.rows());
		#ifdef  QUICC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			double factor = (m==0) ? 1.0 : 2.0;

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				ArrayZ temp = (this->mProjMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).slice(j).col(k)).array()*Mparts[std::make_pair(l,m)].array();
				Array temp2 = (temp.real().array()*Lparts[std::make_pair(l,m)].array()*factor*l*(l+1));
				this->mPolTracer += temp2 ;
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


				ArrayZ temp = (this->mProjMat * vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::POL).slice(j).col(k)).array()*Mparts[std::make_pair(l,m)].array();
				Array temp2 = (temp.real().array()*Lparts[std::make_pair(l,m)].array()*factor*l*(l+1));
				this->mPolTracer += temp2 ;



			}
		}
		#endif //GEOMHDISCC_SPATIALSCHEME_SLFL




	}

	void ShellTorPolTracerWriter::write(){

		// Create file
		this->preWrite();

		// get global values from MPI code
#ifdef GEOMHDISCC_MPI
		Array tracer = this ->mPolTracer;
		MPI_Allreduce(MPI_IN_PLACE, tracer.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		this->mPolTracer = tracer;
#endif // GEOMHDISCC_MPI

		// Check is the workflow allows IO to be performed
		if(FrameworkMacro::allowsIO()){
			this ->mFile << std::setprecision(14) << this->mTime << '\t' << this->mPolTracer.transpose() << std::endl;
		}

		// Close file
		this->postWrite();

	}
}
}



