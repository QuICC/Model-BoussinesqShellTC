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

namespace GeoMHDiSCC{

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
		PythonWrapper::import("geomhdiscc.geometry.spherical.shell_radius");
		PythonWrapper::import("geomhdiscc.projection.shell");
		PythonWrapper::import("geomhdiscc.projection.spherical");

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

		// store arguments and prepare for the proj_radial function call
		PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue,0));
		PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue,1));
		int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>()+2;
		PyTuple_SetItem(pArgs,0,PyLong_FromLong(cols));

		// TODO: decide how the argument r positions is passed
		Array xRadial = mPoints.col(0);
		int m = xRadial.rows();
		int dims[1];
		dims[0]=m;
		pTmp = PyArray_SimpleNewFromData(1,dims,NPY_FLOAT64,xRadial.data());
		PyTuple_SetItem(pArgs,3,pTmp);

		// function call proj_radial
		PythonWrapper::setFunction("proj_radial");
		pValue = PythonWrapper::callFunction(pArgs);

		// Fill matrix and cleanup

		mProjMat = Matrix(m,cols);
		PythonWrapper::getMatrix(mProjMat, pValue);

		// function call for the dT_dr and cleanup
		PythonWrapper::setFunction("proj_dradial_dr");
		pValue = PythonWrapper::callFunction(pArgs);
		mProjDrMat = Matrix(m,cols);
		PythonWrapper::getMatrix(mProjDrMat,pValue);

		// create PyObjects for the 2 vectors theta and phi
		PyObject *vPhi, *vTheta;
		vTheta = PyArray_SimpleNewFromData(1,dims,NPY_FLOAT64,mPoints.col(1).data());
		vPhi = PyArray_SimpleNewFromData(1,dims,NPY_FLOAT64,mPoints.col(2).data());

		// create 2 instances for the execution of eimp and lplm
		PythonWrapper PyWrapL();
		PythonWrapper PyWrapM();
		PyWrapL.import("geomhdiscc.projection.spherical");
		PyWrapM.import("geomhdiscc.projection.spherical");
		PyWrapL.setFunction("lplm");
		PyWrapM.setFunction("eipm");

		// set up the containers for the computed vectors
		Lparts = ArrayMap();
		Mparts = ArrayZMap();


		// precompute the Spherical harmonic part of the projector
		#ifdef GEOMHDISCC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				// prepare the call of lplm
				pArgs = PyTuple_New(3)
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2,vTheta);

				// set the function to lplm
				pValue = PyWrapL.callFunction(pArgs);

				// retrieve the result
				Array lplm;
				PythonWrapper::getVector(lplm,pValue);
				Lparts[std::make_pair(l,m)] = lplm;

				// prepare the call of eimp
				PyTuple_SetItem(pArgs,2, vPhi);

				// set the function to eimp
				pValue = PPyWrapM.callFunction(pArgs);

				// retrieve the result
				Array eimp;
				PythonWrapper::getVector(eimp,pValue);
				Mparts[std::make_pair(l,m)] = eimp;
			}
		}
		#endif // GEOMHDISCC_SPATIALSCHEME_SLFM

		#ifdef GEOMHDISCC_SPATIALSCHEME_SLFL
		for( int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++j){
			// j  is the index to the order l
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(j);

			for(int k=0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(j); ++k){

				// k index is the m degree
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(k, j);

				// prepare the call of lplm
				pArgs = PyTuple_New(3)
				PyTuple_SetItem(pArgs,0,PyLong_FromLong(l));
				PyTuple_SetItem(pArgs,1,PyLong_FromLong(m));
				PyTuple_SetItem(pArgs,2,vTheta);

				// set the function to lplm
				pValue = PyWrapL.callFunction(pArgs);

				// retrieve the result
				Array lplm;
				PythonWrapper::getVector(lplm,pValue);
				Lparts[std::make_pair(l,m)] = lplm;

				// prepare the call of eimp
				PyTuple_SetItem(pArgs,2, vPhi);

				// set the function to eimp
				pValue = PPyWrapM.callFunction(pArgs);

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
		assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::POL)

		// dealias toroidal  variable data
		coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

		// Recover dealiased BWD data
		Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

		// Get FWD storage
		Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

		// Initialize the probes to zero
		this->mPolTracer *= 0.0;

		#ifdef  GEOMHDISCC_SPATIALSCHEME_SLFM
		for( int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k){
			// k  is the index to the degree m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			double factor = (m==0) ? 1.0 : 2.0;

			for(int j=0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j){

				// j index is the l order
				int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				ArrayZ temp = this->mProjMat * rInVarPol.slice(k).col(j);

				this->mPolTracer += (this->Lparts[std::make_pair(l,m)].array()*this->Mparts[std::make_pair(l,m)].array()*temp.array()).real()*factor*l*(l+1);

			}
		}
		#endif //GEOMHDISCC_SPATIALSCHEME_SLFM
		#ifdef  GEOMHDISCC_SPATIALSCHEME_SLFL
		for( int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++j){
			// j  is the index to the order l
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(j);

			for(int k=0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(j); ++k){

				// k index is the m degree
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(k, j);

				double factor = (m==0) ? 1.0 : 2.0;

				ArrayZ temp = this->mProjMat * rInVarPol.slice(j).col(k);

				this->mPolTracer += (this->Lparts[std::make_pair(l,m)].array()*this->Mparts[std::make_pair(l,m)].array()*temp.array()).real()*factor*l*(l+1);



			}
		}
		#endif //GEOMHDISCC_SPATIALSCHEME_SLFL

		// free BWD and TWD storage
		coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPol);
		coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPol);



	}
}
}



