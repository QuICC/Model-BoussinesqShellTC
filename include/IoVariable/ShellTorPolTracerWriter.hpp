/**
 * @file ShellTorPolTracerWriter.hpp
 * @brief Implementation of the ASCII tracer writer for Toroidal-Poloidal decomposed fields in a spherical shell
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

#ifndef SHELLTORPOLTRACERWRITER_HPP
#define SHELLTORPOLTRACERWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//

#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include<utility>
#include<map>

namespace QuICC{

namespace IoVariable{


	/*
	 * @brief Implementation of the ASCII tracer writer for Toroidal-Poloidal decomposed fields in a spherical shell
	 */
	class ShellTorPolTracerWriter : public IVariableAsciiEWriter
	{
	public:

		/*
		 * @brief Constructor of ShellTorPolTracerWriter
		 *
		 * @param prefix Prefix to use for the file name
		 * @param type Type of the file (typically scheme name)
		 */
		ShellTorPolTracerWriter(const std::string& prefix, const std::string& type, const Matrix& Points);

		/*
		 * @brief Destructor
		 */
		~ShellTorPolTracerWriter();

		/*
		 * @brief Initialize the operator, transform and file
		 */
		virtual void init();

		/*
		 * @brief Compute energy for TorPol decomposed field
		 */
		virtual void compute(Transform::TransformCoordinatorType& coord);

		/*
		 * brief Write to state file
		 */
		virtual void write();

		/*
		 *@brief Requires heavy computations?
		 */
		virtual bool isHeavy() const;

	protected:

	private:

		/*
		 * @brief Vector that stores the projection of the Toroidal field
		 */
		Array mTorTracer;

		/*
		 * @brief Vector that stores the projection of the Poloidal field
		 */
		Array mPolTracer;

		// TODO: are other fields that we need?
		// TODO: mPoints (given at constructor time?),

		/*
		 * @briefs Position of the probes, size(mPoints)=[M,3], written in (r,theta,phi)
		 */
		Matrix mPoints;

		/*
		 * @brief Precomputed M part of the spherical harmonic
		 */
		typedef std::map<std::pair<int,int>,Array> ArrayMap;
		typedef std::map<std::pair<int,int>,ArrayZ> ArrayZMap;
		ArrayZMap Mparts;

		/*
		 * @brief Precomputed L part of the spherical harmonic
		 */
		ArrayMap Lparts;

		/*
		 * @brief
		 */
		Matrix mProjMat;

		/*
		 * @brief
		 */
		Matrix mProjDrMat;

		/*
		 * If we decide that we need current Tracers as well
		Matrix mProjTricurlMat;
		 */
	};




	inline bool ShellTorPolTracerWriter::isHeavy() const
	{
		return true;
	}

	typedef SharedPtrMacro<ShellTorPolTracerWriter> SharedShellTorPolTracerWriter;




}

}

#endif
