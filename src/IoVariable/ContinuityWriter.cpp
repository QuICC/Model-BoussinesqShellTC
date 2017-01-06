/** 
 * @file ContinuityWriter.cpp
 * @brief Source of the implementation of the ASCII maximal continuity value writer
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
#include "IoVariable/ContinuityWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/ContinuityTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   ContinuityWriter::ContinuityWriter(std::string type)
      : IVariableAsciiEWriter(ContinuityTags::BASENAME, ContinuityTags::EXTENSION, ContinuityTags::HEADER, type, ContinuityTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   ContinuityWriter::~ContinuityWriter()
   {
   }

   void ContinuityWriter::write()
   {
      // Create file
      this->preWrite();

      ContinuityWriter::vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      ContinuityWriter::vector_iterator  vIt = vRange.first;

      #ifdef QUICC_SPATIALDIMENSION_3D
         MHDFloat continuity = (vIt->second->dom(0).grad(FieldComponents::Spectral::ONE).comp(FieldComponents::Physical::ONE).data() + vIt->second->dom(0).grad(FieldComponents::Spectral::TWO).comp(FieldComponents::Physical::TWO).data() + vIt->second->dom(0).grad(FieldComponents::Spectral::THREE).comp(FieldComponents::Physical::THREE).data()).array().abs().maxCoeff();
      #else
         MHDFloat continuity = (vIt->second->dom(0).grad(FieldComponents::Spectral::ONE).comp(FieldComponents::Physical::ONE).data() + vIt->second->dom(0).grad(FieldComponents::Spectral::TWO).comp(FieldComponents::Physical::TWO).data()).array().abs().maxCoeff();
      #endif //QUICC_SPATIALDIMENSION_3D

      // Get the "global" velocity divergence from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &continuity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << continuity << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if continuity is NaN
      if(std::isnan(continuity))
      {
         FrameworkMacro::abort(99);

         throw Exception("Continuity is NaN!");
      }
   }

}
}
