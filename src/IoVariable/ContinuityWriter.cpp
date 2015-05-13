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
   
      MHDFloat continuity = (vIt->second->dom(0).grad(FieldComponents::Spectral::X).comp(FieldComponents::Physical::X).data() + vIt->second->dom(0).grad(FieldComponents::Spectral::Y).comp(FieldComponents::Physical::Y).data() + vIt->second->dom(0).grad(FieldComponents::Spectral::Z).comp(FieldComponents::Physical::Z).data()).array().abs().maxCoeff();

      // Get the "global" velocity divergence from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &continuity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

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
         #ifdef GEOMHDISCC_MPI
            MPI_Abort(MPI_COMM_WORLD, 99);
         #endif //GEOMHDISCC_MPI

         throw Exception("Continuity is NaN!");
      }
   }

}
}
