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

#include <iostream>

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

      ContinuityWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 3);
      std::vector<FieldComponents::Physical::Id>  comp;
      std::vector<ContinuityWriter::scalar_iterator>  sIt;
      for(ContinuityWriter::scalar_iterator it = sRange.first; it != sRange.second; ++it)
      {
         #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
            if(it->first == PhysicalNames::VELOCITYX)
            {
               comp.push_back(FieldComponents::Physical::TWO);
            } else if(it->first == PhysicalNames::VELOCITYY)
            {
               comp.push_back(FieldComponents::Physical::THREE);
            } else if(it->first == PhysicalNames::VELOCITYZ)
            {
               comp.push_back(FieldComponents::Physical::ONE);
            }
         #else
            if(it->first == PhysicalNames::VELOCITYX)
            {
               comp.push_back(FieldComponents::Physical::ONE);
            } else if(it->first == PhysicalNames::VELOCITYY)
            {
               comp.push_back(FieldComponents::Physical::TWO);
            } else if(it->first == PhysicalNames::VELOCITYZ)
            {
               comp.push_back(FieldComponents::Physical::THREE);
            }
         #endif //GEOMHDISCC_SPATIALSCHEME_TFF

         sIt.push_back(it);
      }
   
      MHDFloat continuity = (sIt.at(0)->second->dom(0).grad().comp(comp.at(0)).data() + sIt.at(1)->second->dom(0).grad().comp(comp.at(1)).data() + sIt.at(2)->second->dom(0).grad().comp(comp.at(2)).data()).array().abs().maxCoeff();

      // Get the "global" Continuity number from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &continuity, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << continuity << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
