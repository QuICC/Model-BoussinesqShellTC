/** \file TruncationPart.cpp
 *  \brief Source of the implementation of the truncation part of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoConfig/ConfigParts/TruncationPart.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string TruncationPart::PARENTTAG = "truncation";

   TruncationPart::TruncationPart(const int dim)
      : IConfigurationPart(TruncationPart::PARENTTAG)
   {
      this->init(dim);
   }

   TruncationPart::~TruncationPart()
   {
   }

   void TruncationPart::init(const int dim)
   {
      // Add first dimension truncation tags
      if(dim > 0)
      {
         this->addIntegerTag("dim1D", -1);
         // this->addIntegerTag("stride1D", -1);
      }

      // Add second dimension truncation tags
      if(dim > 1)
      {
         this->addIntegerTag("dim2D", -1);
         // this->addIntegerTag("stride2D", -1);
      }

      // Add third dimension truncation tags
      if(dim > 2)
      {
         this->addIntegerTag("dim3D", -1);
         // this->addIntegerTag("stride3D", -1);
      }
   }

   void TruncationPart::checkData()
   {
   }

}
}
