/** 
 * @file TruncationPart.cpp
 * @brief Source of the implementation of the truncation part of the configuration
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "IoConfig/ConfigParts/TruncationPart.hpp"

// Project includes
//

namespace QuICC {

namespace IoConfig {

   const std::string TruncationPart::PARENTTAG = "truncation";

   TruncationPart::TruncationPart(const int dim, const std::vector<bool>& isPeriodicBox)
      : IConfigurationPart(TruncationPart::PARENTTAG)
   {
      this->init(dim, isPeriodicBox);
   }

   TruncationPart::~TruncationPart()
   {
   }

   void TruncationPart::init(const int dim, const std::vector<bool>& isPeriodicBox)
   {
      // Safety assert
      assert(isPeriodicBox.size() == static_cast<size_t>(dim));

      // Add first dimension truncation tags
      if(dim > 0)
      {
         this->addIntegerTag("dim1D", -1);

         if(isPeriodicBox.at(0))
         {
            this->addFloatTag("kc1D", -1);
            this->addFloatTag("box1D", -1);
         }
      }

      // Add second dimension truncation tags
      if(dim > 1)
      {
         this->addIntegerTag("dim2D", -1);

         if(isPeriodicBox.at(1))
         {
            this->addFloatTag("kc2D", -1);
            this->addFloatTag("box2D", -1);
         }
      }

      // Add third dimension truncation tags
      if(dim > 2)
      {
         this->addIntegerTag("dim3D", -1);

         if(isPeriodicBox.at(2))
         {
            this->addFloatTag("kc3D", -1);
            this->addFloatTag("box3D", -1);
         }
      }
   }

   void TruncationPart::checkData()
   {
   }

}
}
