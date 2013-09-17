/** 
 * @file HumanToId.cpp
 * @brief Source of human strings to enum id converters
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
// Configuration includes
//
#include "Exceptions/Exception.hpp"

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoTools/HumanToId.hpp"

// Project includes
//
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoTools {

   NonDimensional::Id HumanToId::toNd(const std::string& id)
   {
      if(id == IdToHuman::toTag(NonDimensional::EKMAN))
      {
         return NonDimensional::EKMAN;
      } else if(id == IdToHuman::toTag(NonDimensional::ROBERTS))
      {
         return NonDimensional::ROBERTS;
      } else if(id == IdToHuman::toTag(NonDimensional::RAYLEIGH))
      {
         return NonDimensional::RAYLEIGH;
      } else if(id == IdToHuman::toTag(NonDimensional::MAGEKMAN))
      {
         return NonDimensional::MAGEKMAN;
      } else if(id == IdToHuman::toTag(NonDimensional::PRANDTL))
      {
         return NonDimensional::PRANDTL;
      } else if(id == IdToHuman::toTag(NonDimensional::MAGPRANDTL))
      {
         return NonDimensional::MAGPRANDTL;
      } else if(id == IdToHuman::toTag(NonDimensional::CHI))
      {
         return NonDimensional::CHI;
      } else if(id == IdToHuman::toTag(NonDimensional::GAMMA))
      {
         return NonDimensional::GAMMA;
      } else if(id == IdToHuman::toTag(NonDimensional::GAPWIDTH))
      {
         return NonDimensional::GAPWIDTH;
      } else if(id == IdToHuman::toTag(NonDimensional::RRATIO))
      {
         return NonDimensional::RRATIO;
      } else
      {
         throw Exception("Unknown string to ID conversion requested");
      }
   }
}

}
