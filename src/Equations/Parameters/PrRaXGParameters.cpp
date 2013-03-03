/** \file PrRaXGParameters.cpp
 *  \brief Source of the implementation of the PrRaXG non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/PrRaXGParameters.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   PrRaXGParameters::PrRaXGParameters()
      : IEquationParameters()
   {
   }

   PrRaXGParameters::~PrRaXGParameters()
   {
   }

   std::vector<std::string> PrRaXGParameters::names() const
   {
      // Storage for the names
      std::vector<std::string> names;

      // Setup the names of the 4 parameters
         // Prandtl number
      names.push_back("prandtl");
         // Rayleigh number
      names.push_back("rayleigh");
         // Slope angle chi
      names.push_back("chi");
         // Topographic ratio
      names.push_back("gamma");

      return names;
   }

   void PrRaXGParameters::init(const std::map<std::string,MHDFloat>& parameters)
   {
      // Storage for the name to ID mapping
      std::map<std::string, NonDimensional::Id> namesMap;

      // Setup the name to ID mapping
         // Prandtl number
      namesMap.insert(std::make_pair("prandtl", NonDimensional::PRANDTL));
         // Rayleigh number
      namesMap.insert(std::make_pair("rayleigh", NonDimensional::RAYLEIGH));
         // Slope angle chi 
      namesMap.insert(std::make_pair("chi", NonDimensional::CHI));
         // Topographic ratio
      namesMap.insert(std::make_pair("gamma", NonDimensional::GAMMA));

      // Add non dimensional parameters
      std::map<std::string,MHDFloat>::const_iterator it;
      for(it = parameters.begin(); it != parameters.end(); it++)
      {
         this->mND.insert(std::make_pair(namesMap.at(it->first), it->second));
      }
   }

}
