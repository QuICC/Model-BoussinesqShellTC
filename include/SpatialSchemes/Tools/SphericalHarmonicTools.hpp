/** \file SphericalHarmonicTools.hpp
 *  \brief Implementation of the tools for the spherical harmonics based schemes
 */

#ifndef SPHERICALHARMONICSTOOLS_HPP
#define SPHERICALHARMONICSTOOLS_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <deque>
#include <queue>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of the tools for the spherical harmonics based schemes
    */
   class SphericalHarmonicTools
   {
      public:
         static int nM(const int l, const int nM);

         static int nHarmonics(const int nL, const int nM);

         static void buildLMap(std::multimap<int,int>& harmonics, const int nL, const int nM);

         static void buildLHMap(std::multimap<int,int>& harmonics, const int nL, const int nM, const int h0, const int nH);

         static void initMLLoad(std::deque<int>& list, std::vector<int>& load, std::queue<int>& optimal, const int nL, const int nM, const int bins);

         static void combineMPairs(std::deque<int>& list, std::multimap<int, int>& regular, const int nM, const int bins);

         static void fillMPairsBins(std::deque<int>& list, std::multimap<int, int>& regular, std::vector<int>& load, const int bins);

         static void fillMRestBins(std::deque<int>& list, std::multimap<int, int>& regular, std::vector<int>& load, std::queue<int>& optimal, const int bins);

         static void convertLoadToOrders(std::multimap<int, int>& regular, const int nL);
   };
}
}

#endif // SPHERICALHARMONICSTOOLS_HPP
