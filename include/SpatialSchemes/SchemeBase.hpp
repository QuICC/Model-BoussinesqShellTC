/** \file SchemeBase.hpp
 *  \brief Base implementation of a cost based spatial schemes
 *
 *  \mhdBug Needs test
 */

#ifndef SCHEMEBASE_HPP
#define SCHEMEBASE_HPP

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
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Base implementation of a cost based spatial schemes
    */
   class SchemeBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         SchemeBase(const int dims);

         /**
          * @brief Destructor
          */
         virtual ~SchemeBase();

         /**
          * @brief Get load balancing weights
          */
         virtual Array loadWeights() = 0;

         /**
          * @brief Get memory related score weight
          */
         virtual double memoryScore(SharedResolution spRes) = 0;
         
      protected:
         /**
          * @brief Cost of the transforms
          *
          * A value smaller than 1 increases the negative effect of inbalance.
          * A value bigger than 1 reduces the negative effect of inbalance.
          */
         Array   mCosts;

         /**
          * @brief Scaling of the transforms
          *
          * A value smaller than 1 increases the negative effect of inbalance.
          * A value bigger than 1 reduces the negative effect of inbalance.
          */
         Array   mScalings;

         /**
          * @brief Memory footprint of the transforms
          *
          * A value smaller than 1 increases the negative effect of inbalance.
          * A value bigger than 1 reduces the negative effect of inbalance.
          */
         Array   mMemory;

         /**
          * @brief Storage for the list of loads
          */
         std::deque<int>   mLoadList;

         /**
          * @brief Storage for the list of optimal load
          */
         std::queue<int>   mOptimalLoad;

         /**
          * @brief Storage for summed load
          */
         std::vector<int>   mLoad;

         /**
          * @brief Storage for the regularized load
          */
         std::multimap<int, int>   mRegularLoad;

         /**
          * @brief Set transform costs
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setCosts(const int shift = 0) = 0;

         /**
          * @brief Set transform scalings
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setScalings(const int shift = 0) = 0;

         /**
          * @brief Set transform memory footprint
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setMemory(const int shift = 0) = 0;

         /**
          * @brief Reset all the load related variables
          */
         void resetLoad();

         /**
          * @brief Update the load sum for the regular load
          */
         void updateLoad(const int parts);

      private:
   };
}


#endif // SCHEMEBASE_HPP
