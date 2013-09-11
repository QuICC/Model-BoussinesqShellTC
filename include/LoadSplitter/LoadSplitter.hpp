/** 
 * @file LoadSplitter.hpp
 * @brief Implementation of the workload splitter over the available CPUs
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef LOADSPLITTER_HPP
#define LOADSPLITTER_HPP

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "LoadSplitter/Algorithms/SplittingDescription.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * \brief Implementation of the workload splitter over the available CPUs
    */
   class LoadSplitter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id      ID of the CPU
          * @param nCpu    Number of CPUs used
          */
         LoadSplitter(const int id, const int nCpu);

         /**
          * @brief Destructor
          */
         ~LoadSplitter();

         /**
          * @brief Initialise the algorithms
          *
          * @param dim  Dimensions
          *
          * @tparam TSchemeType  Underlying spatial scheme
          */
         template <typename TSchemeType> void init(const ArrayI& dim);

         /**
          * @brief Get splitting information of the best splitting
          */
         std::pair<SharedResolution,SplittingDescription> bestSplitting() const;

         /**
          * @brief Show description of some splittings
          *
          * @param n Maximum number of splittings to show
          */
         void  showSplittings(const int n) const;
         
      protected:

      private:
         /**
          * @brief ID of CPU
          */
         int mId;

         /**
          * @brief Number of CPUs
          */
         int mNCpu;

         /**
          * @brief Storage for the actual splitting algorithms used
          */
         std::vector<SharedSplittingAlgorithm>  mAlgorithms;

         /**
          * @brief Storage for the scores of the splitting
          */
         std::multimap<int, std::pair<SharedResolution,SplittingDescription> >  mScores;

         /**
          * @brief Initialise the splitting algorithms
          *
          * @param dim     Dimensions (spectral)
          */
         void initAlgorithms(const ArrayI& dim);

         /**
          * @brief Initialise the scores and corresponding splitting
          */
         void initScores();

         /**
          * @brief Describe the obtained splitting
          *
          * @param descr Splitting description object
          */
         void describeSplitting(const SplittingDescription& descr) const;
   };

   template <typename TSchemeType> void LoadSplitter::init(const ArrayI& dim)
   {
      // Assert that scheme and dimensions are compatible
      assert(TSchemeType::DIMENSIONS == dim.size());

      // Initialise the splitting algorithms
      this->initAlgorithms(dim);

      // Iterator over splitting algorithms
      std::vector<SharedSplittingAlgorithm>::iterator it;

      // Loop over all initialised algorithms
      for(it = this->mAlgorithms.begin(); it != this->mAlgorithms.end(); it++)
      {
         (*it)->initScheme<TSchemeType>(dim);
      }

      // Compute the core resolutions and corresponding scores
      this->initScores();
   }
}
}

#endif // LOADSPLITTER_HPP
