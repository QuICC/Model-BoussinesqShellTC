/** \file LoadSplitter.hpp
 *  \brief Implementation of the workload splitter over the available CPUs
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
         virtual ~LoadSplitter() {};

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
          * @brief Check if more splittings are available
          */
         bool  moreSplittings() const;

         /**
          * @brief Get splitting information of the next splitting
          */
         std::pair<SharedResolution,SplittingDescription> nextSplitting();
         
      protected:
         /**
          * @brief Get id of the CPU
          */
         int id() const;

         /**
          * @brief Get Number of CPUs
          */
         int nCpu() const;

      private:
         /**
          * @brief Counter to skip splittings
          */
         int mSkip;

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

#endif // LOADSPLITTER_HPP
