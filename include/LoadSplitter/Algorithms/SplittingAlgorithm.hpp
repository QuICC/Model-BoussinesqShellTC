/** \file SplittingAlgorithm.hpp
 *  \brief Base of the implementation of the load splitting algorithms
 */

#ifndef SPLITTINGALGORITHM_HPP
#define SPLITTINGALGORITHM_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <list>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Resolutions/Resolution.hpp"
#include "SpatialSchemes/SpatialScheme.hpp"
#include "LoadSplitter/Algorithms/SplittingDescription.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Base of the implementation of the load splitting algorithms
    */
   class SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id      ID of the CPU
          * @param nCpu    Number of CPUs used
          * @param dim     Dimensions
          * @param algo    ID of the algorithm
          */
         SplittingAlgorithm(const int id, const int nCpu, const ArrayI& dim, const Splitting::Algorithms::Id algo);

         /**
          * @brief Destructor
          */
         virtual ~SplittingAlgorithm();

         /**
          * @brief Compute splitting and the corresponding score
          */
         std::pair<int, std::pair<SharedResolution, SplittingDescription> > scoreSplitting();

         /**
          * @brief Initialise the spatial scheme
          *
          * @param dim  Dimensions
          */
         template <typename TSchemeType> void initScheme(const ArrayI& dim);

         /**
          * @brief Select next set of splitting factors
          */
         bool useNextFactors();

         /**
          * @brief Check if factorisation is applicable to scheme
          */
         virtual bool applicable() const = 0;
         
      protected:
         /**
          * @brief ID of the algorithm
          */
         Splitting::Algorithms::Id mAlgo;

         /**
          * @brief ID of the grouper
          */
         Splitting::Groupers::Id mGrouper;

         /**
          * @brief Shared spatial scheme
          */
         SharedSpatialScheme  mspScheme;

         /**
          * @brief Split dimension i
          *
          * @param i    Split the ith dimension
          * @param id   ID of the CPU
          */
         virtual SharedTransformResolution splitDimension(const int i, const int id) = 0;

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual int computeScore(SharedResolution spResolution) = 0;

         /**
          * @brief Compute score related to communication structure
          *
          * @param spRes   Shared resolution object
          * @param details Details of the communication structure
          */
         double communicationScore(SharedResolution spRes, ArrayI& details);

         /**
          * @brief Compute score related to load balancing
          *
          * @param spRes   Shared resolution object
          * @param balance Details of the load balancing (on input they contain weights)
          */
         double balancingScore(SharedResolution spRes, Array& balance);

         /**
          * @brief Get id of core
          */
         int id() const;

         /**
          * @brief Get Number of CPUs
          */
         int nCpu() const;

         /**
          * @brief Dimensionality of the fields
          */
         int dims() const;

         /**
          * @brief Get \f$F_{i}\f$ factor in \f$N_{cpu} = \prod_{i} F_{i}\f$
          *
          * @param i index of the factor
          */
         int factor(const int i) const;

         /**
          * @brief Get \f$F_{i}\f$ factors in \f$N_{cpu} = \prod_{i} F_{i}\f$
          */
         const ArrayI& factors() const;

         /**
          * @brief Get maximum \f$F_{i}\f$ factor in \f$N_{cpu} = \prod_{i} F_{i}\f$
          */
         int maxFactor() const;

         /**
          * @brief Decompose \f$N_{cpu}\f$ into possible factors
          *
          * @param factors Number of factors in factorisation
          */
         void factoriseNCpu(const int factors);

         /**
          * @brief Filter out the unusable factors
          */
         void filterFactors();

         /**
          * @brief Test the splitting factors compatibility
          */
         bool confirmFactors() const;

         /**
          * @brief Convert ID to \f$F_{i}\f$ groupd ID
          *
          * @param i    ID of the factorisation group
          * @param id   ID of CPU
          */
         int groupId(const int i, const int id) const;

         /**
          * @brief Compute a simple balanced split of the elements with regard to the given number of parts
          *
          * @param n0      Output start index
          * @param nN      Output number of indexes
          * @param tot     Total number of indexes
          * @param parts   Number of parts to split total into
          * @param id      ID of the CPU
          */
         void balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id) const;

         /**
          * @brief Extract splitting from mapped indexes
          *
          * @param mapped  Mapped indexes to extract from
          * @param rIdx    Output storage to put the indexes in
          * @param id      ID of the CPU
          */
         void splitMapped(const std::multimap<int, int>& mapped, ArrayI &rIdx, const int id) const;

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
          * @brief Number of dimensions
          */
         int mDims;

         /**
          * @brief Storage for the \f$N_{cpu}\f$ factorisation factors
          */
         ArrayI mFactors;

         /**
          * @brief List of \f$N_{cpu}\f$ factorisation factors
          */
         std::list<int> mNCpuFactors;

         /**
          * @brief Storage for the simulation resolution
          */
         ArrayI mSimDim;
   };

   template <typename TSchemeType> void SplittingAlgorithm::initScheme(const ArrayI& dim)
   {
      // Create shared pointer
      this->mspScheme = SharedPtrMacro<TSchemeType>(new TSchemeType(dim));

      // Initialise scheme
      this->mspScheme->init();
   }

   /// Typedef for a shared pointer to a SplittingAlgorithm object
   typedef SharedPtrMacro<SplittingAlgorithm>   SharedSplittingAlgorithm;

}

#endif // SPLITTINGALGORITHM_HPP
