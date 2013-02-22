/** \file SplittingDescription.hpp
 *  \brief Implementation basic description of the obtained splitting (for any algorithm)
 */

#ifndef SPLITTINGDESCRIPTION_HPP
#define SPLITTINGDESCRIPTION_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/Enums/Splittings.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation basic description of the obtained splitting (for any algorithm)
    */
   class SplittingDescription
   {
      public:
         /**
          * @brief Constructor
          *
          * @param algorithm  ID of the algorithm
          * @param grouper    ID of the transform grouper
          * @param dims       Number of dimensions
          * @param factors    CPU splitting factors
          * @param score      Score 
          */
         SplittingDescription(const Splittings::Algorithms::Id algorithm, const Splittings::Groupers::Id grouper, const int dims, const ArrayI& factors, const MHDFloat score);

         /**
          * @brief Destructor
          */
         virtual ~SplittingDescription();

         /**
          * @brief ID of the algorithm
          */
         Splittings::Algorithms::Id algorithm;

         /**
          * @brief ID of the grouper
          */
         Splittings::Groupers::Id grouper;

         /**
          * @brief Number of dimensions
          */
         int dims;

         /**
          * @brief Storage for the \f$N_{cpu}\f$ factorisation factors
          */
         ArrayI factors;

         /**
          * @brief Score of the splitting
          */
         MHDFloat score;
         
      protected:

      private:
   };

}

#endif // SPLITTINGDESCRIPTION_HPP
