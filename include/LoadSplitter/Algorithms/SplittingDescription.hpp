/** 
 * @file SplittingDescription.hpp
 * @brief Implementation basic description of the obtained splitting (for any algorithm)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Enums/Splitting.hpp"
#include "IoXml/VtpWriter.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

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
         SplittingDescription();

         /**
          * @brief Destructor
          */
         ~SplittingDescription();

         /**
          * @brief ID of the algorithm
          */
         Splitting::Algorithms::Id algorithm;

         /**
          * @brief ID of the grouper
          */
         Splitting::Groupers::Id grouper;

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

         /**
          * @brief Communication structure
          */
         std::vector<std::multimap<int,int> >   structure;

         #ifdef GEOMHDISCC_DEBUG
         /**
          * @brief Storage for data distribution visualization files
          */
         std::vector<IoXml::SharedVtpWriter> vtpFiles;
         #endif //GEOMHDISCC_DEBUG
         
      protected:

      private:
   };

}
}

#endif // SPLITTINGDESCRIPTION_HPP
