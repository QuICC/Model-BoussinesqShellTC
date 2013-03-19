/** \file VisualizationFileTags.hpp
 *  \brief Definitions and names use by the visualisation file writers
 */

#ifndef VISUALIZATIONFILEDEFS_HPP
#define VISUALIZATIONFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * \brief Definitions and names use by the visualisation file writers
    */
   class VisualizationFileTags
   {
      public:
         /**
          * @brief HEADER part for visualisation file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for visualisation file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of visualisation file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of visualisation file
          */
         static const std::string   EXTENSION;

         /**
          * @brief Mesh description for visualisation file
          */
         static const std::string   MESH;

         /**
          * @brief Grid description for visualisation file
          */
         static const std::string   GRID;

      private:
         /**
         * @brief Empty destructor
         */
         VisualizationFileTags();

         /**
         * @brief Destructor
         */
         ~VisualizationFileTags();
   };
}
}

#endif // VISUALIZATIONFILEDEFS_HPP
