/** \file VisualizeResolution.hpp
 *  \brief Implementation of a few useful tools to output resolution and load splitting
 *
 *  \mhdBug Implementation needs to be checked
 */

#ifndef VISUALIZERESOLUTION_HPP
#define VISUALIZERESOLUTION_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <iostream>

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace IoTools {

   /**
    * \brief Implementation of a few useful output formating static functions
    */
   class VisualizeResolution
   {
      public:
         /**
          * @brief Show the resolution in usable format
          */
         static void show(std::ostream& stream, std::string graphName, SharedResolution spRes);
         
      protected:

      private:
         /**
          * @brief Show the 1D resolution in usable format
          */
         static void show1D(std::ostream& stream, std::string graphName, SharedResolution spRes);

         /**
          * @brief Show the 2D resolution in usable format
          */
         static void show2D(std::ostream& stream, std::string graphName, SharedResolution spRes);
         
         /**
          * @brief Show the 3D resolution in usable format
          */
         static void show3D(std::ostream& stream, std::string graphName, SharedResolution spRes);
         
         /**
         * @brief Constructor 
         */
         VisualizeResolution();

         /**
         * @brief Destructor
         */
         ~VisualizeResolution();
   };

}
}

#endif // VISUALIZERESOLUTION_HPP
