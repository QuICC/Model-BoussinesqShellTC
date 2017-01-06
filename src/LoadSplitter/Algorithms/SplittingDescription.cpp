/** 
 * @file SplittingDescription.cpp
 * @brief Source of the base of the implementation of the load splitting algorithm description
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "LoadSplitter/Algorithms/SplittingDescription.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   SplittingDescription::SplittingDescription()
   {
      #ifdef QUICC_DEBUG
         IoXml::SharedVtpWriter pVtp(new IoXml::VtpWriter("Data_distribution_TRAB1D"));
         pVtp->init();
         this->vtpFiles.push_back(pVtp);
         pVtp = IoXml::SharedVtpWriter(new IoXml::VtpWriter("Data_distribution_TRAB2D"));
         pVtp->init();
         this->vtpFiles.push_back(pVtp);
         pVtp = IoXml::SharedVtpWriter(new IoXml::VtpWriter("Data_distribution_TRAB3D"));
         pVtp->init();
         this->vtpFiles.push_back(pVtp);
      #endif //QUICC_DEBUG
   }

   SplittingDescription::~SplittingDescription()
   {
   }

}
}
