/** 
 * @file CflWriter.hpp
 * @brief Implementation of a CFL writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CFLWRITER_HPP
#define CFLWRITER_HPP

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "IoAscii/IAsciiEWriter.hpp"

namespace GeoMHDiSCC {

namespace IoAscii {

   /**
    * @brief Implementation of a CFL writer
    */
   class CflWriter: public IoAscii::IAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          */
         CflWriter();

         /**
          * @brief Destructor
          */
         virtual ~CflWriter();

         /**
          * @brief Set the simulation time parameters
          *
          * @param time       Reached simulation time
          * @param timestep   Last timestep size
          */
         void setSimTime(const MHDFloat time, const MHDFloat timestep);

         /**
          * @brief Write State to file
          */
         virtual void write();

      protected:
         /**
          * @brief Time
          */
         MHDFloat mTime;

         /**
          * @brief Timestep
          */
         MHDFloat mTimestep;

         /**
          * @brief Flag for timestep change
          */
         bool mChanged;

      private:
   };

   /// Typedef for a smart reference counting pointer of a CflWriter writer
   typedef SharedPtrMacro<CflWriter>   SharedCflWriter;

}
}

#endif // CFLWRITER_HPP
