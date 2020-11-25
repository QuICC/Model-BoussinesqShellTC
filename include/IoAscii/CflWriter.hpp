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
#include "IoAscii/IAsciiWriter.hpp"

namespace QuICC {

namespace IoAscii {

   /**
    * @brief Implementation of a CFL writer
    */
   class CflWriter: public IoAscii::IAsciiWriter
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
          * @param cfl        CFL timesteps
          * @param steps      Number of steps with previous timestep
          */
         void setSimTime(const MHDFloat time, const Matrix& cfl, const MHDFloat steps);

         /**
          * @brief Make header
          */

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
          * @brief Number of steps with same timestep
          */
         MHDFloat mSteps;

         /**
          * @brief Error
          */
         MHDFloat mError;

         /**
          * @brief Flag for timestep change
          */
         bool mChanged;

         /**
          * @brief Timestep details
          */
         Matrix mDt;

      private:
         /**
          * @brief Write fancy header
          */
         virtual void fancyHeader();

         /**
          * @brief High precision output
          */
         const int mcIoHigh;

         /**
          * @brief Low precision output
          */
         const int mcIoLow;

         /**
          * @brief Width of exponent
          */
         const int mcIoExpW;

         /**
          * @brief Width of fixed
          */
         const int mcIoFixW;

         /**
          * @brief Width of integer
          */
         const int mcIoIntW;

         /**
          * @brief Need fancy header?
          */
         bool mNeedFancy;
   };

   /// Typedef for a smart reference counting pointer of a CflWriter writer
   typedef SharedPtrMacro<CflWriter>   SharedCflWriter;

}
}

#endif // CFLWRITER_HPP
