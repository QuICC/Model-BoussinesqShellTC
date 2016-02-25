/** 
 * @file ImExRK2CN.hpp
 * @brief Implementation of an implicit/explicit Runge-Kutta 2nd order with Crank-Nicolson scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IMEXRK2CN_HPP
#define IMEXRK2CN_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * @brief Implementation of an implicit/explicit Runge-Kutta 2nd ordre with Crank-Nicolson scheme
    */
   class ImExRK2CN
   {
      public:
         /**
          * @brief Butcher's tableau a_ij factor for implicit scheme
          */
         static MHDFloat aIm(const int i, const int j);

         /**
          * @brief Butcher's tableau b_i factor for implicit scheme
          */
         static MHDFloat bIm(const int i);

         /**
          * @brief Butcher's tableau a_ij factor for explicit scheme
          */
         static MHDFloat aEx(const int i, const int j);

         /**
          * @brief Butcher's tableau b_i factor for explicit scheme
          */
         static MHDFloat bEx(const int i);

         /**
          * @brief Butcher's tableau b_i factor for implicit embedded lower order scheme
          */
         static MHDFloat bImErr(const int i);

         /**
          * @brief Butcher's tableau b_i factor for explicit embedded lower order scheme
          */
         static MHDFloat bExErr(const int i);

         /**
          * @brief Number of substeps for final step (this is +1 compared to theoretical value due to implementation)
          */
         static const int STEPS; 

         /**
          * @brief Order of the scheme
          */
         static const int ORDER;

         /**
          * @brief Scheme has embedded lower order scheme?
          */
         static const bool HAS_EMBEDDED;

         /**
          * @brief Name of the scheme
          */
         static const std::string NAME;

         /**
          * @brief Initialize Butcher's tableau
          */
         static void init();
         
      protected:
         /**
          * @brief Storage for the implicit a factors
          */
         static Eigen::Array<MHDFloat,3,3> mAIm;

         /**
          * @brief Storage for the implicit b factors
          */
         static Eigen::Array<MHDFloat,3,1> mBIm;

         /**
          * @brief Storage for the explicit a factors
          */
         static Eigen::Array<MHDFloat,3,3> mAEx;

         /**
          * @brief Storage for the explicit b factors
          */
         static Eigen::Array<MHDFloat,3,1> mBEx;

         /**
          * @brief Storage for the implicit embedded scheme b factors
          */
         static Eigen::Array<MHDFloat,3,1> mBImErr;

         /**
          * @brief Storage for the explicit embedded scheme b factors
          */
         static Eigen::Array<MHDFloat,3,1> mBExErr;

      private:
         /**
          * @brief Constructor
          */
         ImExRK2CN();

         /**
          * @brief Destructor
          */
         ~ImExRK2CN();

   };
}
}

#endif // IMEXRK2CN_HPP
