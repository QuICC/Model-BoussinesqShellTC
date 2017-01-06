/** 
 * @file ImExRKCB3e.hpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3e (Cavaglieri & Bewley, 2015)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IMEXRKCB3E_HPP
#define IMEXRKCB3E_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3e
    */
   class ImExRKCB3e
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
          * @brief Butcher's tableau c_i factor for explicit scheme
          */
         static MHDFloat cEx(const int i);

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
          * @brief Use scheme's embedded lower order scheme?
          */
         static bool USE_EMBEDDED;

         /**
          * @brief Name of the scheme
          */
         static const std::string NAME;

         /**
          * @brief Initialize Butcher's tableau
          */
         static void init();

         /**
          * @brief Activate embedded scheme
          */
         static void useEmbedded();
         
      protected:
         /**
          * @brief Storage for the implicit a factors
          */
         static Eigen::Array<MHDFloat,4,4> mAIm;

         /**
          * @brief Storage for the implicit b factors
          */
         static Eigen::Array<MHDFloat,4,1> mBIm;

         /**
          * @brief Storage for the explicit a factors
          */
         static Eigen::Array<MHDFloat,4,4> mAEx;

         /**
          * @brief Storage for the explicit b factors
          */
         static Eigen::Array<MHDFloat,4,1> mBEx;

         /**
          * @brief Storage for the explicit c factors
          */
         static Eigen::Array<MHDFloat,4,1> mCEx;

         /**
          * @brief Storage for the implicit embedded scheme b factors
          */
         static Eigen::Array<MHDFloat,4,1> mBImErr;

         /**
          * @brief Storage for the explicit embedded scheme b factors
          */
         static Eigen::Array<MHDFloat,4,1> mBExErr;

      private:
         /**
          * @brief Constructor
          */
         ImExRKCB3e();

         /**
          * @brief Destructor
          */
         ~ImExRKCB3e();

   };
}
}

#endif // IMEXRKCB3E_HPP
