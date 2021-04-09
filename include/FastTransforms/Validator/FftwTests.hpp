/**
 * @file Fftw.hpp
 * @brief Validation tests for FFTW library
 */

#ifndef VALIDATOR_FFTWTESTS_HPP
#define VALIDATOR_FFTWTESTS_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Validator {

   /**
    * @brief Validation tests for FFTW 
    */
   class FftwTests
   {
      public:
         /**
          * @brief Perform validation test for R2C transform
          */
         static void validateR2C();
         
      protected:

      private:

         /**
          * @brief Empty constructor
          */
         FftwTests() = default;

         /**
          * @brief Empty Destructor
          */
         ~FftwTests() = default;

   };

}
}
}

#endif // VALIDATOR_FFTWTESTS_HPP
