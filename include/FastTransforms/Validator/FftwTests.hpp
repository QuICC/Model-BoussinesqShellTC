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

         /**
          * @brief Perform validation test for C2R transform
          */
         static void validateC2R();

         /**
          * @brief Perform validation test for forward DFT transform
          */
         static void validateDftFwd();

         /**
          * @brief Perform validation test for backward DFT transform
          */
         static void validateDftBwd();

         /**
          * @brief Perform validation test for real R2R type 01 transform
          */
         static void validateR2R01();

         /**
          * @brief Perform validation test for real R2R type 10 transform
          */
         static void validateR2R10();

         /**
          * @brief Perform validation test for real R2R type 11 transform
          */
         static void validateR2R11();

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
