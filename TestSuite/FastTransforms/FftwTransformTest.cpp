/** \file FftwTransformTest.cpp
 *  \brief Implementation of a test case for FftwTransform
 */

#include "FastTransforms/FftwTransform.hpp"
#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the FftwTransform implementation
    */
   class FftwTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FftwTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~FftwTransformTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
         
         int mSize;
         int mHowmany;
         int mSpec;
         double mError;
   };

   FftwTransformTest::FftwTransformTest()
      : mSize(45), mHowmany(25), mSpec(15), mError(1e-13)
   {
   }

   FftwTransformTest::~FftwTransformTest()
   {
   }

   /**
    * @brief Test mesh grid
    */
   TEST_F(FftwTransformTest, MeshGrid) {
      // Create FftwTransform
      Transform::FftwTransform fft;

      Transform::SharedFftSetup spSetup(new Transform::FftSetup(this->mSize, this->mHowmany, this->mSpec, true));
      fft.init(spSetup);

      EXPECT_EQ(this->mSize, fft.meshGrid().size());
   }

   /**
    * @brief Test accuracy
    */
   TEST_F(FftwTransformTest, ForwardMixedAccuracy) {
      // Create FftwTransform
      Transform::FftwTransform fft;

      Transform::SharedFftSetup spSetup(new Transform::FftSetup(this->mSize, this->mHowmany, this->mSpec, true));
      fft.init(spSetup);

      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      Array phi = fft.meshGrid();

      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         phys.col(i) = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
      }
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftwTransform::IntegratorType::INTG);

      EXPECT_NEAR(spec(0,0).real(), 1.0, this->mError);
      EXPECT_NEAR(spec(0,0).imag(), 0.0, this->mError);
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         EXPECT_NEAR(spec(i,i).real(), 0.5, this->mError);
         EXPECT_NEAR(spec(i,i).imag(), -0.5, this->mError);
      }
   }

   /**
    * @brief Test accuracy
    */
   TEST_F(FftwTransformTest, BackwardMixedAccuracy) {
      // Create FftwTransform
      Transform::FftwTransform fft;

      Transform::SharedFftSetup spSetup(new Transform::FftSetup(this->mSize, this->mHowmany, this->mSpec, true));
      fft.init(spSetup);

      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      spec(0,0) = MHDComplex(1.0,0.0);
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         spec(i,i) = MHDComplex(0.5,-0.5);
      }
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftwTransform::ProjectorType::PROJ);

      Array phi = fft.meshGrid();

      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         for(int j = 0; j < phi.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), std::cos(static_cast<MHDFloat>(i)*phi(j)) + std::sin(static_cast<MHDFloat>(i)*phi(j)), this->mError);
         }
      }
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
