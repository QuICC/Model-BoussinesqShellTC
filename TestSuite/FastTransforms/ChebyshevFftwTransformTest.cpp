/** \file ChebyshevFftwTransformTest.cpp
 *  \brief Implementation of a test case for ChebyshevFftwTransform
 */

#include "FastTransforms/ChebyshevFftwTransform.hpp"
#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ChebyshevFftwTransform implementation
    */
   class ChebyshevFftwTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ChebyshevFftwTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevFftwTransformTest();

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

   ChebyshevFftwTransformTest::ChebyshevFftwTransformTest()
      : mSize(45), mHowmany(25), mSpec(15), mError(1e-13)
   {
   }

   ChebyshevFftwTransformTest::~ChebyshevFftwTransformTest()
   {
   }

   /**
    * @brief Test mesh grid
    */
   TEST_F(ChebyshevFftwTransformTest, MeshGrid) {
      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      Transform::SharedFftSetup spSetup(new Transform::FftSetup(this->mSize, this->mHowmany, this->mSpec, true));
      fft.init(spSetup);

      EXPECT_EQ(this->mSize, fft.meshGrid().size());
   }

   /**
    * @brief Test accuracy
    */
   TEST_F(ChebyshevFftwTransformTest, ForwardRealAccuracy) {
      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      Transform::SharedFftSetup spSetup(new Transform::FftSetup(this->mSize, this->mHowmany, this->mSpec, true));
      fft.init(spSetup);

      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      Array phi = fft.meshGrid();

      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         phys.col(i) = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
      }
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevFftwTransform::IntegratorType::INTG);

      EXPECT_NEAR(spec(0,0), 0.0, this->mError);
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         EXPECT_NEAR(spec(i,i), 0.5, this->mError);
      }
   }

   /**
    * @brief Test accuracy
    */
   TEST_F(ChebyshevFftwTransformTest, BackwardRealAccuracy) {
      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      Transform::SharedFftSetup spSetup(new Transform::FftSetup(this->mSize, this->mHowmany, this->mSpec, true));
      fft.init(spSetup);

      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      spec(0,0) = 1.0;
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         spec(i,i) = 0.5;
      }
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

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
