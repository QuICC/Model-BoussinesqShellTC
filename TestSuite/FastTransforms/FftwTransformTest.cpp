/** \file FftwTransformTest.cpp
 *  \brief Implementation of a test case for FftwTransform
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "FastTransforms/FftwTransform.hpp"
#include "FastTransforms/FftwTools.hpp"

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
         
         /// Acceptable error
         double mError;
   };

   FftwTransformTest::FftwTransformTest()
      : mError(1e-13)
   {
      // Initilise framework
      FrameworkMacro::init();

      // Set nCpu for serial run
      int nCpu = 1;

      // Set ID and nCpu in MPI case
      #ifdef GEOMHDISCC_MPI
         // Get MPI size
         int size;
         MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
      #endif //GEOMHDISCC_MPI

      // Setup framework
      FrameworkMacro::setup(nCpu);
   }

   FftwTransformTest::~FftwTransformTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

   /**
    * @brief Test mesh grid
    */
   TEST_F(FftwTransformTest, MeshGrid)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasMixedFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::MIXED));

      // Create FftwTransform
      Transform::FftwTransform fft;

      // Initialise FFT
      fft.init(spSetup);

      EXPECT_EQ(xN, fft.meshGrid().size());
   }

   /**
    * @brief Test accuracy for forward mixed transform
    */
   TEST_F(FftwTransformTest, ForwardMixedAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasMixedFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::MIXED));

      // Create FftwTransform
      Transform::FftwTransform fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get fft grid
      Array phi = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         phys.col(i) = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftwTransform::IntegratorType::INTG);

      // Check solution
      EXPECT_NEAR(spec(0,0).real(), 1.0, this->mError);
      EXPECT_NEAR(spec(0,0).imag(), 0.0, this->mError);
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         EXPECT_NEAR(spec(i,i).real(), 0.5, this->mError);
         EXPECT_NEAR(spec(i,i).imag(), -0.5, this->mError);
      }
   }

   /**
    * @brief Test accuracy for backward mixed transform
    */
   TEST_F(FftwTransformTest, BackwardMixedAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasMixedFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::MIXED));

      // Create FftwTransform
      Transform::FftwTransform fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data as 0.5-0.5i up to maxN (1 for n = 0)
      spec(0,0) = MHDComplex(1.0,0.0);
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         spec(i,i) = MHDComplex(0.5,-0.5);
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftwTransform::ProjectorType::PROJ);

      // Get phi grid
      Array phi = fft.meshGrid();

      // Check the solution
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         for(int j = 0; j < phi.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), std::cos(static_cast<MHDFloat>(i)*phi(j)) + std::sin(static_cast<MHDFloat>(i)*phi(j)), this->mError);
         }
      }
   }

   /**
    * @brief Test accuracy for backward mixed derivative transform
    */
   TEST_F(FftwTransformTest, BackwardMixedDiffAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasMixedFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::MIXED));

      // Create FftwTransform
      Transform::FftwTransform fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get phi grid
      Array phi = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         phys.col(i) = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
      }

      // Compute backward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftwTransform::IntegratorType::INTG);

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftwTransform::ProjectorType::DIFF);

      // Check the solution
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         for(int j = 0; j < phi.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), -static_cast<MHDFloat>(i)*std::sin(static_cast<MHDFloat>(i)*phi(j)) + static_cast<MHDFloat>(i)*std::cos(static_cast<MHDFloat>(i)*phi(j)), this->mError);
         }
      }
   }

   /**
    * @brief Test accuracy for forward transform
    */
   TEST_F(FftwTransformTest, ForwardAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::EQUAL));

      // Create FftwTransform
      Transform::FftwTransform fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      MatrixZ  phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get fft grid
      Array phi = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         phys.col(i).real() = (static_cast<MHDFloat>(i)*phi).array().cos() - (static_cast<MHDFloat>(i)*phi).array().sin();
         phys.col(i).imag() = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
      }
      phys.col(0).imag().setZero();

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftwTransform::IntegratorType::INTG);

      // Check solution
      EXPECT_NEAR(spec(0,0).real(), 1.0, this->mError);
      EXPECT_NEAR(spec(0,0).imag(), 0.0, this->mError);
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         EXPECT_NEAR(spec(i,i).real(), 1.0, this->mError);
         EXPECT_NEAR(spec(i,i).imag(), 1.0, this->mError);
      }
   }

   /**
    * @brief Test accuracy for backward transform
    */
   TEST_F(FftwTransformTest, BackwardAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::EQUAL));

      // Create FftwTransform
      Transform::FftwTransform fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      MatrixZ  phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data as 0.5-0.5i up to maxN (1 for n = 0)
      spec.setZero();
      spec(0,0) = MHDComplex(1.0,0.0);
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         spec(i,i) = MHDComplex(1.0,1.0);
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftwTransform::ProjectorType::PROJ);

      // Get phi grid
      Array phi = fft.meshGrid();

      // Check the solution
      for(int j = 0; j < phi.size(); ++j)
      {
         EXPECT_NEAR(phys(j,0).real(), 1.0, this->mError);
         EXPECT_NEAR(phys(j,0).imag(), 0.0, this->mError);
      }
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         for(int j = 0; j < phi.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i).real(), std::cos(static_cast<MHDFloat>(i)*phi(j)) - std::sin(static_cast<MHDFloat>(i)*phi(j)), this->mError);
            EXPECT_NEAR(phys(j,i).imag(), std::cos(static_cast<MHDFloat>(i)*phi(j)) + std::sin(static_cast<MHDFloat>(i)*phi(j)), this->mError);
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
