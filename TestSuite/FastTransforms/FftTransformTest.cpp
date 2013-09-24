/** 
 * @file FftTransformTest.cpp
 * @brief Implementation of test cases for FFT transforms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "TypeSelectors/FftSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Googletest fixture for the FftTransform implementation
    */
   class FftTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FftTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~FftTransformTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
         
         /**
          * @brief Spectral size of transform
          */
         int mMaxN;

         /**
          * @brief How many identical transforms to compute
          */
         int mHowmany;
         
         /**
          * @brief Acceptable absolute error
          */
         double mError;
         
         /**
          * @brief Acceptable relative error
          */
         double mRelError;
   };

   FftTransformTest::FftTransformTest()
      //: mMaxN(512), mHowmany(100), mError(1e-10), mRelError(1e-10)
      : mMaxN(15), mHowmany(24), mError(1e-10), mRelError(1e-10)
   {
   }

   FftTransformTest::~FftTransformTest()
   {
   }

   /**
    * @brief Test the mesh grid
    *
    * @param FftTransformTest   Test fixture ID
    * @param MeshGrid            Test ID
    */
   TEST_F(FftTransformTest, MeshGrid)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasMixedFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::MIXED));

      // Create FFT transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      EXPECT_EQ(xN, fft.meshGrid().size());
      EXPECT_LT(fft.meshGrid()(xN-1), 2.0*MathConstants::PI);
      EXPECT_GE(fft.meshGrid()(0), 0);
   }

   /**
    * @brief Accuracy test for mixed forward transform
    *
    * @param FftTransformTest       Test fixture ID
    * @param ForwardMixedAccuracy   Test ID
    */
   TEST_F(FftTransformTest, ForwardMixedAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasMixedFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::MIXED));

      // Create FFT transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get fft grid
      Array phi = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            phys.col(i) = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftTransformType::IntegratorType::INTG);

      // Check solution
      EXPECT_NEAR(1.0, spec(0,0).real(), this->mError);
      EXPECT_NEAR(0.0, spec(0,0).imag(), this->mError);
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 1; i < spSetup->specSize(); ++i)
         {
            if(i == j && j < spSetup->specSize())
            {
               EXPECT_NEAR(0.5, spec(i,i).real(), this->mError);
               EXPECT_NEAR(-0.5, spec(i,i).imag(), this->mError);
            } else
            {
               EXPECT_NEAR(0.0, spec(i,j).real(), this->mError);
               EXPECT_NEAR(0.0, spec(i,j).imag(), this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for mixed backward transform
    *
    * @param FftTransformTest       Test fixture ID
    * @param BackwardMixedAccuracy  Test ID
    */
   TEST_F(FftTransformTest, BackwardMixedAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasMixedFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::MIXED));

      // Create FFF transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data as 0.5-0.5i up to maxN (1 for n = 0)
      spec.setZero();
      spec(0,0) = MHDComplex(1.0,0.0);
      for(int i = 1; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            spec(i,i) = MHDComplex(0.5,-0.5);
         } 
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftTransformType::ProjectorType::PROJ);

      // Get phi grid
      Array phi = fft.meshGrid();

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         if(j < spSetup->specSize())
         {
            for(int i = 0; i < phi.size(); ++i)
            {
               EXPECT_NEAR(phys(i,j), std::cos(static_cast<MHDFloat>(j)*phi(i)) + std::sin(static_cast<MHDFloat>(j)*phi(i)), this->mError);
            }
         } else
         {
            for(int i = 0; i < phi.size(); ++i)
            {
               EXPECT_NEAR(phys(i,j), 0.0, this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for mixed backward derivative transform
    *
    * @param FftTransformTest         Test fixture ID
    * @param BackwardMixedDiffAccuracy Test ID
    */
   TEST_F(FftTransformTest, BackwardMixedDiffAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasMixedFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::MIXED));

      // Create FFT transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get phi grid
      Array phi = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            phys.col(i) = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
         }
      }

      // Compute backward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftTransformType::IntegratorType::INTG);

      // Compute derivative projection transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftTransformType::ProjectorType::DIFF);

      // Check the solution
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            for(int j = 0; j < phi.size(); ++j)
            {
               EXPECT_NEAR(-static_cast<MHDFloat>(i)*std::sin(static_cast<MHDFloat>(i)*phi(j)) + static_cast<MHDFloat>(i)*std::cos(static_cast<MHDFloat>(i)*phi(j)), phys(j,i), this->mError);
            }
         } else
         {
            for(int j = 0; j < phi.size(); ++j)
            {
               EXPECT_NEAR(0.0, phys(j,i), this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for backward-forward mixed transform loop
    *
    * @param FftTransformTest    Test fixture ID
    * @param ForwardAccuracy     Test ID
    */
   TEST_F(FftTransformTest, BackwardMixedLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasMixedFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::MIXED));

      // Create FFF transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data as 1-i up to maxN (1 for n = 0)
      spec.setConstant(0.0);
      spec.block(0,0,spSetup->specSize(),spSetup->howmany()).setConstant(MHDComplex(1.0,1.0));
      spec.row(0).imag().setConstant(0.0);

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftTransformType::IntegratorType::INTG);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i != 0)
            {
               EXPECT_NEAR(spec(i,j).real(), 1.0, this->mError);
               EXPECT_NEAR(spec(i,j).imag(), 1.0, this->mError);
            } else
            {
               EXPECT_NEAR(spec(i,j).real(), 1.0, this->mError);
               EXPECT_NEAR(spec(i,j).imag(), 0.0, this->mError);
            }
         }
         for(int i = spSetup->specSize(); i < spSetup->bwdSize(); ++i)
         {
            EXPECT_NEAR(spec(i,j).real(), 0.0, this->mError);
            EXPECT_NEAR(spec(i,j).imag(), 0.0, this->mError);
         }
      }
   }

   /**
    * @brief Accuracy test for forward transform
    *
    * @param FftTransformTest   Test fixture ID
    * @param ForwardAccuracy     Test ID
    */
   TEST_F(FftTransformTest, ForwardAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create FFT transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      MatrixZ  phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get fft grid
      Array phi = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) -/+ sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            phys.col(i).real() = (static_cast<MHDFloat>(i)*phi).array().cos() - (static_cast<MHDFloat>(i)*phi).array().sin();
            phys.col(i).imag() = (static_cast<MHDFloat>(i)*phi).array().cos() + (static_cast<MHDFloat>(i)*phi).array().sin();
         } else
         {
            phys.col(i).setZero();
         }
      }
      phys.col(0).imag().setZero();

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftTransformType::IntegratorType::INTG);

      // Check solution
      EXPECT_NEAR(1.0, spec(0,0).real(), this->mError);
      EXPECT_NEAR(0.0, spec(0,0).imag(), this->mError);
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 1; i < spSetup->specSize(); ++i)
         {
            if(i == j && j < spSetup->specSize())
            {
               EXPECT_NEAR(1.0, spec(i,i).real(), this->mError);
               EXPECT_NEAR(1.0, spec(i,i).imag(), this->mError);
            } else
            {
               EXPECT_NEAR(0.0, spec(i,j).real(), this->mError);
               EXPECT_NEAR(0.0, spec(i,j).imag(), this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for backward transform
    *
    * @param FftTransformTest   Test fixture ID
    * @param BackwardAccuracy    Test ID
    */
   TEST_F(FftTransformTest, BackwardAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create FFT transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      MatrixZ  phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data as 1.0-1.0i up to maxN (1 for n = 0)
      spec.setZero();
      spec(0,0) = MHDComplex(1.0,0.0);
      for(int i = 1; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            spec(i,i) = MHDComplex(1.0,1.0);
         }
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftTransformType::ProjectorType::PROJ);

      // Get phi grid
      Array phi = fft.meshGrid();

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         if(j < spSetup->specSize())
         {
            for(int i = 0; i < phi.size(); ++i)
            {
               EXPECT_NEAR(phys(i,j).real(), std::cos(static_cast<MHDFloat>(j)*phi(i)) - std::sin(static_cast<MHDFloat>(j)*phi(i)), this->mError);
               EXPECT_NEAR(phys(i,j).imag(), std::cos(static_cast<MHDFloat>(j)*phi(i)) + std::sin(static_cast<MHDFloat>(j)*phi(i)), this->mError);
            }
         } else
         {
            for(int i = 0; i < phi.size(); ++i)
            {
               EXPECT_NEAR(phys(i,j).real(), 0.0, this->mError);
               EXPECT_NEAR(phys(i,j).imag(), 0.0, this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for backward-forward transform loop
    *
    * @param FftTransformTest   Test fixture ID
    * @param ForwardAccuracy     Test ID
    */
   TEST_F(FftTransformTest, BackwardLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create FFT transform
      Transform::FftTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Create test data storage
      MatrixZ  phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ  spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get fft grid
      Array phi = fft.meshGrid();

      // Initialise the spectral test data as 1-i up to maxN (1 for n = 0)
      spec.setZero();
      spec.block(0,0,spSetup->specSize(),spSetup->howmany()).setConstant(MHDComplex(1.0,1.0));

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::FftTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::FftTransformType::IntegratorType::INTG);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i != 0)
            {
               EXPECT_NEAR(spec(i,j).real(), 1.0, this->mError);
               EXPECT_NEAR(spec(i,j).imag(), 1.0, this->mError);
            } else
            {
               EXPECT_NEAR(spec(i,j).real(), 1.0, this->mError);
               EXPECT_NEAR(spec(i,j).imag(), 0.0, this->mError);
            }
         }
         for(int i = spSetup->specSize(); i < spSetup->bwdSize(); ++i)
         {
            EXPECT_NEAR(spec(i,j).real(), 0.0, this->mError);
            EXPECT_NEAR(spec(i,j).imag(), 0.0, this->mError);
         }
      }
   }

}
}

/**
 * @brief Main function to execute all test cases
 *
 * @param argc Number of arguments
 * @param argv Arguments
 */
int main(int argc, char **argv)
{
   // Initilise framework
   GeoMHDiSCC::FrameworkMacro::init();

   // Set nCpu for serial run
   int nCpu = 1;

   // Set ID and nCpu in MPI case
   #ifdef GEOMHDISCC_MPI
      // Get MPI size
      MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
   #endif //GEOMHDISCC_MPI

   // Setup framework
   GeoMHDiSCC::FrameworkMacro::setup(nCpu);

   ::testing::InitGoogleTest(&argc, argv);
   int status = RUN_ALL_TESTS();

   // Finalise framework
   GeoMHDiSCC::FrameworkMacro::finalize();

   return status;
}
