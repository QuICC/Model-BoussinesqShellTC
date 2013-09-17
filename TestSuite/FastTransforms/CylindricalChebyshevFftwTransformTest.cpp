/** 
 * @file CylindricalChebyshevFftwTransformTest.cpp
 * @brief Implementation of test cases for cylindrical fast Chebyshev transform (FFTW backend)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "FastTransforms/CylindricalChebyshevFftwTransform.hpp"
#include "FastTransforms/FftwTools.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the CylindricalChebyshevFftwTransform implementation
    */
   class CylindricalChebyshevFftwTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         CylindricalChebyshevFftwTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~CylindricalChebyshevFftwTransformTest();

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

   CylindricalChebyshevFftwTransformTest::CylindricalChebyshevFftwTransformTest()
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

   CylindricalChebyshevFftwTransformTest::~CylindricalChebyshevFftwTransformTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

   /**
    * @brief Test the mesh grid 
    *
    * @param CylindricalChebyshevFftwTransformTest Test fixture ID
    * @param MeshGrid                              Test ID
    */
   TEST_F(CylindricalChebyshevFftwTransformTest, MeshGrid)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::EQUAL));

      // Create CylindricalChebyshevFftwTransform
      Transform::CylindricalChebyshevFftwTransform fft;

      // Initialise FFT
      fft.init(spSetup);

      // Check size
      EXPECT_EQ(xN, fft.meshGrid().size());
   }

   /**
    * @brief Accuracy test for a real forward transform 
    *
    * @param CylindricalChebyshevFftwTransformTest Test fixture ID
    * @param ForwardRealAccuracy                   Test ID
    */
   TEST_F(CylindricalChebyshevFftwTransformTest, ForwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::EQUAL));

      // Create CylindricalChebyshevFftwTransform
      Transform::CylindricalChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         phys.col(i) = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::CylindricalChebyshevFftwTransform::IntegratorType::INTG);

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == 0 && j == 0)
            {
               EXPECT_NEAR(spec(i,j), 1.0, this->mError);
            } else if(j == i)
            {
               EXPECT_NEAR(spec(i,j), 0.5, this->mError);
            } else
            {
               EXPECT_NEAR(spec(i,j), 0, this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for a real backward transform 
    *
    * @param CylindricalChebyshevFftwTransformTest Test fixture ID
    * @param BackwardRealAccuracy                  Test ID
    */
   TEST_F(CylindricalChebyshevFftwTransformTest, BackwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::EQUAL));

      // Create CylindricalChebyshevFftwTransform
      Transform::CylindricalChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data
      spec(0,0) = 1.0;
      for(int i = 1; i < spSetup->specSize(); ++i)
      {
         spec(i,i) = 0.5;
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::CylindricalChebyshevFftwTransform::ProjectorType::PROJ);

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Check the solution
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         Array cheb = (static_cast<MHDFloat>(i)*x.array().acos()).cos();

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), cheb(j), this->mError);
         }
      }
   }

   /**
    * @brief Accuracy test for a complex forward transform 
    *
    * @param CylindricalChebyshevFftwTransformTest Test fixture ID
    * @param ForwardComplexAccuracy                Test ID
    */
   TEST_F(CylindricalChebyshevFftwTransformTest, ForwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN+1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::COMPONENT));

      // Create CylindricalChebyshevFftwTransform
      Transform::CylindricalChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      MatrixZ   phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ   spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         phys.col(i).real() = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
         phys.col(i).imag() = (static_cast<MHDFloat>(spSetup->specSize()-1-i)*x.array().acos()).cos();
      }
      phys.col(0).imag().setZero();

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::CylindricalChebyshevFftwTransform::IntegratorType::INTG);

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == 0 && j == 0)
            {
               EXPECT_NEAR(spec(i,j).real(), 1.0, this->mError);
            } else if(i == j)
            {
               EXPECT_NEAR(spec(i,j).real(), 0.5, this->mError);
            } else
            {
               EXPECT_NEAR(spec(i,j).real(), 0, this->mError);
            }

            if(i == 0 && spSetup->specSize()-1-j == 0)
            {
               EXPECT_NEAR(spec(i,j).imag(), 1.0, this->mError);
            } else if(j == spSetup->specSize()-1-i && j > 0)
            {
               EXPECT_NEAR(spec(i,j).imag(), 0.5, this->mError);
            } else
            {
               EXPECT_NEAR(spec(i,j).imag(), 0, this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for a complex backward transform 
    *
    * @param CylindricalChebyshevFftwTransformTest    Test fixture ID
    * @param BackwardComplexAccuracy                  Test ID
    */
   TEST_F(CylindricalChebyshevFftwTransformTest, BackwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int maxN = 11;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 13;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::COMPONENT));

      // Create CylindricalChebyshevFftwTransform
      Transform::CylindricalChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      MatrixZ   phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ   spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data
      spec.setZero();
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == 0 && j == 0)
            {
               spec(i,j) = MHDComplex(1.0,0.0);
            } else if(i == j)
            {
               spec(i,j) = MHDComplex(0.5,0.0);
            }
            
            if(i == 0 && spSetup->specSize()-1-j == 0)
            {
               spec(i,j) = MHDComplex(0.0,1.0);
            } else if(i == spSetup->specSize()-1-j && j > 0)
            {
               spec(i,j) = MHDComplex(0.0,0.5);
            }
         }
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::CylindricalChebyshevFftwTransform::ProjectorType::PROJ);

      // Get chebyshev grid
      Array x = fft.meshGrid();


      // Check the solution
      for(int i = 0; i < spSetup->specSize(); ++i)
      {
         Array chebReal = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
         Array chebImag = (static_cast<MHDFloat>(spSetup->specSize()-i-1)*x.array().acos()).cos();

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i).real(), chebReal(j), this->mError);
            if(i > 0)
            {
               EXPECT_NEAR(phys(j,i).imag(), chebImag(j), this->mError);
            } else
            {
               EXPECT_NEAR(phys(j,i).imag(), 0.0, this->mError);
            }
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
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
