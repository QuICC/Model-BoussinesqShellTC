/** 
 * @file ChebyshevFftTransformTest.cpp
 * @brief Implementation of test cases for Chebyshev transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "TypeSelectors/FftSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the Chebyshev transform implementation
    */
   class ChebyshevFftTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ChebyshevFftTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevFftTransformTest();

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

   ChebyshevFftTransformTest::ChebyshevFftTransformTest()
      : mMaxN(512), mHowmany(100), mError(1e-10), mRelError(1e-10)
   {
   }

   ChebyshevFftTransformTest::~ChebyshevFftTransformTest()
   {
   }

   /**
    * @brief Test the mesh grid 
    *
    * @param ChebyshevFftTransformTest   Test fixture ID
    * @param MeshGrid                    Test ID
    */
   TEST_F(ChebyshevFftTransformTest, MeshGrid)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      EXPECT_EQ(xN, fft.meshGrid().size());
      EXPECT_LT(fft.meshGrid()(xN-1), 1.0);
      EXPECT_GE(fft.meshGrid()(0), -1.0);
   }

   /**
    * @brief Accuracy test for real forward transform
    *
    * @param ChebyshevFftTransformTest   Test fixture ID
    * @param ForwardRealAccuracy          Test ID
    */
   TEST_F(ChebyshevFftTransformTest, ForwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data as T_i(x)
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         phys.col(i) = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == 0 && j == 0)
            {
               EXPECT_NEAR(1.0, spec(i,j), this->mError);
            } else if(j == i)
            {
               EXPECT_NEAR(0.5, spec(i,j), this->mError);
            } else
            {
               EXPECT_NEAR(0, spec(i,j), this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for real backward transform
    *
    * @param ChebyshevFftTransformTest   Test fixture ID
    * @param BackwardRealAccuracy         Test ID
    */
   TEST_F(ChebyshevFftTransformTest, BackwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data
      spec(0,0) = 1.0;
      for(int i = 1; i < spSetup->howmany(); ++i)
      {
         spec(i,i) = 0.5;
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::PROJ);

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Check the solution
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         Array cheb = (static_cast<MHDFloat>(i)*x.array().acos()).cos();

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(cheb(j), phys(j,i), this->mError);
         }
      }
   }

   /**
    * @brief Accuracy test for a complex forward transform 
    *
    * @param ChebyshevFftTransformTest Test fixture ID
    * @param ForwardComplexAccuracy     Test ID
    */
   TEST_F(ChebyshevFftTransformTest, ForwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      MatrixZ   phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ   spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data as T_n(n*phi) + T_n(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         phys.col(i).real() = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
         phys.col(i).imag() = (static_cast<MHDFloat>(spSetup->specSize()-1-i)*x.array().acos()).cos();
      }
      phys.col(0).imag().setZero();

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == 0 && j == 0)
            {
               EXPECT_NEAR(1.0, spec(i,j).real(), this->mError);
            } else if(i == j)
            {
               EXPECT_NEAR(0.5, spec(i,j).real(), this->mError);
            } else
            {
               EXPECT_NEAR(0, spec(i,j).real(), this->mError);
            }

            if(i == 0 && spSetup->specSize()-1-j == 0)
            {
               EXPECT_NEAR(1.0, spec(i,j).imag(), this->mError);
            } else if(j == spSetup->specSize()-1-i && j > 0)
            {
               EXPECT_NEAR(0.5, spec(i,j).imag(), this->mError);
            } else
            {
               EXPECT_NEAR(0, spec(i,j).imag(), this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for a complex backward transform 
    *
    * @param ChebyshevFftTransformTest   Test fixture ID
    * @param BackwardComplexAccuracy      Test ID
    */
   TEST_F(ChebyshevFftTransformTest, BackwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

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
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::PROJ);

      // Get chebyshev grid
      Array x = fft.meshGrid();


      // Check the solution
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         Array chebReal = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
         Array chebImag = (static_cast<MHDFloat>(spSetup->specSize()-i-1)*x.array().acos()).cos();

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(chebReal(j), phys(j,i).real(), this->mError);
            if(i > 0)
            {
               EXPECT_NEAR(chebImag(j), phys(j,i).imag(), this->mError);
            } else
            {
               EXPECT_NEAR(0.0, phys(j,i).imag(), this->mError);
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
