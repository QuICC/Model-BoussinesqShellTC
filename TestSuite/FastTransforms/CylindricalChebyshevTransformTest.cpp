/** 
 * @file CylindricalChebyshevTransformTest.cpp
 * @brief Implementation of test cases for cylindrical Chebyshev transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "TypeSelectors/FftSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the cylindrical Chebyshev transform implementation
    */
   class CylindricalChebyshevTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         CylindricalChebyshevTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~CylindricalChebyshevTransformTest();

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

   CylindricalChebyshevTransformTest::CylindricalChebyshevTransformTest()
      : mMaxN(512), mHowmany(100), mError(1e-10), mRelError(1e-10)
   {
   }

   CylindricalChebyshevTransformTest::~CylindricalChebyshevTransformTest()
   {
   }

   /**
    * @brief Test the mesh grid 
    *
    * @param CylindricalChebyshevTransformTest Test fixture ID
    * @param MeshGrid                              Test ID
    */
   TEST_F(CylindricalChebyshevTransformTest, MeshGrid)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create cylindrical Chebyshev transform
      Transform::CylindricalChebyshevTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      // Check size
      EXPECT_EQ(xN, fft.meshGrid().size());
      EXPECT_LT(fft.meshGrid()(xN-1), 1.0);
      EXPECT_GE(fft.meshGrid()(0), 0.0);
   }

   /**
    * @brief Accuracy test for a real forward transform 
    *
    * @param CylindricalChebyshevTransformTest Test fixture ID
    * @param ForwardRealAccuracy                   Test ID
    */
   TEST_F(CylindricalChebyshevTransformTest, ForwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create cylindrical Chebyshev transform
      Transform::CylindricalChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            phys.col(i) = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::CylindricalChebyshevTransformType::IntegratorType::INTG);

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == 0 && j == 0 && j < spSetup->specSize())
            {
               EXPECT_NEAR(1.0, spec(i,j), this->mError);
            } else if(j == i && j < spSetup->specSize())
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
    * @brief Accuracy test for a real backward transform 
    *
    * @param CylindricalChebyshevTransformTest Test fixture ID
    * @param BackwardRealAccuracy                  Test ID
    */
   TEST_F(CylindricalChebyshevTransformTest, BackwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create cylindrical Chebyshev transform
      Transform::CylindricalChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data
      spec.setZero();
      spec(0,0) = 1.0;
      for(int i = 1; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            spec(i,i) = 0.5;
         }
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::CylindricalChebyshevTransformType::ProjectorType::PROJ);

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         if(j < spSetup->specSize())
         {
            Array cheb = (static_cast<MHDFloat>(j)*x.array().acos()).cos();

            for(int i = 0; i < x.size(); ++i)
            {
               EXPECT_NEAR(cheb(i), phys(i,j), this->mError);
            }
         } else
         {
            for(int i = 0; i < x.size(); ++i)
            {
               EXPECT_NEAR(phys(i,j), 0.0, this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for real backward-forward transform loop
    *
    * @param CylindricalChebyshevTransformTest   Test fixture ID
    * @param BackwardRealAccuracy         Test ID
    */
   TEST_F(CylindricalChebyshevTransformTest, BackwardRealLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::EQUAL));

      // Create Chebyshev transform
      Transform::CylindricalChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data
      spec.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            spec(i,i) = 1.0;
         }
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::CylindricalChebyshevTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::CylindricalChebyshevTransformType::IntegratorType::INTG);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         if(j < spSetup->specSize())
         {
            EXPECT_NEAR(spec(j,j), 1.0, this->mError);
         } else
         {
            EXPECT_NEAR(spec(j,j), 0.0, this->mError);
         }
      }
   }

   /**
    * @brief Accuracy test for a complex forward transform 
    *
    * @param CylindricalChebyshevTransformTest Test fixture ID
    * @param ForwardComplexAccuracy                Test ID
    */
   TEST_F(CylindricalChebyshevTransformTest, ForwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create cylindrical Chebyshev transform
      Transform::CylindricalChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      MatrixZ   phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ   spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         if(i < spSetup->specSize())
         {
            phys.col(i).real() = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
            phys.col(i).imag() = (static_cast<MHDFloat>(spSetup->specSize()-1-i)*x.array().acos()).cos();
         } else
         {
            phys.col(i).setZero();
         }
      }
      phys.col(0).imag().setZero();

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::CylindricalChebyshevTransformType::IntegratorType::INTG);

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == 0 && j == 0 && j < spSetup->specSize())
            {
               EXPECT_NEAR(1.0, spec(i,j).real(), this->mError);
            } else if(i == j && j < spSetup->specSize())
            {
               EXPECT_NEAR(0.5, spec(i,j).real(), this->mError);
            } else
            {
               EXPECT_NEAR(0, spec(i,j).real(), this->mError);
            }

            if(i == 0 && spSetup->specSize()-1-j == 0 && j < spSetup->specSize())
            {
               EXPECT_NEAR(1.0, spec(i,j).imag(), this->mError);
            } else if(j == spSetup->specSize()-1-i && j > 0 && j < spSetup->specSize())
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
    * @param CylindricalChebyshevTransformTest    Test fixture ID
    * @param BackwardComplexAccuracy                  Test ID
    */
   TEST_F(CylindricalChebyshevTransformTest, BackwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create cylindrical Chebyshev transform
      Transform::CylindricalChebyshevTransformType fft;

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
            if(i == 0 && j == 0 && j < spSetup->specSize())
            {
               spec(i,j) = MHDComplex(1.0,0.0);
            } else if(i == j && j < spSetup->specSize() && j < spSetup->specSize())
            {
               spec(i,j) = MHDComplex(0.5,0.0);
            } else
            {
               spec(i,j) = 0.0;
            }
            
            if(i == 0 && spSetup->specSize()-1-j == 0 && j < spSetup->specSize())
            {
               spec(i,j) = MHDComplex(0.0,1.0);
            } else if(i == spSetup->specSize()-1-j && j > 0 && j < spSetup->specSize())
            {
               spec(i,j) = MHDComplex(0.0,0.5);
            } else
            {
               spec(i,j) = 0.0;
            }
         }
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::CylindricalChebyshevTransformType::ProjectorType::PROJ);

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         if(j < spSetup->specSize())
         {
            Array chebReal = (static_cast<MHDFloat>(j)*x.array().acos()).cos();
            Array chebImag = (static_cast<MHDFloat>(spSetup->specSize()-j-1)*x.array().acos()).cos();

            for(int i = 0; i < x.size(); ++i)
            {
               EXPECT_NEAR(chebReal(i), phys(i,j).real(), this->mError);
               if(j > 0)
               {
                  EXPECT_NEAR(chebImag(i), phys(i,j).imag(), this->mError);
               } else
               {
                  EXPECT_NEAR(0.0, phys(i,j).imag(), this->mError);
               }
            }
         } else
         {
            for(int i = 0; i < x.size(); ++i)
            {
               EXPECT_NEAR(phys(i,j).real(), 0.0, this->mError);
               EXPECT_NEAR(phys(i,j).imag(), 0.0, this->mError);
            }
         }
      }
   }

   /**
    * @brief Accuracy test for a complex backward-forward transform  loop
    *
    * @param CylindricalChebyshevTransformTest   Test fixture ID
    * @param BackwardComplexAccuracy      Test ID
    */
   TEST_F(CylindricalChebyshevTransformTest, BackwardComplexLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create Chebyshev transform
      Transform::CylindricalChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      MatrixZ   phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ   spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data
      spec.setZero();
      spec.block(0,0,spSetup->specSize(),spSetup->howmany()).setConstant(MHDComplex(1.0,1.0));

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::CylindricalChebyshevTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::CylindricalChebyshevTransformType::IntegratorType::INTG);

      // Get chebyshev grid
      Array x = fft.meshGrid();

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
