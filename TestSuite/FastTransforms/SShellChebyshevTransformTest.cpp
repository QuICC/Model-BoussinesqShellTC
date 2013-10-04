/** 
 * @file SShellChebyshevTransformTest.cpp
 * @brief Implementation of test cases for spherical Chebyshev transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "TypeSelectors/FftSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the spherical Chebyshev transform implementation
    */
   class SShellChebyshevTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SShellChebyshevTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~SShellChebyshevTransformTest();

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
   };

   SShellChebyshevTransformTest::SShellChebyshevTransformTest()
      : mMaxN(512), mHowmany(100), mError(1e-10)
   {
   }

   SShellChebyshevTransformTest::~SShellChebyshevTransformTest()
   {
   }

   /**
    * @brief Test the mesh grid 
    *
    * @param SShellChebyshevTransformTest   Test fixture ID
    * @param MeshGrid                              Test ID
    */
   TEST_F(SShellChebyshevTransformTest, MeshGrid)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create spherical Chebyshev transform
      Transform::SShellChebyshevTransformType fft;

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
    * @param SShellChebyshevTransformTest   Test fixture ID
    * @param ForwardRealAccuracy                   Test ID
    */
   TEST_F(SShellChebyshevTransformTest, ForwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create spherical Chebyshev transform
      Transform::SShellChebyshevTransformType fft;

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
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::SShellChebyshevTransformType::IntegratorType::INTG);

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
    * @param SShellChebyshevTransformTest   Test fixture ID
    * @param BackwardRealAccuracy                  Test ID
    */
   TEST_F(SShellChebyshevTransformTest, BackwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create spherical Chebyshev transform
      Transform::SShellChebyshevTransformType fft;

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
      fft.project<Arithmetics::SET>(phys, spec, Transform::SShellChebyshevTransformType::ProjectorType::PROJ);

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
    * @param SShellChebyshevTransformTest   Test fixture ID
    * @param BackwardRealAccuracy         Test ID
    */
   TEST_F(SShellChebyshevTransformTest, BackwardRealLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create Chebyshev transform
      Transform::SShellChebyshevTransformType fft;

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
      fft.project<Arithmetics::SET>(phys, spec, Transform::SShellChebyshevTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::SShellChebyshevTransformType::IntegratorType::INTG);

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
    * @param SShellChebyshevTransformTest
    * @param ForwardComplexAccuracy
    */
   TEST_F(SShellChebyshevTransformTest, ForwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create spherical Chebyshev transform
      Transform::SShellChebyshevTransformType fft;

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
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::SShellChebyshevTransformType::IntegratorType::INTG);

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
    * @param SShellChebyshevTransformTest   Test fixture ID
    * @param BackwardComplexAccuracy               Test ID
    */
   TEST_F(SShellChebyshevTransformTest, BackwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create spherical Chebyshev transform
      Transform::SShellChebyshevTransformType fft;

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
      fft.project<Arithmetics::SET>(phys, spec, Transform::SShellChebyshevTransformType::ProjectorType::PROJ);

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
    * @param SShellChebyshevTransformTest   Test fixture ID
    * @param BackwardComplexAccuracy      Test ID
    */
   TEST_F(SShellChebyshevTransformTest, BackwardComplexLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::COMPONENT));

      // Create Chebyshev transform
      Transform::SShellChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      MatrixZ   phys = MatrixZ::Zero(spSetup->fwdSize(), spSetup->howmany());
      MatrixZ   spec = MatrixZ::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Initialise the spectral test data
      spec.setZero();
      spec.block(0,0,spSetup->specSize(),spSetup->howmany()).setConstant(MHDComplex(1.0,1.0));

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::SShellChebyshevTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::SShellChebyshevTransformType::IntegratorType::INTG);

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
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
