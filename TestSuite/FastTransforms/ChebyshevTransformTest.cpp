/** 
 * @file ChebyshevTransformTest.cpp
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
   class ChebyshevTransformTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ChebyshevTransformTest();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevTransformTest();

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

   ChebyshevTransformTest::ChebyshevTransformTest()
      : mMaxN(512), mHowmany(100), mError(1e-10)
   {
   }

   ChebyshevTransformTest::~ChebyshevTransformTest()
   {
   }

   /**
    * @brief Test the mesh grid 
    *
    * @param ChebyshevTransformTest   Test fixture ID
    * @param MeshGrid                    Test ID
    */
   TEST_F(ChebyshevTransformTest, MeshGrid)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

      // Initialise FFT
      fft.init(spSetup);

      EXPECT_EQ(xN, fft.meshGrid().size());
      EXPECT_LT(fft.meshGrid()(xN-1), 1.0);
      EXPECT_GE(fft.meshGrid()(0), -1.0);
   }

   /**
    * @brief Test transform of chebyshev polynomials
    *
    * @param ChebyshevTransformTest Test fixture ID
    * @param Polynomials            Test ID
    */
   TEST_F(ChebyshevTransformTest, Polynomials)
   {
      // Set spectral and physical sizes
      int nN = 8;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, 8, nN, Transform::FftSetup::REAL));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      ACoeff x = fft.meshGrid();

      phys.col(0).setConstant(1.0);
      phys.col(1) = x;
      phys.col(2) = 2.0*x.pow(2) - 1.0;
      phys.col(3) = 4.0*x.pow(3) - 3.0*x;
      phys.col(4) = 8.0*x.pow(4) - 8.0*x.pow(2) + 1.0;
      phys.col(5) = 16.0*x.pow(5) - 20.0*x.pow(3) + 5.0*x;
      phys.col(6) = 32.0*x.pow(6) - 48.0*x.pow(4) + 18.0*x.pow(2) - 1.0;
      phys.col(7) = 64.0*x.pow(7) - 112.0*x.pow(5) + 56.0*x.pow(3) - 7.0*x;

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

      std::cerr << spec << std::endl;

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            if(i == j)
            {
               EXPECT_NEAR(1.0, spec(i,j), this->mError);
            } else
            {
               EXPECT_NEAR(0, spec(i,j), this->mError);
            }
         }
      }

   }

   /**
    * @brief Accuracy test for real forward transform
    *
    * @param ChebyshevTransformTest   Test fixture ID
    * @param ForwardRealAccuracy          Test ID
    */
   TEST_F(ChebyshevTransformTest, ForwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

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
         if(i < spSetup->specSize())
         {
            phys.col(i) = (static_cast<MHDFloat>(i)*x.array().acos()).cos();
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

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
    * @brief Accuracy test for real backward transform
    *
    * @param ChebyshevTransformTest   Test fixture ID
    * @param BackwardRealAccuracy         Test ID
    */
   TEST_F(ChebyshevTransformTest, BackwardRealAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

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
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::PROJ);

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
    * @brief Accuracy test for real backward derivative transform
    *
    * @param ChebyshevTransformTest   Test fixture ID
    * @param ForwardRealAccuracy          Test ID
    */
   TEST_F(ChebyshevTransformTest, BackwardRealDiffAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

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
         if(i < spSetup->specSize())
         {
            phys.col(i) = (5-i) - x.array() + x.array().pow(2) + -1.3*x.array().pow(3) + i*x.array().pow(4) -3.0*x.array().pow(i);
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

      // Compute backward derivative transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::DIFF);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         if(j < spSetup->specSize())
         {
            Array cheb;
            if(j > 0)
            {
               cheb = - 1 + 2.0*x.array() + -3.9*x.array().pow(2) + 4.0*j*x.array().pow(3) -3.0*j*x.array().pow(j-1);
            } else
            {
               cheb = - 1 + 2.0*x.array().array() + -3.9*x.array().pow(2) + 4.0*j*x.array().pow(3);
            }

            for(int i = 0; i < x.size(); ++i)
            {
               MHDFloat eta = std::max(10.0, std::abs(cheb(i)));
               EXPECT_NEAR(cheb(i)/eta, phys(i,j)/eta, this->mError);
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
    * @param ChebyshevTransformTest   Test fixture ID
    * @param BackwardRealAccuracy         Test ID
    */
   TEST_F(ChebyshevTransformTest, BackwardRealLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create Chebyshev transform
      Transform::ChebyshevTransformType fft;

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
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

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
    * @param ChebyshevTransformTest Test fixture ID
    * @param ForwardComplexAccuracy     Test ID
    */
   TEST_F(ChebyshevTransformTest, ForwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

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
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

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
    * @param ChebyshevTransformTest   Test fixture ID
    * @param BackwardComplexAccuracy      Test ID
    */
   TEST_F(ChebyshevTransformTest, BackwardComplexAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

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
            if(i == 0 && j == 0 && j < spSetup->specSize())
            {
               spec(i,j).real() = 1.0;
            } else if(i == j && j < spSetup->specSize() && j < spSetup->specSize())
            {
               spec(i,j).real() = 0.5;
            } else
            {
               spec(i,j).real() = 0.0;
            }
            
            if(i == 0 && spSetup->specSize()-1-j == 0 && j < spSetup->specSize())
            {
               spec(i,j).imag() = 1.0;
            } else if(i == spSetup->specSize()-1-j && j > 0 && j < spSetup->specSize())
            {
               spec(i,j).imag() = 0.5;
            } else
            {
               spec(i,j).imag() = 0.0;
            }
         }
      }

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::PROJ);

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
    * @brief Accuracy test for a complex forward derivative transform 
    *
    * @param ChebyshevTransformTest Test fixture ID
    * @param ForwardComplexAccuracy     Test ID
    */
   TEST_F(ChebyshevTransformTest, ForwardComplexDiffAccuracy)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

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
         if(i < spSetup->specSize())
         {
            phys.col(i).real() = (5-i) - x.array() + x.array().pow(2) - 1.3*x.array().pow(3) + i*x.array().pow(4) -3.0*x.array().pow(i);
            phys.col(i).imag() = 3.0*x.array() + 0.2*x.array().pow(3) - 2.1*i*x.array().pow(5) + x.array().pow(6) + 1.7*x.array().pow(i);
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

      // Compute backward derivative transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::DIFF);

      // Check solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         if(j < spSetup->specSize())
         {
            Array chebReal;
            Array chebImag;
            if(j > 0)
            {
               chebReal = -1.0 + 2.0*x.array() - 3.9*x.array().pow(2) + 4.0*j*x.array().pow(3) - 3.0*j*x.array().pow(j-1);
               chebImag =  3.0 + 0.6*x.array().pow(2) - 10.5*j*x.array().pow(4) + 6.0*x.array().pow(5) + 1.7*j*x.array().pow(j-1);
            } else
            {
               chebReal = -1.0 + 2.0*x.array() - 3.9*x.array().pow(2) + 4.0*j*x.array().pow(3);
               chebImag =  3.0 + 0.6*x.array().pow(2) - 10.5*j*x.array().pow(4) + 6.0*x.array().pow(5);
            }

            for(int i = 0; i < x.size(); ++i)
            {
               MHDFloat etaReal = std::max(10.0, std::abs(chebReal(i)));
               MHDFloat etaImag = std::max(10.0, std::abs(chebImag(i)));
               EXPECT_NEAR(chebReal(i)/etaReal, phys(i,j).real()/etaReal, this->mError);
               EXPECT_NEAR(chebImag(i)/etaImag, phys(i,j).imag()/etaImag, this->mError);
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
    * @param ChebyshevTransformTest   Test fixture ID
    * @param BackwardComplexAccuracy      Test ID
    */
   TEST_F(ChebyshevTransformTest, BackwardComplexLoop)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

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
      spec.block(0,0,spSetup->specSize(),spSetup->howmany()).setConstant(MHDComplex(1.0,1.0));

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevTransformType::ProjectorType::PROJ);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevTransformType::IntegratorType::INTG);

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         for(int i = 0; i < spSetup->specSize(); ++i)
         {
            EXPECT_NEAR(spec(i,j).real(), 1.0, this->mError);
            EXPECT_NEAR(spec(i,j).imag(), 1.0, this->mError);
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
