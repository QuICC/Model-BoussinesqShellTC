/** 
 * @file ChebyshevOperatorTest.cpp
 * @brief Implementation of test case for ChebyshevOperator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "SpectralOperators/ChebyshevOperator.hpp"
#include "FastTransforms/ChebyshevFftwTransform.hpp"
#include "TypeSelectors/FftSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ChebyshevOperator implementation
    */
   class ChebyshevOperatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ChebyshevOperatorTest();

         /**
          * @brief Destructor
          */
         virtual ~ChebyshevOperatorTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
         
         /**
          * @brief Acceptable absolute error
          */
         double mError;
         
         /**
          * @brief Spectral size of transform
          */
         int mMaxN;

         /**
          * @brief How many identical transforms to compute
          */
         int mHowmany;
   };

   ChebyshevOperatorTest::ChebyshevOperatorTest()
      : mError(1e-10), mMaxN(31), mHowmany(10)
   {
   }

   ChebyshevOperatorTest::~ChebyshevOperatorTest()
   {
   }

//   void ChebyshevOperatorTest::SetUp()
//   {
//   }

//   void ChebyshevOperatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test the identity and pseudo identity
    */
   TEST_F(ChebyshevOperatorTest, Identity)
   {
      int nN = this->mMaxN + 1;
      Spectral::ChebyshevOperator   spec(nN);

      // Identity
      SparseMatrix qId = spec.id();
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i==j)
            {
               EXPECT_EQ(qId.coeff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qId.coeff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }

      // Positive q
      for(int q = 0; q < 5; ++q)
      {
         qId = spec.id(q);
         for(int j = 0; j < nN; ++j)
         {
            for(int i = 0; i < nN; ++i)
            {
               if(i > q-1 && i==j)
               {
                  EXPECT_EQ(qId.coeff(i,j), 1.0) << "i = " << i << " j = " << j;
               } else
               {
                  EXPECT_EQ(qId.coeff(i,j), 0.0) << "i = " << i << " j = " << j;
               }
            }
         }
      }

      // Negative q
      for(int q = -4; q < 0; ++q)
      {
         qId = spec.id(q);
         for(int j = 0; j < nN; ++j)
         {
            for(int i = 0; i < nN; ++i)
            {
               if(i < (nN + q) && i==j)
               {
                  EXPECT_EQ(qId.coeff(i,j), 1.0) << "i = " << i << " j = " << j;
               } else
               {
                  EXPECT_EQ(qId.coeff(i,j), 0.0) << "i = " << i << " j = " << j;
               }
            }
         }
      }
   }

   /**
    * @brief Test the shifted identity
    */
   TEST_F(ChebyshevOperatorTest, ShiftedIdentity)
   {
      ASSERT_TRUE(false) << "Test is not implemented yet!";
   }

   /**
    * @brief Test \f$D^{1}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, D1)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::Fft::ToolsSelector::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   tmpSpec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->specSize(), spSetup->howmany());

      // Get chebyshev grid
      ACoeff x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         MHDFloat di = static_cast<MHDFloat>(i);

         if(i < spSetup->specSize())
         {
            phys.col(i) = (5.0-di) - x + x.pow(2) - 1.3*x.pow(3) + di*x.pow(4) -3.0*x.pow(i);
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(tmpSpec, phys, Transform::ChebyshevFftwTransform::IntegratorType::INTG);

      // Dealias
      spec = tmpSpec.topRows(spSetup->specSize());

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   op(spSetup->specSize());

      // Compute first derivative
      spec = op.diff(0,1)*spec;

      // Compute backward transform
      tmpSpec.topRows(spSetup->specSize()) = spec;
      tmpSpec.bottomRows(spSetup->padSize()).setZero();
      fft.project<Arithmetics::SET>(phys, tmpSpec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         if(j < spSetup->specSize())
         {
            Array cheb = -1.0 + 2.0*x - 3.9*x.pow(2) + 4.0*dj*x.pow(3) - 3.0*dj*x.pow(j-1);

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
    * @brief Test \f$D^{2}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, D2)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::Fft::ToolsSelector::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   tmpSpec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->specSize(), spSetup->howmany());

      // Get chebyshev grid
      ACoeff x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         MHDFloat di = static_cast<MHDFloat>(i);

         if(i < spSetup->specSize())
         {
            phys.col(i) = (5.0-di) - x + x.pow(2) - 1.3*x.pow(3) + di*x.pow(4) -3.0*x.pow(i);
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(tmpSpec, phys, Transform::ChebyshevFftwTransform::IntegratorType::INTG);

      // Dealias
      spec = tmpSpec.topRows(spSetup->specSize());

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   op(spSetup->specSize());

      // Compute second derivative
      spec = op.diff(0,2)*spec;

      // Compute backward transform
      tmpSpec.topRows(spSetup->specSize()) = spec;
      tmpSpec.bottomRows(spSetup->padSize()).setZero();
      fft.project<Arithmetics::SET>(phys, tmpSpec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         if(j < spSetup->specSize())
         {
            Array cheb = 2.0 - 7.8*x + 12.0*dj*x.pow(2) - 3.0*dj*(dj-1.0)*x.pow(j-2);

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
    * @brief Test \f$D^{3}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, D3)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::Fft::ToolsSelector::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   tmpSpec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->specSize(), spSetup->howmany());

      // Get chebyshev grid
      ACoeff x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         MHDFloat di = static_cast<MHDFloat>(i);

         if(i < spSetup->specSize())
         {
            phys.col(i) = (5.0-di) - x + x.pow(2) - 1.3*x.pow(3) + di*x.pow(4) -3.0*x.pow(i);
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(tmpSpec, phys, Transform::ChebyshevFftwTransform::IntegratorType::INTG);

      // Dealias
      spec = tmpSpec.topRows(spSetup->specSize());

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   op(spSetup->specSize());

      // Compute second derivative
      spec = op.diff(0,3)*spec;

      std::cerr << op.diff(0,3) << std::endl;

      // Compute backward transform
      tmpSpec.topRows(spSetup->specSize()) = spec;
      tmpSpec.bottomRows(spSetup->padSize()).setZero();
      fft.project<Arithmetics::SET>(phys, tmpSpec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);

         if(j < spSetup->specSize())
         {
            Array cheb = - 7.8 + 24.0*dj*x - 3.0*dj*(dj-1.0)*(dj-2.0)*x.pow(j-3);

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
    * @brief Test \f$D^{4}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, D4)
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::Fft::ToolsSelector::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   tmpSpec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->specSize(), spSetup->howmany());

      // Get chebyshev grid
      ACoeff x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      for(int i = 0; i < spSetup->howmany(); ++i)
      {
         MHDFloat di = static_cast<MHDFloat>(i);

         if(i < spSetup->specSize())
         {
            phys.col(i) = (5.0-di) - x + x.pow(2) - 1.3*x.pow(3) + di*x.pow(4) -3.0*x.pow(i);
         } else
         {
            phys.col(i).setZero();
         }
      }

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(tmpSpec, phys, Transform::ChebyshevFftwTransform::IntegratorType::INTG);

      // Dealias
      spec = tmpSpec.topRows(spSetup->specSize());

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   op(spSetup->specSize());

      // Compute second derivative
      spec = op.diff(0,4)*spec;

      // Compute backward transform
      tmpSpec.topRows(spSetup->specSize()) = spec;
      tmpSpec.bottomRows(spSetup->padSize()).setZero();
      fft.project<Arithmetics::SET>(phys, tmpSpec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int j = 0; j < spSetup->howmany(); ++j)
      {
         MHDFloat dj = static_cast<MHDFloat>(j);
         if(j < spSetup->specSize())
         {
            Array cheb = 24.0*dj - 3.0*dj*(dj-1.0)*(dj-2.0)*(dj-3.0)*x.pow(j-4);

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
    * @brief Test \f$D^{-1}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q1)
   {
      ASSERT_TRUE(false) << "Test is not implemented yet!";
   }

   /**
    * @brief Test \f$D^{-2}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q2)
   {
      ASSERT_TRUE(false) << "Test is not implemented yet!";
   }

   /**
    * @brief Test \f$D^{-4}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q4)
   {
      ASSERT_TRUE(false) << "Test is not implemented yet!";
   }

   /**
    * @brief Test \f$D^{-1}D^{1}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q1D1)
   {
      int nN = this->mMaxN + 1;
      Spectral::ChebyshevOperator   spec(nN);

      // Compute product explicitly
      SparseMatrix qId = spec.qDiff(1,0)*spec.diff(0,1);

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 0 && i == j)
            {
               EXPECT_NEAR(qId.coeff(i,j), 1.0, this->mError) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_NEAR(qId.coeff(i,j), 0.0, this->mError) << "i = " << i << " j = " << j;
            }
         }
      }

      // Get implicit product (no computation involved)
      qId = spec.qDiff(1,1);

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 0 && i==j)
            {
               EXPECT_EQ(qId.coeff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qId.coeff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test \f$D^{-2}D^{2}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q2D2)
   {
      int nN = this->mMaxN + 1;
      Spectral::ChebyshevOperator   spec(nN);

      // Compute product explicitly
      SparseMatrix qId = spec.qDiff(2,0)*spec.diff(0,2);

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 1 && i == j)
            {
               EXPECT_NEAR(qId.coeff(i,j), 1.0, this->mError) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_NEAR(qId.coeff(i,j), 0.0, this->mError) << "i = " << i << " j = " << j;
            }
         }
      }

      // Get implicit product (no computation involved)
      qId = spec.qDiff(2,2);

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 1 && i==j)
            {
               EXPECT_EQ(qId.coeff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qId.coeff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test \f$D^{-4}D^{4}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q4D4)
   {
      // Set size of the basis
      int nN = this->mMaxN + 1;
      Spectral::ChebyshevOperator   spec(nN);

      // Compute product explicitly
      SparseMatrix qId = spec.qDiff(4,0)*spec.diff(0,4);

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 3 && i == j)
            {
               EXPECT_NEAR(qId.coeff(i,j), 1.0, this->mError) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_NEAR(qId.coeff(i,j), 0.0, this->mError) << "i = " << i << " j = " << j;
            }
         }
      }

      // Get implicit product (no computation involved)
      qId = spec.qDiff(4,4);

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 3 && i==j)
            {
               EXPECT_EQ(qId.coeff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qId.coeff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test \f$D^{-2}D^{1}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q2D1)
   {
      // Set size of the basis
      int nN = this->mMaxN + 1;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      Matrix qExact = spec.qDiff(2,1);
      SparseMatrix tmp = spec.id(2)*spec.qDiff(2,0)*spec.id(-2)*spec.diff(0,1);
      Matrix qProd = tmp;

      // Compare exact and product
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            EXPECT_NEAR(qExact(i,j), qProd(i,j), this->mError) << "i = " << i << " j = " << j;
         }
      }

      // Check operator
      tmp = spec.id(2)*spec.qDiff(1,0);
      Matrix corr = tmp;
      for(int j = 0; j < nN-1; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            EXPECT_EQ(qExact(i,j), corr(i,j)) << "i = " << i << " j = " << j;
         }
      }
   }

   /**
    * @brief Test \f$D^{-4}D^{2}\f$ operator
    */
   TEST_F(ChebyshevOperatorTest, Q4D2)
   {
      // Set size of the basis
      int nN = this->mMaxN + 1;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);


      // Get quasi-inverse of order 1
      Matrix qExact = spec.qDiff(4,2);
      SparseMatrix tmp = spec.id(4)*spec.qDiff(4,0)*spec.id(-4)*spec.diff(0,2);
      Matrix qProd = tmp;

      // Compare exact and product
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            EXPECT_NEAR(qExact(i,j), qProd(i,j), this->mError) << "i = " << i << " j = " << j;
         }
      }

      // Check operator
      tmp = spec.id(4)*spec.qDiff(2,0);
      Matrix corr = tmp;
      for(int j = 0; j < nN-2; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            EXPECT_EQ(qExact(i,j), corr(i,j)) << "i = " << i << " j = " << j;
         }
      }
   }

}
}

/// Main to execute all test from test case
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
