/** 
 * @file ChebyshevOperatorTest.cpp
 * @brief Implementation of test case for ChebyshevOperator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "SpectralOperators/ChebyshevOperator.hpp"
#include "FastTransforms/ChebyshevFftwTransform.hpp"
#include "FastTransforms/FftwTools.hpp"

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
   };

   ChebyshevOperatorTest::ChebyshevOperatorTest()
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

   ChebyshevOperatorTest::~ChebyshevOperatorTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void ChebyshevOperatorTest::SetUp()
//   {
//   }

//   void ChebyshevOperatorTest::TearDown()
//   {
//   }

   /**
    * @brief Test c constant in by computing first derivative (fft transform loop)
    */
   TEST_F(ChebyshevOperatorTest, Diff1)
   {
      // Set spectral and physical sizes
      int maxN = 6;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 5;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::EQUAL));

      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      phys.col(0).setConstant(1.0);
      phys.col(1) = x;
      phys.col(2) = 1+x.array()*x.array();
      phys.col(3) = 1+x.array().pow(2) - 3*x.array().pow(3)+ 4*x.array().pow(4)-2*x.array().pow(5)+x.array().pow(6);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevFftwTransform::IntegratorType::INTG);

      // Dealias
      spec.bottomRows(xN-nN).setZero();

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   op(nN);

      // Compute first derivative
      spec.topRows(nN) = op.diff(0,1)*spec.topRows(nN);

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int i = 0; i < howmany; ++i)
      {
         Array cheb = Array::Zero(xN);
         if(i == 0)
         {
            cheb.setConstant(0);
         } else if(i == 1)
         {
            cheb.setConstant(1);
         } else if(i == 2)
         {
            cheb = 2*x;
         } else if(i == 3)
         {
            cheb = 2*x.array() - 9*x.array().pow(2)+ 16*x.array().pow(3)-10*x.array().pow(4)+6*x.array().pow(5);
         }

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), cheb(j), 1e-13) << "i = " << i << " j = " << j;
         }
      }
   }

   /**
    * @brief Test c constant in by computing second derivative (fft transform loop)
    */
   TEST_F(ChebyshevOperatorTest, Diff2)
   {
      // Set spectral and physical sizes
      int maxN = 6;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 5;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::EQUAL));

      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   phys = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   spec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      // Initialise the physical test data as cos(n*phi) + sin(n*phi) up to the highest maxN
      phys.setZero();
      phys.col(0).setConstant(1.0);
      phys.col(1) = x;
      phys.col(2) = 1+x.array()*x.array();
      phys.col(3) = 1+x.array().pow(2) - 3*x.array().pow(3)+ 4*x.array().pow(4)-2*x.array().pow(5)+x.array().pow(6);

      // Compute forward transform
      fft.integrate<Arithmetics::SET>(spec, phys, Transform::ChebyshevFftwTransform::IntegratorType::INTG);

      // Dealias
      spec.bottomRows(xN-nN).setZero();

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   op(nN);

      // Compute first derivative
      spec.topRows(nN) = op.diff(0,2)*spec.topRows(nN);

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int i = 0; i < howmany; ++i)
      {
         Array cheb = Array::Zero(xN);
         if(i == 0)
         {
            cheb.setConstant(0);
         } else if(i == 1)
         {
            cheb.setConstant(0);
         } else if(i == 2)
         {
            cheb.setConstant(2);
         } else if(i == 3)
         {
            cheb = 2 - 18*x.array()+ 48*x.array().pow(2)-40*x.array().pow(3)+30*x.array().pow(4);
         }

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), cheb(j), 1e-12) << "i = " << i << " j = " << j;
         }
      }
   }

   /**
    * @brief Test quasi identity through product
    */
   TEST_F(ChebyshevOperatorTest, QuasiIdentity)
   {
      // Set size of the basis
      int nN = 16;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      SparseMatrix sparse = spec.qDiff(1,0)*spec.diff(0,1);
      Matrix qid = sparse;

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 0 && i == j)
            {
               EXPECT_EQ(qid(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qid(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test quasi inverse of order 1
    */
   TEST_F(ChebyshevOperatorTest, QuasiInverse10)
   {
      // Set size of the basis
      int nN = 16;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      Matrix qdiff = spec.qDiff(1,0);

      //std::cerr << qdiff << std::endl;

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 0 && j < nN - 1)
            {
               //EXPECT_EQ(qdiff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qdiff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test quasi inverse of order 2
    */
   TEST_F(ChebyshevOperatorTest, QuasiInverse20)
   {
      // Set size of the basis
      int nN = 16;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      Matrix qdiff = spec.qDiff(2,0);

      //std::cerr << qdiff << std::endl;

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 1 && j < nN - 2)
            {
               //EXPECT_EQ(qdiff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qdiff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test quasi inverse of order 4
    */
   TEST_F(ChebyshevOperatorTest, QuasiInverse40)
   {
      // Set size of the basis
      int nN = 16;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      Matrix qdiff = spec.qDiff(4,0);

      //std::cerr << spec.qDiff(4,0) << std::endl;

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 3 && j < nN - 4)
            {
               //EXPECT_EQ(qdiff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qdiff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test quasi inverse of order 1 on derivative of order 1
    */
   TEST_F(ChebyshevOperatorTest, QuasiInverse11)
   {
      // Set size of the basis
      int nN = 16;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      Matrix qdiff = spec.qDiff(1,1);

      // Check operator
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 0 && i==j)
            {
               EXPECT_EQ(qdiff(i,j), 1.0) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qdiff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test quasi inverse of order 2 on derivative of order 1
    */
   TEST_F(ChebyshevOperatorTest, QuasiInverse21)
   {
      // Set size of the basis
      int nN = 16;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      Matrix qdiff = spec.qDiff(2,1);

      // Check operator
      Matrix corr = spec.qDiff(1,0);
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 1)
            {
               EXPECT_EQ(qdiff(i,j), corr(i,j)) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qdiff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
         }
      }
   }

   /**
    * @brief Test quasi inverse of order 4 on derivative of order 2
    */
   TEST_F(ChebyshevOperatorTest, QuasiInverse42)
   {
      // Set size of the basis
      int nN = 16;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   spec(nN);

      // Get quasi-inverse of order 1
      Matrix qdiff = spec.qDiff(4,2);

      // Check operator
      Matrix corr = spec.qDiff(2,0);
      for(int j = 0; j < nN; ++j)
      {
         for(int i = 0; i < nN; ++i)
         {
            if(i > 3)
            {
               EXPECT_EQ(qdiff(i,j), corr(i,j)) << "i = " << i << " j = " << j;
            } else
            {
               EXPECT_EQ(qdiff(i,j), 0.0) << "i = " << i << " j = " << j;
            }
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
