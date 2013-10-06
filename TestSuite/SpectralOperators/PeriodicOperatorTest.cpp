/** 
 * @file PeriodicOperatorTest.cpp
 * @brief Implementation of test case for PeriodicOperator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "SpectralOperators/ChebyshevOperator.hpp"
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"
#include "FastTransforms/ChebyshevFftwTransform.hpp"
#include "FastTransforms/FftwTools.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the PeriodicOperator implementation
    */
   class PeriodicOperatorTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         PeriodicOperatorTest();

         /**
          * @brief Destructor
          */
         virtual ~PeriodicOperatorTest();

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
   };

   PeriodicOperatorTest::PeriodicOperatorTest()
   {
   }

   PeriodicOperatorTest::~PeriodicOperatorTest()
   {
   }

//   void PeriodicOperatorTest::SetUp()
//   {
//   }

//   void PeriodicOperatorTest::TearDown()
//   {
//   }

   TEST_F(PeriodicOperatorTest, Laplacian2D)
   TEST_F(PeriodicOperatorTest, Laplacian2D_F)
   TEST_F(PeriodicOperatorTest, Laplacian2D_FF)
   TEST_F(PeriodicOperatorTest, Bilaplacian2D)
   TEST_F(PeriodicOperatorTest, Bilaplacian2D_F)
   TEST_F(PeriodicOperatorTest, Bilaplacian2D_FF)
   TEST_F(PeriodicOperatorTest, Laplacian3D)
   TEST_F(PeriodicOperatorTest, Laplacian3D_F)
   TEST_F(PeriodicOperatorTest, Laplacian3D_FF)
   TEST_F(PeriodicOperatorTest, Laplacian3D_FFF)
   TEST_F(PeriodicOperatorTest, Bilaplacian3D)
   TEST_F(PeriodicOperatorTest, Bilaplacian3D_F)
   TEST_F(PeriodicOperatorTest, Bilaplacian3D_FF)
   TEST_F(PeriodicOperatorTest, Bilaplacian3D_FFF)
   TEST_F(PeriodicOperatorTest, QLaplacian2D)
   TEST_F(PeriodicOperatorTest, QLaplacian2D_F)
   TEST_F(PeriodicOperatorTest, QLaplacian2D_FF)
   TEST_F(PeriodicOperatorTest, QBilaplacian2D)
   TEST_F(PeriodicOperatorTest, QBilaplacian2D_F)
   TEST_F(PeriodicOperatorTest, QBilaplacian2D_FF)
   TEST_F(PeriodicOperatorTest, QLaplacian3D)
   TEST_F(PeriodicOperatorTest, QLaplacian3D_F)
   TEST_F(PeriodicOperatorTest, QLaplacian3D_FF)
   TEST_F(PeriodicOperatorTest, QLaplacian3D_FFF)
   TEST_F(PeriodicOperatorTest, QBilaplacian3D)
   TEST_F(PeriodicOperatorTest, QBilaplacian3D_F)
   TEST_F(PeriodicOperatorTest, QBilaplacian3D_FF)
   TEST_F(PeriodicOperatorTest, QBilaplacian3D_FFF)

   /**
    * @brief Test default constructor
    */
   TEST_F(PeriodicOperatorTest, BasicM0)
   {
      // Set spectral and physical sizes
      int maxN = 12;
      int nN = maxN + 1;

      // Create Chebyshev spectral operator
      Spectral::ChebyshevOperator   op(nN);
   }

   /**
    * @brief Test default constructor
    */
   TEST_F(PeriodicOperatorTest, Laplacian2D)
   {
      // Set spectral and physical sizes
      int maxN = 12;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 5;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::REAL));

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
      int m = 4;
      spec.topRows(nN) = Spectral::BoxTools::laplacian2D(op,m,0)*spec.topRows(nN);

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int i = 0; i < howmany; ++i)
      {
         Array cheb = Array::Zero(xN);
         if(i == 0)
         {
            cheb.setConstant(-m*m);
         } else if(i == 1)
         {
            cheb = -m*m*x;
         } else if(i == 2)
         {
            cheb = 2 - m*m*(1+x.array().pow(2));
         } else if(i == 3)
         {
            cheb = 2 - 18*x.array()+ 48*x.array().pow(2)-40*x.array().pow(3)+30*x.array().pow(4) - m*m*(1+x.array().pow(2) - 3*x.array().pow(3)+ 4*x.array().pow(4)-2*x.array().pow(5)+x.array().pow(6));
         }

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), cheb(j), 2e-10) << "i = " << i << " j = " << j;
         }
      }

   }

   /**
    * @brief Test default constructor
    */
   TEST_F(PeriodicOperatorTest, Bilaplacian2D)
   {
      // Set spectral and physical sizes
      int maxN = 12;
      int nN = maxN + 1;
      int xN = Transform::FftwTools::dealiasFft(nN);
      int howmany = 5;

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, howmany, nN, Transform::FftSetup::REAL));

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
      int m = 4;
      spec.topRows(nN) = Spectral::BoxTools::bilaplacian2D(op,m,0)*spec.topRows(nN);

      // Compute backward transform
      fft.project<Arithmetics::SET>(phys, spec, Transform::ChebyshevFftwTransform::ProjectorType::PROJ);

      // Check the solution
      for(int i = 0; i < howmany; ++i)
      {
         Array cheb = Array::Zero(xN);
         if(i == 0)
         {
            cheb.setConstant(static_cast<MHDFloat>(std::pow(m,4)));
         } else if(i == 1)
         {
            cheb = std::pow(m,4)*x;
         } else if(i == 2)
         {
            cheb = -2*std::pow(m,2)-std::pow(m,2)*(2-std::pow(m,2)*(1+x.array().pow(2)));
         } else if(i == 3)
         {
            cheb = 96 - 4*std::pow(m,2) + std::pow(m,4) + (-240 + 36*std::pow(m,2))*x.array() + (360 - 96*std::pow(m,2) + std::pow(m,4))*x.array().pow(2)+(80*std::pow(m,2)-3*std::pow(m,4))*x.array().pow(3) + (-60*std::pow(m,2)+4*std::pow(m,4))*x.array().pow(4)-2*std::pow(m,4)*x.array().pow(5) + std::pow(m,4)*x.array().pow(6);
         }

         for(int j = 0; j < x.size(); ++j)
         {
            EXPECT_NEAR(phys(j,i), cheb(j), 1e-8) << "i = " << i << " j = " << j;
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
