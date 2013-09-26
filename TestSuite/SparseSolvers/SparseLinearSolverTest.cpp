/** 
 * @file SparseLinearSolverTest.cpp
 * @brief Implementation of test cases for a generic sparse linear solver
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "Base/Typedefs.hpp"
#include "TypeSelectors/SparseSolverSelector.hpp"
#include "SpectralOperators/ChebyshevOperator.hpp"
#include "FastTransforms/ChebyshevFftwTransform.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/FftSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   typedef Eigen::Array<MHDFloat, Eigen::Dynamic, 1>   CoeffArray;

   namespace internal {
      void setSolution1(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const int bcId);
   }

   /**
    * @brief Test fixture for the sparse direct solvers
    */
   class SparseLinearSolverTest : public ::testing::Test {
      public:

         void solveProblem(Matrix& sol, const Matrix& rhs, const SparseMatrix& matA);
         void setupProblem(Matrix& sol, Matrix& rhs, const Array& param, const int bcId, void (*pSolFct)(Matrix&, Matrix&, const CoeffArray&, const Array&, const int));

      protected:
         /**
          * @brief Constructor
          */
         SparseLinearSolverTest();

         /**
          * @brief Destructor
          */
         virtual ~SparseLinearSolverTest();

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

   SparseLinearSolverTest::SparseLinearSolverTest()
      : mMaxN(31), mHowmany(1), mError(1e-10), mRelError(1e-10)
   {
   }

   SparseLinearSolverTest::~SparseLinearSolverTest()
   {
   }

//   void SparseLinearSolverTest::SetUp()
//   {
//   }

//   void SparseLinearSolverTest::TearDown()
//   {
//   }

   /**
    * @brief Small and simple real problem taken from Pardiso example 
    *
    * @param SparseLinearSolverTest Test fixture ID
    * @param SimpleReal       Test ID
    */
   TEST_F(SparseLinearSolverTest, SimpleReal)
   {
      // Build test matrix triplets
      typedef Eigen::Triplet<MHDFloat> T;
      std::vector<T> tripletList;
      tripletList.reserve(20);
      tripletList.push_back(T(0,0,7.0));
      tripletList.push_back(T(0,2,1.0));
      tripletList.push_back(T(0,5,2.0));
      tripletList.push_back(T(0,6,7.0));
      tripletList.push_back(T(1,1,-4.0));
      tripletList.push_back(T(1,2,8.0));
      tripletList.push_back(T(1,4,2.0));
      tripletList.push_back(T(2,2,1.0));
      tripletList.push_back(T(2,7,5.0));
      tripletList.push_back(T(3,3,7.0));
      tripletList.push_back(T(3,6,9.0));
      tripletList.push_back(T(4,1,-4.0));
      tripletList.push_back(T(5,2,7.0));
      tripletList.push_back(T(5,5,3.0));
      tripletList.push_back(T(5,7,8.0));
      tripletList.push_back(T(6,1,1.0));
      tripletList.push_back(T(6,6,11.0));
      tripletList.push_back(T(7,2,-3.0));
      tripletList.push_back(T(7,6,2.0));
      tripletList.push_back(T(7,7,5.0));

      // Build test matrix
      SparseMatrix   matA(8,8);
      matA.setFromTriplets(tripletList.begin(), tripletList.end());

      // Build RHS
      Array rhs(8);
      rhs.setLinSpaced(8,0.,7.);

      // Create solver
      SparseSolverSelector<SparseMatrix>::SolverType   solver;

      // Compute factorisation
      solver.compute(matA);

      // Compute solution
      Array sol(8);
      sol = solver.solve(rhs);

      // Exact solution
      Array exact(8);
      exact(0) = -1.153896103896104;exact(1) = -1;exact(2) = -0.9318181818181819;exact(3) = -0.3896103896103896;exact(4) = 2.227272727272728;exact(5) = 2.277272727272727;exact(6) = 0.6363636363636364;exact(7) = 0.5863636363636366;

      // Check solution
      for(int i = 0; i < rhs.size(); ++i)
      {
         EXPECT_NEAR(sol(i), exact(i), 1e-14) << "i = " << i;
      }
   }

   /**
    * @brief Small and simple complex problem taken from Pardiso example 
    *
    * @param SparseLinearSolverTest Test fixture ID
    * @param SimpleComplex    Test ID
    */
   TEST_F(SparseLinearSolverTest, SimpleComplex)
   {
      // Build test matrix triplets
      typedef Eigen::Triplet<MHDComplex> T;
      std::vector<T> tripletList;
      tripletList.reserve(20);
      tripletList.push_back(T(0,0,MHDComplex( 7.0 , 1.0)));
      tripletList.push_back(T(0,2,MHDComplex( 1.0 , 1.0)));
      tripletList.push_back(T(0,5,MHDComplex( 2.0 , 1.0)));
      tripletList.push_back(T(0,6,MHDComplex( 7.0 , 1.0)));
      tripletList.push_back(T(1,1,MHDComplex(-4.0 , 0.0)));
      tripletList.push_back(T(1,2,MHDComplex( 8.0 , 1.0)));
      tripletList.push_back(T(1,4,MHDComplex( 2.0 , 1.0)));
      tripletList.push_back(T(2,2,MHDComplex( 1.0 , 1.0)));
      tripletList.push_back(T(2,7,MHDComplex( 5.0 , 1.0)));
      tripletList.push_back(T(3,3,MHDComplex( 7.0 , 0.0)));
      tripletList.push_back(T(3,6,MHDComplex( 9.0 , 1.0)));
      tripletList.push_back(T(4,1,MHDComplex(-4.0 , 1.0)));
      tripletList.push_back(T(5,2,MHDComplex( 7.0 , 1.0)));
      tripletList.push_back(T(5,5,MHDComplex( 3.0 , 1.0)));
      tripletList.push_back(T(5,7,MHDComplex( 8.0 , 0.0)));
      tripletList.push_back(T(6,1,MHDComplex( 1.0 , 1.0)));
      tripletList.push_back(T(6,6,MHDComplex( 11.0, 1.0)));
      tripletList.push_back(T(7,2,MHDComplex(-3.0 , 1.0)));
      tripletList.push_back(T(7,6,MHDComplex( 2.0 , 1.0)));
      tripletList.push_back(T(7,7,MHDComplex( 5.0 , 0.0)));

      // Build test matrix
      SparseMatrixZ   matA(8,8);
      matA.setFromTriplets(tripletList.begin(), tripletList.end());

      // Build RHS
      ArrayZ rhs(8);
      rhs.setConstant(MHDComplex(1.0, 1.0));

      // Create solver
      SparseSolverSelector<SparseMatrixZ>::SolverType   solver;

      // Compute factorisation
      solver.compute(matA);

      // Compute solution
      ArrayZ sol(8);
      sol = solver.solve(rhs);

      // Exact solution
      ArrayZ exact(8);
      exact(0) = MHDComplex(0.174767984570878, 0.021177242044359);exact(1) = MHDComplex(-0.176470588235294, -0.294117647058824);exact(2) = MHDComplex(0.049321761491482, 0.029598199935712);exact(3) = MHDComplex(0.042981126876980, -0.031409285025486);exact(4) = MHDComplex(-0.120858887817422, -0.170859530697525);exact(5) = MHDComplex(-0.369347476695596, -0.000861459337833);exact(6) = MHDComplex(0.091610414657666, 0.125361620057859);exact(7) = MHDComplex(0.223940855030537, 0.139427836708454);

      // Check solution
      for(int i = 0; i < rhs.size(); ++i)
      {
         EXPECT_NEAR(sol(i).real(), exact(i).real(), 1e-14) << "Re i = " << i;
         EXPECT_NEAR(sol(i).imag(), exact(i).imag(), 1e-14) << "Im i = " << i;
      }
   }

   void SparseLinearSolverTest::solveProblem(Matrix& sol, const Matrix& rhs, const SparseMatrix& matA)
   {
      // Create solver
      SparseSolverSelector<SparseMatrix>::SolverType   solver;

      solver.compute(matA);

      sol.resize(rhs.rows(),1);
      sol.setConstant(-4242);
      sol = solver.solve(rhs);
   }

   void SparseLinearSolverTest::setupProblem(Matrix& sol, Matrix& rhs, const int bcId, void (*pSolFct)(Matrix&, Matrix&, const CoeffArray&, const Array&, const int))
   {
      // Set spectral and physical sizes
      int nN = this->mMaxN + 1;
      int xN = Transform::FftToolsType::dealiasCosFft(nN);

      // Create setup
      Transform::SharedFftSetup spSetup(new Transform::FftSetup(xN, this->mHowmany, nN, Transform::FftSetup::REAL));

      // Create ChebyshevFftwTransform
      Transform::ChebyshevFftwTransform fft;

      // Initialise fft
      fft.init(spSetup);

      // Create test data storage
      Matrix   physSol = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   physRhs = Matrix::Zero(spSetup->fwdSize(), spSetup->howmany());
      Matrix   tmpSpec = Matrix::Zero(spSetup->bwdSize(), spSetup->howmany());

      sol.resize(spSetup->specSize(), spSetup->howmany());
      rhs.resize(spSetup->specSize(), spSetup->howmany());

      // Get chebyshev grid
      Array x = fft.meshGrid();

      (*pSolFct)(physSol, physRhs, x, param, bcId);

      fft.integrate<Arithmetics::SET>(tmpSpec, physSol, Transform::ChebyshevFftwTransform::IntegratorType::INTG);
      sol = tmpSpec.topRows(spSetup->specSize());
      fft.integrate<Arithmetics::SET>(tmpSpec, physRhs, Transform::ChebyshevFftwTransform::IntegratorType::INTG);
      rhs = tmpSpec.topRows(spSetup->specSize());
   }

   namespace internal 
   {
      void setSolution1(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const int bcId)
      {
         MHDFloat k = param(0);

         if(bcId == 21)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = -(-1 + x.pow(2)).pow(2)*(8 - 56*x.pow(2) + k*k*(-1 + x.pow(2)).pow(2));
         } else if(bcId == 22)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = -(-1 + x.pow(2)).pow(2)*(8 - 56*x.pow(2) + k*k*(-1 + x.pow(2)).pow(2));
         } else if(bcId == 23)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = -(-1 + x.pow(2)).pow(2)*(8 - 56*x.pow(2) + k*k*(-1 + x.pow(2)).pow(2));
         }
      }
   }

   /**
    * @brief Cartesian 1D \f$\nabla^2 x = b\f$ with Dirichlet boundary conditions 
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param HarmonicDirichlet1D Test ID
    */
   TEST_F(SparseLinearSolverTest, HarmonicDirichlet1D)
   {
      int nN = this->mMaxN + 1;
      Array param(1);
      param(0) = 5;

      // Setup problem
      Matrix exactSol;
      Matrix rhs;
      this->setupProblem(exactSol, rhs, param, 21, &internal::setSolution1);

      // Create matrix
      Spectral::ChebyshevOperator   op(nN);
      SparseMatrix  matA(nN,nN);

      // Solve
      Matrix sol;
      this->solveProblem(sol, rhs, matA);

      // Test solution
      for(int i = 0; i < nN; ++i)
      {
         MHDFloat eta = std::max(1.0, std::abs(exactSol(i)));
         EXPECT_NEAR(sol(i)/eta, exactSol(i)/eta, this->mError);
      }

   }

   /**
    * @brief Cartesian 1D \f$\nabla^2 x = b\f$ with Neumann boundary conditions 
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param HarmonicNeumann1D   Test ID
    */
   TEST_F(SparseLinearSolverTest, HarmonicNeumann1D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief  Cartesian 1D \f$\nabla^4 x = b\f$ with Dirichlet boundary conditions
    *
    * @param SparseLinearSolverTest       Test fixture ID
    * @param BiHarmonicDirichlet1D  Test ID
    */
   TEST_F(SparseLinearSolverTest, BiHarmonicDirichlet1D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 1D \f$\nabla^4 x = b\f$ with Neumann boundary conditions 
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param BiHarmonicNeumann1D Test ID
    */
   TEST_F(SparseLinearSolverTest, BiHarmonicNeumann1D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^2 \otimes\nabla_y^2 x = b\f$ with Dirichlet/Dirichlet boundary conditions 
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param HarmonicDirichlet2D Test ID
    */
   TEST_F(SparseLinearSolverTest, HarmonicDirichlet2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^2 \otimes\nabla_y^2 x = b\f$ with Neumann/Neumann boundary conditions 
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param HarmonicNeumann2D   Test ID
    */
   TEST_F(SparseLinearSolverTest, HarmonicNeumann2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^2 \otimes\nabla_y^2 x = b\f$ with mixed Dirichlet/Neumann boundary conditions
    *
    * @param SparseLinearSolverTest Test fixture ID
    * @param HarmonicMixed2D  Test ID
    */
   TEST_F(SparseLinearSolverTest, HarmonicMixed2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^4 x = b\f$ with Dirichlet/Dirichlet boundary conditions 
    *
    * @param SparseLinearSolverTest       Test fixture ID
    * @param BiHarmonicDirichlet2D  Test ID
    */
   TEST_F(SparseLinearSolverTest, BiHarmonicDirichlet2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^4 x = b\f$ with Neumann/Neumann boundary conditions 
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param BiHarmonicNeumann2D Test ID
    */
   TEST_F(SparseLinearSolverTest, BiHarmonicNeumann2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^4 x = b\f$ with mixed Dirichlet/Neumann boundary conditions 
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param BiHarmonicMixed2D   Test ID
    */
   TEST_F(SparseLinearSolverTest, BiHarmonicMixed2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f$\nabla_x^4 \otimes\nabla_y^2 x = b\f$ with mixed Dirichlet/Neumann boundary conditions 
    *
    * @param SparseLinearSolverTest       Test fixture ID
    * @param MixedHarmonicMixed2D   Test ID
    */
   TEST_F(SparseLinearSolverTest, MixedHarmonicMixed2D)
   {
      ASSERT_TRUE(false);
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
