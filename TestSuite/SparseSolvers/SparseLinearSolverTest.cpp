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
#include "BoundaryCondition/BoundaryCondition.hpp"
#include "FastTransforms/ChebyshevFftwTransform.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/FftSelector.hpp"
#include "SpectralOperators/TauChebyshev.hpp"
#include "SpectralOperators/GalerkinChebyshev.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   typedef Eigen::Array<MHDFloat, Eigen::Dynamic, 1>   CoeffArray;

   namespace internal {
      void setSolution1DA(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId);
      void setSolution1DB(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId);
      void setSolution1DC(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId);
      void setSolution1DD(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId);
      void setSolution1DE(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId);
      void setSolution1DF(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId);
   }

   /**
    * @brief Test fixture for the sparse direct solvers
    */
   class SparseLinearSolverTest : public ::testing::Test {
      public:

         void solveProblem(Matrix& sol, const Matrix& rhs, const SparseMatrix& matA);
         void setupProblem(Matrix& sol, Matrix& rhs, const Array& param, const ArrayI& bcId, void (*pSolFct)(Matrix&, Matrix&, const CoeffArray&, const Array&, const ArrayI&));

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

   SparseLinearSolverTest::SparseLinearSolverTest()
      : mError(1e-12), mMaxN(63), mHowmany(1)
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

   void SparseLinearSolverTest::solveProblem(Matrix& sol, const Matrix& rhs, const SparseMatrix& matA)
   {
      // Create solver
      Solver::SparseSelector<SparseMatrix>::Type   solver;

      solver.compute(matA);

      sol.resize(rhs.rows(),1);
      sol.setConstant(-4242);
      sol = solver.solve(rhs);
   }

   void SparseLinearSolverTest::setupProblem(Matrix& sol, Matrix& rhs, const Array& param, const ArrayI& bcId, void (*pSolFct)(Matrix&, Matrix&, const CoeffArray&, const Array&, const ArrayI&))
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
      Solver::SparseSelector<SparseMatrix>::Type   solver;

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
      Solver::SparseSelector<SparseMatrixZ>::Type   solver;

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

   /**
    * @brief Cartesian 1D \f[\nabla^2 x = b\f] with QI approach
    *
    * @param SparseLinearSolverTest Test fixture ID
    * @param HarmonicDirichlet1D    Test ID
    */
   TEST_F(SparseLinearSolverTest, HarmonicQI1D)
   {
      // Loop over boundary condition ids
      for(int ibc = 21; ibc <= 23; ++ibc)
      {
         int nN = this->mMaxN + 1;
         int minN = 9;
         Array param(1);
         param(0) = 5;
         ArrayI bcId(1);
         bcId(0) = ibc;

         // Setup problem
         Matrix exactSol;
         Matrix exactRhs;
         this->setupProblem(exactSol, exactRhs, param, bcId, &internal::setSolution1DA);

         // Define boundary condition flags
         Boundary::BCVector bcs;
         if(ibc == 21)
         {
            bcs.push_back(Boundary::BoundaryCondition(Boundary::VALUE, Boundary::LEFT));
            bcs.push_back(Boundary::BoundaryCondition(Boundary::VALUE, Boundary::RIGHT));
         } else if(ibc == 22)
         {
            bcs.push_back(Boundary::BoundaryCondition(Boundary::D1, Boundary::LEFT));
            bcs.push_back(Boundary::BoundaryCondition(Boundary::D1, Boundary::RIGHT));
         } else if(ibc == 23)
         {
            bcs.push_back(Boundary::BoundaryCondition(Boundary::D2, Boundary::LEFT));
            bcs.push_back(Boundary::BoundaryCondition(Boundary::D2, Boundary::RIGHT));
         }

         // Create matrix
         Spectral::ChebyshevOperator   op(nN);
         SparseMatrix  matA(nN,nN);
         matA = Spectral::PeriodicOperator::qLaplacian2D(op, param(0), 2);
         Spectral::TauChebyshev tauBc(matA.rows(), bcs, 2);
         Spectral::GalerkinChebyshev galBc(matA.rows(), bcs, 2);

         SparseMatrix matT = tauBc.constrain(matA);
         SparseMatrix matG = galBc.constrain(matA);

         SparseMatrix matQ = op.qDiff(2,0);
         SparseMatrix matQT = tauBc.constrain(matQ);
         SparseMatrix matQG = galBc.constrain(matQ);

         Matrix rhs = matQ*exactRhs;
         Matrix rhsT = tauBc.restrict(rhs);
         Matrix rhsG = galBc.restrict(rhs);

         // Solve Tau problem
         Matrix solT;
         this->solveProblem(solT, rhsT, matT);
         solT = tauBc.extend(solT);

         // Check Tau solution
         for(int j = 0; j < exactSol.cols(); ++j)
         {
            for(int i = 0; i < exactSol.rows(); ++i)
            {
               MHDFloat eta = std::max(10.0, std::abs(exactSol(i,j)));
               EXPECT_NEAR(solT(i)/eta, exactSol(i)/eta, this->mError) << "Tau" << ibc;
            }
         }

         // Solve Galerkin problem
         Matrix solG;
         this->solveProblem(solG, rhsG, matG);
         solG = galBc.extend(solG);

         // Check Galerkin solution
         for(int j = 0; j < exactSol.cols(); ++j)
         {
            for(int i = 0; i < exactSol.rows(); ++i)
            {
               MHDFloat eta = std::max(10.0, std::abs(exactSol(i,j)));
               EXPECT_NEAR(solG(i)/eta, exactSol(i)/eta, this->mError) << "Galerkin: " << ibc;
            }
         }
//
//         // Test minimal resolution for Tau
//         op = Spectral::ChebyshevOperator(minN);
//         matA = Spectral::PeriodicOperator::qLaplacian2D(op, param(0), 2) + Spectral::BoundaryConditions::tauMatrix(bc, ids).first;
//         matQ = op.qDiff(2,0);
//         rhs = matQ*exactRhs.block(0,0,minN, exactRhs.cols());
//         this->solveProblem(sol, rhs, matA);
//
//         // Test solution
//         for(int j = 0; j < sol.cols(); ++j)
//         {
//            for(int i = 0; i < sol.rows(); ++i)
//            {
//               MHDFloat eta = std::max(10.0, std::abs(exactSol(i,j)));
//               EXPECT_NEAR(sol(i,j)/eta, exactSol(i,j)/eta, this->mError) << "Minial Tau" << ibc;
//            }
//         }
      }
   }

   /**
    * @brief  Cartesian 1D \f[\nabla^4 x = b\f] with QI approach
    *
    * @param SparseLinearSolverTest       Test fixture ID
    * @param BiHarmonicDirichlet1D  Test ID
    */
   TEST_F(SparseLinearSolverTest, BiHarmonicQI1D)
   {
//      // Loop over boundary condition ids
//      for(int ibc = 42; ibc <= 43; ++ibc)
//      {
//         int nN = this->mMaxN + 1;
//         int minN = 9;
//         Array param(1);
//         param(0) = 5;
//         ArrayI bcId(1);
//         bcId(0) = ibc;
//
//         // Setup problem
//         Matrix exactSol;
//         Matrix exactRhs;
//         Matrix rhs;
//         this->setupProblem(exactSol, exactRhs, param, bcId, &internal::setSolution1DB);
//
//         // Create matrix
//         Spectral::ChebyshevOperator   op(nN);
//         SparseMatrix  matA(nN,nN);
//         SparseMatrix  matQ(nN,nN);
//
//         std::vector<std::pair<Boundary::BCType,Boundary::BCPosition> > ids;
//         if(ibc == 42)
//         {
//            ids.push_back(std::make_pair(Boundary::VALUE, Boundary::LEFT));
//            ids.push_back(std::make_pair(Boundary::VALUE, Boundary::RIGHT));
//            ids.push_back(std::make_pair(Boundary::D1, Boundary::LEFT));
//            ids.push_back(std::make_pair(Boundary::D1, Boundary::RIGHT));
//         } else if(ibc == 43)
//         {
//            ids.push_back(std::make_pair(Boundary::VALUE, Boundary::LEFT));
//            ids.push_back(std::make_pair(Boundary::VALUE, Boundary::RIGHT));
//            ids.push_back(std::make_pair(Boundary::D2, Boundary::LEFT));
//            ids.push_back(std::make_pair(Boundary::D2, Boundary::RIGHT));
//         }
//         matA = Spectral::PeriodicOperator::qBilaplacian2D(op, param(0), 4) + Spectral::BoundaryConditions::tauMatrix(bc, ids).first;
//         matQ = op.qDiff(4,0);
//
//         // Solve
//         Matrix sol;
//         rhs = matQ*exactRhs;
//         this->solveProblem(sol, rhs, matA);
//
//         // Test solution
//         for(int j = 0; j < exactSol.cols(); ++j)
//         {
//            for(int i = 0; i < exactSol.rows(); ++i)
//            {
//               MHDFloat eta = std::max(10.0, std::abs(exactSol(i,j)));
//               EXPECT_NEAR(sol(i,j)/eta, exactSol(i,j)/eta, this->mError);
//            }
//         }
//
//         // Test minimal resolution
//         op = Spectral::ChebyshevOperator(minN);
//         matA = Spectral::PeriodicOperator::qBilaplacian2D(op, param(0), 4) + Spectral::BoundaryConditions::tauMatrix(bc, ids).first;
//         matQ = op.qDiff(4,0);
//         rhs = matQ*exactRhs.block(0,0,minN, exactRhs.cols());
//         this->solveProblem(sol, rhs, matA);
//
//         // Test solution
//         for(int j = 0; j < sol.cols(); ++j)
//         {
//            for(int i = 0; i < sol.rows(); ++i)
//            {
//               MHDFloat eta = std::max(10.0, std::abs(exactSol(i,j)));
//               EXPECT_NEAR(sol(i,j)/eta, exactSol(i,j)/eta, this->mError);
//            }
//         }
//      }
   }

   /**
    * @brief Cartesian 1D coupled system \f{eqnarray*}D_x u + v & = & b_1\\ -u + D_x v & = & b_2\f} with QI approach
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param HarmonicDirichlet2D       Test ID
    */
   TEST_F(SparseLinearSolverTest, CoupledDQI1D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 1D coupled system \f{eqnarray*}\nabla^2 u + v & = & b_1\\ -u + \nabla^2 v & = & b_2\f} with QI approach
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param HarmonicDirichlet2D       Test ID
    */
   TEST_F(SparseLinearSolverTest, CoupledLaplacianQI1D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 1D coupled system \f{eqnarray*}\nabla^4 u + v & = & b_1\\ -u + \nabla^2 v & = & b_2\f} with QI approach
    *
    * @param SparseLinearSolverTest    Test fixture ID
    * @param HarmonicDirichlet2D       Test ID
    */
   TEST_F(SparseLinearSolverTest, CoupledMixedLaplacianQI1D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f[\nabla_x^2 \otimes\nabla_y^2 x = b\f] with QI approach 
    *
    * @param SparseLinearSolverTest Test fixture ID
    * @param HarmonicDirichlet2D    Test ID
    */
   TEST_F(SparseLinearSolverTest, HarmonicQI2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f[\nabla_x^4 \otimes\nabla_y^4 x = b\f] with QI approach 
    *
    * @param SparseLinearSolverTest Test fixture ID
    * @param BiHarmonicDirichlet2D  Test ID
    */
   TEST_F(SparseLinearSolverTest, BiHarmonicQI2D)
   {
      ASSERT_TRUE(false);
   }

   /**
    * @brief Cartesian 2D \f[\nabla_x^4 \otimes\nabla_y^2 x = b\f] with QI approach 
    *
    * @param SparseLinearSolverTest Test fixture ID
    * @param MixedHarmonicMixed2D   Test ID
    */
   TEST_F(SparseLinearSolverTest, MixedHarmonicQI2D)
   {
      ASSERT_TRUE(false);
   }

   namespace internal 
   {
      void setSolution1DA(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId)
      {
         MHDFloat k = param(0);

         if(bcId(0) == 21)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = -(-1 + x.pow(2)).pow(2)*(8 - 56*x.pow(2) + k*k*(-1 + x.pow(2)).pow(2));
         } else if(bcId(0) == 22)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = -(-1 + x.pow(2)).pow(2)*(8 - 56*x.pow(2) + k*k*(-1 + x.pow(2)).pow(2));
         } else if(bcId(0) == 23)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = -(-1 + x.pow(2)).pow(2)*(8 - 56*x.pow(2) + k*k*(-1 + x.pow(2)).pow(2));
         }
      }

      void setSolution1DB(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId)
      {
         MHDFloat k = param(0);

         if(bcId(0) == 42)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = std::pow(k,4)*(-1 + x.pow(2)).pow(4)-16*k*k*(-1 + x.pow(2)).pow(2)*(-1 + 7.*x.pow(2)) + 48.*(3. - 30.*x.pow(2) + 35.*x.pow(4));
         } else if(bcId(0) == 43)
         {
            sol = ((1 + x)*(1 - x)).pow(4);
            rhs = std::pow(k,4)*(-1 + x.pow(2)).pow(4)-16.*k*k*(-1 + x.pow(2)).pow(2)*(-1 + 7.*x.pow(2)) + 48.*(3. - 30.*x.pow(2) + 35.*x.pow(4));
         }
      }

      void setSolution1DC(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId)
      {
         sol.resize(2*x.size(), 1);
         rhs.resize(2*x.size(), 1);

         if(bcId(0) == 1 && bcId(1) == -1)
         {
            sol.topRows(x.size()) = (1 - x).pow(8);
            sol.bottomRows(x.size()) = (1 + x).pow(8);
            rhs.topRows(x.size()) = 8.*(-1 + x).pow(7) + (1 + x).pow(8);
            rhs.bottomRows(x.size()) = -(-1 + x).pow(8) + 8.*(1 + x).pow(7);
         } else if(bcId(0) == 2 && bcId(1) == -2)
         {
            sol.topRows(x.size()) = (1 - x).pow(8);
            sol.bottomRows(x.size()) = (1 + x).pow(8);
            rhs.topRows(x.size()) = 8.*(-1 + x).pow(7) + (1 + x).pow(8);
            rhs.bottomRows(x.size()) = -(-1 + x).pow(8) + 8.*(1 + x).pow(7);
         } else if(bcId(0) == 3 && bcId(1) == -3)
         {
            sol.topRows(x.size()) = (1 - x).pow(8);
            sol.bottomRows(x.size()) = (1 + x).pow(8);
            rhs.topRows(x.size()) = 8.*(-1 + x).pow(7) + (1 + x).pow(8);
            rhs.bottomRows(x.size()) = -(-1 + x).pow(8) + 8.*(1 + x).pow(7);
         }
      }

      void setSolution1DD(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId)
      {
         sol.resize(2*x.size(), 1);
         rhs.resize(2*x.size(), 1);

         if(bcId(0) == 1 && bcId(1) == -1)
         {
            sol.topRows(x.size()) = (1-x).pow(8) + (1+x).pow(8);
            sol.bottomRows(x.size()) = (1-x).pow(8) - (1+x).pow(8);
            rhs.topRows(x.size()) = 32.*x*(3+7.*x.pow(2)*(2+x.pow(2)));
            rhs.bottomRows(x.size()) = -2.*(9+196.*x.pow(2)+350.*x.pow(4)+84.*x.pow(6)+x.pow(8));
         } else if(bcId(0) == 2 && bcId(1) == -2)
         {
            sol.topRows(x.size()) = (1-x).pow(8);
            sol.bottomRows(x.size()) = (1+x).pow(8);
            rhs.topRows(x.size()) = 8.*(-1+x).pow(7)+(1+x).pow(8);
            rhs.bottomRows(x.size()) = -(-1+x).pow(8)+8.*(1+x).pow(7);
         } else if(bcId(0) == 3 && bcId(1) == -3)
         {
            sol.topRows(x.size()) = (1-x).pow(8);
            sol.bottomRows(x.size()) = (1+x).pow(8);
            rhs.topRows(x.size()) = 8.*(-1+x).pow(7)+(1+x).pow(8);
            rhs.bottomRows(x.size()) = -(-1+x).pow(8)+8.*(1+x).pow(7);
         }
      }

      void setSolution1DE(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId)
      {
         MHDFloat k = param(0);
         sol.resize(2*x.size(), 1);
         rhs.resize(2*x.size(), 1);

         if(bcId(0) == 21 && bcId(1) == 21)
         {
            sol.topRows(x.size()) = (1-x.pow(2)).pow(4);
            sol.bottomRows(x.size()) = (1-x.pow(2)).pow(4);
            rhs.topRows(x.size()) = -(-1+x.pow(2)).pow(2)*(7+k*k*(-1+x.pow(2)).pow(2)-x.pow(2)*(54+x.pow(2)));
            rhs.bottomRows(x.size()) = -(-1+x.pow(2)).pow(2)*(9-58.*x.pow(2)+x.pow(4)+k*k*(-1+x.pow(2)).pow(2));
         } else if(bcId(0) == 21 && bcId(1) == 22)
         {
            sol.topRows(x.size()) = (1-x.pow(2)).pow(4);
            sol.bottomRows(x.size()) = (1-x.pow(2)).pow(4);
            rhs.topRows(x.size()) = -(-1+x.pow(2)).pow(2)*(7+k*k*(-1+x.pow(2)).pow(2)-x.pow(2)*(54+x.pow(2)));
            rhs.bottomRows(x.size()) = -(-1+x.pow(2)).pow(2)*(9-58.*x.pow(2)+x.pow(4)+k*k*(-1+x.pow(2)).pow(2));
         } else if(bcId(0) == 22 && bcId(1) == 22)
         {
            sol.topRows(x.size()) = (1-x.pow(2)).pow(4);
            sol.bottomRows(x.size()) = (1-x.pow(2)).pow(4);
            rhs.topRows(x.size()) = -(-1+x.pow(2)).pow(2)*(7+k*k*(-1+x.pow(2)).pow(2)-x.pow(2)*(54+x.pow(2)));
            rhs.bottomRows(x.size()) = -(-1+x.pow(2)).pow(2)*(9-58.*x.pow(2)+x.pow(4)+k*k*(-1+x.pow(2)).pow(2));
         }
      }

      void setSolution1DF(Matrix& sol, Matrix& rhs, const CoeffArray& x, const Array& param, const ArrayI& bcId)
      {
         MHDFloat k = param(0);

         if(bcId(0) == 21 && bcId(1) == 42)
         {
            sol.topRows(x.size()) = (1-x.pow(2)).pow(4);
            sol.bottomRows(x.size()) = (1-x.pow(2)).pow(4);
            rhs.topRows(x.size()) = -(-1+x.pow(2)).pow(2)*(7+k*k*(-1+x.pow(2)).pow(2)-x.pow(2)*(54+x.pow(2)));
            rhs.bottomRows(x.size()) = 143 - 1436.*x.pow(2) + 1674.*x.pow(4) + 4.*x.pow(6) - x.pow(8) + std::pow(k,4)*(-1+x.pow(2)).pow(4) - 16.*k*k*(-1+x.pow(2)).pow(2)*(-1+7.*x.pow(2));
         } else if(bcId(0) == 21 && bcId(1) == 43)
         {
            sol.topRows(x.size()) = (1-x.pow(2)).pow(4);
            sol.bottomRows(x.size()) = (1-x.pow(2)).pow(4);
            rhs.topRows(x.size()) = -(-1+x.pow(2)).pow(2)*(7+k*k*(-1+x.pow(2)).pow(2)-x.pow(2)*(54+x.pow(2)));
            rhs.bottomRows(x.size()) = 143 - 1436.*x.pow(2) + 1674.*x.pow(4) + 4.*x.pow(6) - x.pow(8) + std::pow(k,4)*(-1+x.pow(2)).pow(4) - 16.*k*k*(-1+x.pow(2)).pow(2)*(-1+7.*x.pow(2));
         } else if(bcId(0) == 22 && bcId(1) == 42)
         {
            sol.topRows(x.size()) = (1-x.pow(2)).pow(4);
            sol.bottomRows(x.size()) = (1-x.pow(2)).pow(4);
            rhs.topRows(x.size()) = -(-1+x.pow(2)).pow(2)*(7+k*k*(-1+x.pow(2)).pow(2)-x.pow(2)*(54+x.pow(2)));
            rhs.bottomRows(x.size()) = 143 - 1436.*x.pow(2) + 1674.*x.pow(4) + 4.*x.pow(6) - x.pow(8) + std::pow(k,4)*(-1+x.pow(2)).pow(4) - 16.*k*k*(-1+x.pow(2)).pow(2)*(-1+7.*x.pow(2));
         } else if(bcId(0) == 22 && bcId(1) == 43)
         {
            sol.topRows(x.size()) = (1-x.pow(2)).pow(4);
            sol.bottomRows(x.size()) = (1-x.pow(2)).pow(4);
            rhs.topRows(x.size()) = -(-1+x.pow(2)).pow(2)*(7+k*k*(-1+x.pow(2)).pow(2)-x.pow(2)*(54+x.pow(2)));
            rhs.bottomRows(x.size()) = 143 - 1436.*x.pow(2) + 1674.*x.pow(4) + 4.*x.pow(6) - x.pow(8) + std::pow(k,4)*(-1+x.pow(2)).pow(4) - 16.*k*k*(-1+x.pow(2)).pow(2)*(-1+7.*x.pow(2));
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
