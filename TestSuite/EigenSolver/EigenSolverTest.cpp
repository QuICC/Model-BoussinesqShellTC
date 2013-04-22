/** \file EigenSolverTest.cpp
 *  \brief Implementation of test case for EigenSolver
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>
#include "Base/MathConstants.hpp"
#include "Timesteppers/SparseSolverMacro.h"
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGSystem.hpp"
#include "../External/Interfaces/ARPACK_Interface.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the EigenSolver implementation
    */
   class EigenSolverTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         EigenSolverTest();

         /**
          * @brief Destructor
          */
         virtual ~EigenSolverTest();

         void arpackSolver(const SparseMatrixZ& lhs, const SparseMatrixZ& rhs);

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   EigenSolverTest::EigenSolverTest()
   {
   }

   EigenSolverTest::~EigenSolverTest()
   {
   }

//   void EigenSolverTest::SetUp()
//   {
//   }

//   void EigenSolverTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(EigenSolverTest, Play)
   {
      // Create shared simulation boundary
      SharedSimulationBoundary  spBcs(new SimulationBoundary());

      // Storage for the dimension ID
      Dimensions::Simulation::Id dimId;

      // Create equation and field keys
      SpectralFieldId eqId;
      SpectralFieldId fieldId;

      // Temperature equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT);
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT);

      // Streamfunction equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      // No-slip boundary conditions
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::LEFT); 
         spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::RIGHT); 
      // Stress-free boundary conditions
      // spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::LEFT); 
      // spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::SECOND_DERIVATIVE, Spectral::IBoundary::RIGHT);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::LEFT); 
      //    ... coupled boundary conditions
      fieldId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 

      // Axial velocity equation
      //    ... boundary conditions
      eqId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      spBcs->initStorage(eqId);
      dimId = Dimensions::Simulation::SIM1D;
      fieldId = std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR);
      spBcs->initBcStorage(eqId, fieldId, dimId);
      // No-slip boundary conditions
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::LEFT); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      // Stress-free boundary conditions
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::LEFT); 
      //spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::FIRST_DERIVATIVE, Spectral::IBoundary::RIGHT); 
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId);
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::VALUE, Spectral::IBoundary::RIGHT); 
      //    ... coupled boundary conditions
      fieldId = std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR);
      dimId = Dimensions::Simulation::SIM3D;
      spBcs->initBcStorage(eqId, fieldId, dimId); 
      spBcs->addBc(eqId, fieldId, dimId, Spectral::BoundaryConditions::BETA_SLOPE, Spectral::IBoundary::RIGHT); 

      // Get X and Z dimensions
      int nx = 7;
      int nz = 7;

      // Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = 800.0;
      MHDFloat Pr = 1.0;
      MHDFloat Gamma = .1;
      MHDFloat chi = 45;

      // Wave number
      MHDFloat k_ = 2.75;

      // Create vector of fields
      std::vector<SpectralFieldId>  fieldIds;
      fieldIds.push_back(std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR));
      fieldIds.push_back(std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR));
      fieldIds.push_back(std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR));

      // Storage for the matrix row
      DecoupledZSparse  block = std::make_pair(SparseMatrix(nx*nz, nx*nz),SparseMatrix(nx*nz, nx*nz));

      SparseMatrix  tmp;
      SparseMatrixZ  sparseLhs(nx*nz*fieldIds.size(), nx*nz*fieldIds.size());
      SparseMatrixZ  sparseRhs(nx*nz*fieldIds.size(), nx*nz*fieldIds.size());

      // Loop over all coupled fields
      int rowIdx = 0;
      std::vector<SpectralFieldId>::const_iterator rowIt;
      std::vector<SpectralFieldId>::const_iterator colIt;
      for(rowIt = fieldIds.begin(); rowIt != fieldIds.end(); ++rowIt)
      {
         int colIdx = 0;
         for(colIt = fieldIds.begin(); colIt != fieldIds.end(); ++colIt)
         {
            SparseMatrix   blockMatrix(fieldIds.size(),fieldIds.size());
            blockMatrix.insert(rowIdx, colIdx) = 1;

            Equations::Beta3DQGSystem::linearBlock(block, *rowIt, *colIt, nx, nz, k_, Ra, Pr, Gamma, chi);
            Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
            sparseLhs += tmp.cast<MHDComplex>();
            Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
            sparseLhs += MathConstants::cI*tmp;

            std::cerr << sparseLhs << std::endl;

            Equations::Beta3DQGSystem::boundaryBlock(block, *rowIt, *colIt, spBcs, nx, nz, k_, Ra, Pr, Gamma, chi);
            Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
            sparseLhs += tmp.cast<MHDComplex>();
            Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
            sparseLhs += MathConstants::cI*tmp;

            if(*rowIt == *colIt)
            {
               Equations::Beta3DQGSystem::timeBlock(block, *rowIt, nx, nz, k_, Ra, Pr, Gamma, chi);
               Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
               sparseRhs += tmp.cast<MHDComplex>();
               Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
               sparseRhs += MathConstants::cI*tmp;
            }

            ++colIdx;
         }
         ++rowIdx;
      }

      // Make sure matrices are in compressed format
      sparseLhs.makeCompressed();
      sparseRhs.makeCompressed();

      std::cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
      std::cerr << "SPARSE LHS OPERATOR OF GENERALIZED EIGENPROBLEM: " << sparseLhs.rows() << " x " << sparseLhs.cols() << std::endl;
      std::cerr << sparseLhs << std::endl;
      std::cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
      std::cerr << "SPARSE RHS OPERATOR OF GENERALIZED EIGENPROBLEM: " << sparseRhs.rows() << " x " << sparseRhs.cols() << std::endl;
      std::cerr << sparseRhs << std::endl;

      this->arpackSolver(sparseLhs, sparseRhs);


   }

   void EigenSolverTest::arpackSolver(const SparseMatrixZ& lhs, const SparseMatrixZ& rhs)
   {
      // Create solver operator and RHS
      MHDComplex  sigma(0.0,0.0);
      SparseMatrixZ solveOp = lhs - sigma*rhs;
      SparseSolverMacro<SparseMatrixZ> solver;
      solver.compute(solveOp);

      int n = rhs.rows();
      int nev = 4;
      int ncv = 20;
      std::string bmat = "G";
      std::string which = "LM";

      int lworkl = std::pow(3*ncv,2) + 5*ncv;
      double tol = 0.0;
      int ido = 0; 
      int info = 0;

      ArrayI   iparam(11);
      iparam.setZero(); // ishfts
      iparam(0) = 1; // ishfts
      iparam(2) = 300; // maxitr
      iparam(6) = 3; // mode


      int ldv = n;
      ArrayZ resid(n);
      resid.setZero();
      MatrixZ v(n,ncv);
      v.setZero();
      ArrayI ipntr(14);
      ipntr.setZero();
      ArrayZ workd(3*n);
      workd.setZero();
      ArrayZ workl(lworkl);
      workl.setZero();
      Array rwork(ncv);
      rwork.setZero();

      while(ido != 99)
      {
         znaupd_(&ido, bmat.c_str(), &n, which.c_str(), &nev, &tol, resid.data(), &ncv, v.data(), &ldv, iparam.data(), ipntr.data(), workd.data(), workl.data(), &lworkl, rwork.data(), &info);

         if(ido == -1)
         {
            std::cerr << "IDO == -1" << std::endl;
            workd.segment(ipntr(1),n) = rhs*workd.segment(ipntr(0), n);
            workd.segment(ipntr(1),n) = solver.solve(workd.segment(ipntr(1),n));
         } else if(ido == 1)
         {
            std::cerr << "IDO == 1" << std::endl;
            workd.segment(ipntr(1),n) = workd.segment(ipntr(2), n);
            workd.segment(ipntr(1),n) = solver.solve(workd.segment(ipntr(1),n));
         } else if(ido == 2)
         {
            std::cerr << "IDO == 2" << std::endl;
            workd.segment(ipntr(1),n) = rhs*workd.segment(ipntr(0), n);
         } else if(ido == 3)
         {
            std::cerr << "IDO == 3" << std::endl;
         } else if(ido == 4)
         {
            std::cerr << "IDO == 4" << std::endl;
         } else if(ido == 99)
         {
            std::cerr << "IDO == 99" << std::endl;
         } else
         {
            std::cerr << "THIS SHOULD NOT HAVE HAPPENED!" << std::endl;
         }
      }

      if(info < 0)
      {
         std::cerr << "ZNAUPD failed!" << std::endl;
      } else
      {
         std::cerr << "ZNAUPD successful!" << std::endl;

         int rvec = 1;
         std::string howmny = "A";
         ArrayI select(ncv);
         select.setZero();
         ArrayZ d(nev);
         d.setZero();
         MatrixZ z(n,nev);
         z.setZero();
         int ldz = n;
         ArrayZ workev(2*ncv);
         workev.setZero();

         zneupd_(&rvec, howmny.c_str(), select.data(), d.data(), z.data(), &ldz, &sigma, workev.data(), bmat.c_str(), &n, which.c_str(), &nev, &tol, resid.data(), &ncv, v.data(), &ldv, iparam.data(), ipntr.data(), workd.data(), workl.data(), &lworkl, rwork.data(), &info);

         if(info < 0)
         {
            std::cerr << "ZNEUPD failed!" << std::endl;
            std::cerr << info << std::endl;
         } else
         {
            std::cerr << d.transpose() << std::endl;
         }
      }
   }

}
}

/// Main to execute all test from test case
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
