/** 
 * @file EigenSolverTest.cpp
 * @brief Implementation of test cases for the generic eigen solver
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"

#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>
#include "Base/MathConstants.hpp"
#include "Exceptions/Exception.hpp"
#include "Timesteppers/SparseSolverMacro.h"
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGSystem.hpp"
#include "../External/Interfaces/ARPACK_Interface.h"
#include "EigenSolver/ArpackEigenSolver.hpp"
#include <cmath>

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
    * @brief Not a real test yet, just a toy
    *
    * @param EigenSolverTest  Test fixture ID
    * @param Play             Test ID
    */
   TEST_F(EigenSolverTest, Play)
   {
      // Create an ARPACK eigen solver
      ArpackEigenSolver eigenSolver;

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

      // Create vector of fields
      std::vector<SpectralFieldId>  fieldIds;
      fieldIds.push_back(std::make_pair(PhysicalNames::STREAMFUNCTION, FieldComponents::Spectral::SCALAR));
      fieldIds.push_back(std::make_pair(PhysicalNames::VELOCITYZ, FieldComponents::Spectral::SCALAR));
      fieldIds.push_back(std::make_pair(PhysicalNames::TEMPERATURE, FieldComponents::Spectral::SCALAR));

      // Get X and Z dimensions
      int nx = 15;
      int nz = 15;

      int nnzLhs = static_cast<int>(0.1*nx*nx);
      int nnzRhs = static_cast<int>(0.1*nx*nx);

      // Storage for the matrix row
      DecoupledZSparse  block = std::make_pair(SparseMatrix(nx*nz, nx*nz),SparseMatrix(nx*nz, nx*nz));
      SparseMatrixZ  sparseLhs(nx*nz*fieldIds.size(), nx*nz*fieldIds.size());
      SparseMatrixZ  sparseRhs(nx*nz*fieldIds.size(), nx*nz*fieldIds.size());

      // Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = 1800.0;
      MHDFloat Pr = 1.0;
      MHDFloat Gamma = .1;
      MHDFloat chi = 45;
      // Wave number
      MHDFloat k_ = 3.0;

      MHDFloat lowRa = 600;
      MHDFloat highRa = 3000;
      Ra = (lowRa+highRa)/2;

      eigenSolver.setWhich("LM");
      eigenSolver.setSigma(MHDComplex(1.0,0.0));

      for(int i = 0; i < 20; ++i)
      {
         std::cerr << "Ra = " << Ra << std::endl;
         SparseMatrix  tmp;
         sparseLhs.reserve(nnzLhs);
         sparseRhs.reserve(nnzRhs);
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

         nnzLhs = sparseLhs.nonZeros();
         nnzRhs = sparseRhs.nonZeros();

         // Factorize eigenvalue problem operators
         eigenSolver.compute(sparseLhs, sparseRhs);

         // Compute eigenvalues
         ArrayZ eigenValues(5*nx+5*nz);
         //MatrixZ eigenVectors(sparseLhs.rows(), eigenValues.size());
         MatrixZ eigenVectors;
         eigenSolver.solve(eigenValues, eigenVectors);

         if(eigenSolver.info() == 0)
         {
            // Show eigenValues
            std::cerr << std::setprecision(16) << "Ra = " << Ra << " EVs = " << eigenValues.transpose() << std::endl;

            if(eigenValues(0).real() > 0)
            {
               if(eigenValues(0).real() < 1e-2)
               {
                  break;
               }
               highRa = Ra;
               Ra = (lowRa+Ra)/2;
            } else
            {
               lowRa = Ra;
               Ra = (highRa+Ra)/2;
            }
         } else
         {
            std::cerr << "ARPACK failed!" << std::endl;
            //throw Exception("ARPACK failed with -14");
            eigenValues.resize(eigenValues.size()*2);
         }
      }

      std::cerr << "Critical Ra = " << Ra << std::endl;
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
