/** 
 * @file ScalarFieldSetupTest.cpp
 * @brief Implementation of test cases for the scalar field setup object
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "Enums/Dimensions.hpp"
#include "ScalarFields/ScalarFieldSetup.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ScalarFieldSetup implementation
    */
   class ScalarFieldSetupTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ScalarFieldSetupTest();

         /**
          * @brief Destructor
          */
         virtual ~ScalarFieldSetupTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ScalarFieldSetupTest::ScalarFieldSetupTest()
   {
   }

   ScalarFieldSetupTest::~ScalarFieldSetupTest()
   {
   }

//   void ScalarFieldSetupTest::SetUp()
//   {
//   }

//   void ScalarFieldSetupTest::TearDown()
//   {
//   }

   /**
    * @brief Test 1D 
    */
   /**
    * @brief Test 1D scalar field 
    *
    * @param ScalarFieldSetupTest   Test fixture ID
    * @param OneDimensional         Test ID
    */
   TEST_F(ScalarFieldSetupTest, OneDimensional)
   {
      // Create 1D scalar field setup
      int dim1D = 10;
      Datatypes::ScalarFieldSetup<Dimensions::ONED>  setup(dim1D);

      // Check the data sizes
      ASSERT_EQ(setup.dataRows(), dim1D);
      ASSERT_EQ(setup.dataCols(), 1);
   }

   /**
    * @brief Test 2D scalar field
    *
    * @param ScalarFieldSetupTest   Test fixture ID
    * @param TwoDimensional         Test ID
    */
   TEST_F(ScalarFieldSetupTest, TwoDimensional)
   {
      // Create 2D scalar field setup
      int dim2D = 10;
      SharedArrayI   spDim1D(new ArrayI(dim2D));
      spDim1D->setConstant(12);
      Datatypes::ScalarFieldSetup<Dimensions::TWOD>  setup(spDim1D, dim2D);

      // Check the data sizes
      ASSERT_EQ(setup.dataRows(), spDim1D->maxCoeff());
      ASSERT_EQ(setup.dataCols(), dim2D);

      // Check column indexes
      for(int i = 0; i < dim2D; ++i)
      {
         ASSERT_EQ(setup.colIdx(i), i);
      }
   }

   /**
    * @brief Test 3D scalar field 
    *
    * @param ScalarFieldSetupTest   Test fixture ID
    * @param ThreeDimensional       Test ID
    */
   TEST_F(ScalarFieldSetupTest, ThreeDimensional)
   {
      // Create 3D scalar field setup
      int dim3D = 10;
      SharedArrayI   spDim2D(new ArrayI(dim3D));
      spDim2D->setConstant(12);
      SharedArrayI   spDim1D(new ArrayI(dim3D));
      spDim1D->setConstant(15);
      Datatypes::ScalarFieldSetup<Dimensions::THREED>  setup(spDim1D, spDim2D, dim3D);

      // Check the data sizes
      ASSERT_EQ(setup.dataRows(), spDim1D->maxCoeff());
      ASSERT_EQ(setup.dataCols(), spDim2D->sum());

      // Check column indexes
      int tot = 0;
      for(int j = 0; j < dim3D; ++j)
      {
         for(int i = 0; i < (*spDim2D)(j); ++i)
         {
            ASSERT_EQ(setup.colIdx(i,j), tot);
            tot++;
         }
      }

      // Check block information
      tot = 0;
      for(int i = 0; i < dim3D; ++i)
      {
         ASSERT_EQ(setup.blockIdx(i), tot);
         ASSERT_EQ(setup.blockRows(i), (*spDim1D)(i));
         ASSERT_EQ(setup.blockCols(i), (*spDim2D)(i));
         tot += (*spDim2D)(i);
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
