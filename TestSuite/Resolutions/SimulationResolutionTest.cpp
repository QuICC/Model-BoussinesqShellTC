/** \file SimulationResolutionTest.cpp
 *  \brief Implementation of test case for SimulationResolution
 */

#include "Resolutions/SimulationResolution.hpp"
#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SimulationResolution implementation
    */
   class SimulationResolutionTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SimulationResolutionTest();

         /**
          * @brief Destructor
          */
         virtual ~SimulationResolutionTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SimulationResolutionTest::SimulationResolutionTest()
   {
   }

   SimulationResolutionTest::~SimulationResolutionTest()
   {
   }

//   void SimulationResolutionTest::SetUp()
//   {
//   }

//   void SimulationResolutionTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(SimulationResolutionTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
