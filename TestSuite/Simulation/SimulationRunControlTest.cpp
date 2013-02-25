/** \file SimulationRunControlTest.cpp
 *  \brief Implementation of test case for SimulationRunControl
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SimulationRunControl implementation
    */
   class SimulationRunControlTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SimulationRunControlTest();

         /**
          * @brief Destructor
          */
         virtual ~SimulationRunControlTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SimulationRunControlTest::SimulationRunControlTest()
   {
   }

   SimulationRunControlTest::~SimulationRunControlTest()
   {
   }

//   void SimulationRunControlTest::SetUp()
//   {
//   }

//   void SimulationRunControlTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(SimulationRunControlTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
