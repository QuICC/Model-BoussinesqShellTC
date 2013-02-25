/** \file SimulationIoControlTest.cpp
 *  \brief Implementation of test case for SimulationIoControl
 */

#include "gtest/gtest.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SimulationIoControl implementation
    */
   class SimulationIoControlTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SimulationIoControlTest();

         /**
          * @brief Destructor
          */
         virtual ~SimulationIoControlTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   SimulationIoControlTest::SimulationIoControlTest()
   {
   }

   SimulationIoControlTest::~SimulationIoControlTest()
   {
   }

//   void SimulationIoControlTest::SetUp()
//   {
//   }

//   void SimulationIoControlTest::TearDown()
//   {
//   }

   /**
    * @brief Test default constructor
    */
   TEST_F(SimulationIoControlTest, Constructor) {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
