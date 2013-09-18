/** 
 * @file TimerTest.cpp
 * @brief Implementation of test case for a generic timer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include <unistd.h>

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "Timers/TimerMacro.h"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for a generic timer implementation
    */
   class TimerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         TimerTest();

         /**
          * @brief Destructor
          */
         virtual ~TimerTest() {};

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};

         /// Number of milliseconds to wait
         int mMSec;

         /// Wait timespec
         timespec mWait;

         /// Remaining timespec
         timespec mRem;
   };

   TimerTest::TimerTest()
   {
      this->mMSec = 500;

      // Setup timespecs to wait for 500 milliseconds
      this->mWait.tv_sec = 0;
      this->mWait.tv_nsec = this->mMSec*1000000;

      // Init
      this->mRem.tv_sec = 0;
      this->mRem.tv_nsec = 0;
   }

   /**
    * @brief Test default constructor 
    *
    * @param TimerTest     Test fixture ID
    * @param DefaultConstructor  Test ID
    */
   TEST_F(TimerTest, DefaultConstructor) {
      // Create timer with default constructor
      TimerMacro timer;

      // Wait without starting the timer
      nanosleep(&this->mWait, &this->mRem);
      // Expect zero time in timer
      EXPECT_EQ(0.0, timer.time());

      // Time and wait
      timer.start();
      nanosleep(&this->mWait, &this->mRem);
      timer.stop();
      // Expect timer to be accurate to 1 milliseconds
      EXPECT_EQ(this->mMSec, static_cast<int>(1000.*timer.time()));
   }

   /**
    * @brief Test explicit noautostart 
    *
    * @param TimerTest  Test fixture ID
    * @param NoAutostart      Test ID
    */
   TEST_F(TimerTest, NoAutostart) {
      // Create timer with no autostart
      TimerMacro timer(false);

      // Wait without starting the timer
      nanosleep(&this->mWait, &this->mRem);
      // Expect zero time in timer
      EXPECT_EQ(0.0, timer.time());

      // Time and wait
      timer.start();
      nanosleep(&this->mWait, &this->mRem);
      timer.stop();
      // Expect timer to be accurate to 1 milliseconds
      EXPECT_EQ(this->mMSec, static_cast<int>(1000.*timer.time()));
   }

   /**
    * @brief Test autostart 
    *
    * @param TimerTest  Test fixture ID
    * @param Autostart        Test ID
    */
   TEST_F(TimerTest, Autostart) {
      // Create timer with autostart
      TimerMacro timer(true);

      // Wait without starting the timer
      nanosleep(&this->mWait, &this->mRem);
      timer.stop();
      // Expect timer to be accurate to 1 milliseconds
      EXPECT_EQ(this->mMSec, static_cast<int>(1000.*timer.time()));

      // Time and wait
      timer.start();
      nanosleep(&this->mWait, &this->mRem);
      timer.stop();
      // Expect timer to be accurate to 1 milliseconds
      EXPECT_EQ(this->mMSec, static_cast<int>(1000.*timer.time()));
   }

}
}

/**
 * @brief Main function to execute all test cases
 *
 * @param argc Number of arguments
 * @param argv Arguments
 */
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
