/** 
 * @file SerialTimerTest.cpp
 * @brief Implementation of test case for SerialTimer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#include "Timers/SerialTimer.hpp"
#include "gtest/gtest.h"
#include <unistd.h>

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SerialTimer implementation
    */
   class SerialTimerTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SerialTimerTest();

         /**
          * @brief Destructor
          */
         //virtual ~SerialTimerTest() {};

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

   SerialTimerTest::SerialTimerTest()
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
    */
   TEST_F(SerialTimerTest, DefaultConstructor) {
      // Create timer with default constructor
      SerialTimer timer;

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
    * @brief Test with explicit noautostart
    */
   TEST_F(SerialTimerTest, NoAutostart) {
      // Create timer with no autostart
      SerialTimer timer(false);

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
    * @brief Test the autostart feature
    */
   TEST_F(SerialTimerTest, Autostart) {
      // Create timer with autostart
      SerialTimer timer(true);

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

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
