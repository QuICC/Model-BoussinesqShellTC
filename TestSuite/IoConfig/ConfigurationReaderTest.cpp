/** 
 * @file ConfigurationReaderTest.cpp
 * @brief Implementation of test case for ConfigurationReader
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoConfig/ConfigurationReader.hpp"
#include "IoConfig/ConfigParts/PhysicalPart.hpp"
#include "IoConfig/ConfigParts/BoundaryPart.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ConfigurationReader implementation
    */
   class ConfigurationReaderTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ConfigurationReaderTest();

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationReaderTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ConfigurationReaderTest::ConfigurationReaderTest()
   {
      // Initilise framework
      FrameworkMacro::init();

      // Set nCpu for serial run
      int nCpu = 1;

      // Set ID and nCpu in MPI case
      #ifdef GEOMHDISCC_MPI
         // Get MPI size
         int size;
         MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
      #endif //GEOMHDISCC_MPI

      // Setup framework
      FrameworkMacro::setup(nCpu);
   }

   ConfigurationReaderTest::~ConfigurationReaderTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void ConfigurationReaderTest::SetUp()
//   {
//   }

//   void ConfigurationReaderTest::TearDown()
//   {
//   }

   /**
    * @brief Test data setup
    */
   TEST_F(ConfigurationReaderTest, ConfigData)
   {
      // Set dimension
      int dim = 3;

      // Set type string
      std::string type = "TestData";

      // Create configuration writer
      IoConfig::ConfigurationReader reader(dim, type);

      // Add the physical part
      std::vector<std::string>   names;
      names.push_back("ekman");
      names.push_back("prandtl");
      names.push_back("rayleigh");
      IoConfig::SharedPhysicalPart   spPhys(new IoConfig::PhysicalPart(names));
      reader.addPart(IoConfig::SimulationBlocks::PHYSICAL, spPhys);

      // Add the boundary part
      names.clear();
      names.push_back("temperature");
      names.push_back("velocityz");
      IoConfig::SharedBoundaryPart   spBound(new IoConfig::BoundaryPart(names));
      reader.addPart(IoConfig::SimulationBlocks::BOUNDARY, spBound);

      // Initialise reader
      reader.init();

      // Read configuration
      reader.read();

      // Finalise reader
      reader.finalize();

      // Check truncation data read from file
      EXPECT_EQ(reader.spTruncation()->iValue("dim1D"), 9);
      EXPECT_EQ(reader.spTruncation()->iValue("dim2D"), 11);
      EXPECT_EQ(reader.spTruncation()->iValue("dim3D"), 13);

      // Check parallel data read from file
      EXPECT_EQ(reader.spParallel()->iValue("cpus"), 42);

      // Check timestepping data read from file
      EXPECT_EQ(reader.spTimestepping()->fValue("time"), 1.5);
      EXPECT_EQ(reader.spTimestepping()->fValue("timestep"), 1e-5);

      // Check run data read from file
      EXPECT_EQ(reader.spRun()->fValue("sim"), 3.0);
      EXPECT_EQ(reader.spRun()->fValue("wall"), 5.23);

      // Check IO data read from file
      EXPECT_EQ(reader.spIo()->fValue("ascii"), 100);
      EXPECT_EQ(reader.spIo()->fValue("hdf5"), 501);

      // Check physical data read from file
      EXPECT_EQ(reader.spPhysical()->fValue("prandtl"), 1.1);
      EXPECT_EQ(reader.spPhysical()->fValue("rayleigh"), 1728.23);
      EXPECT_EQ(reader.spPhysical()->fValue("ekman"), 1.5e-7);

      // Check boundary data read from file
      EXPECT_EQ(reader.spBoundary()->iValue("temperature"), 0);
      EXPECT_EQ(reader.spBoundary()->iValue("velocityz"), 2);
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
