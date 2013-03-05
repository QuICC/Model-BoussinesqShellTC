/** \file ConfigurationWriterTest.cpp
 *  \brief Implementation of test case for ConfigurationWriter
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoConfig/ConfigurationWriter.hpp"
#include "IoConfig/ConfigParts/PhysicalPart.hpp"
#include "IoConfig/ConfigParts/BoundaryPart.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the ConfigurationWriter implementation
    */
   class ConfigurationWriterTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         ConfigurationWriterTest();

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationWriterTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   ConfigurationWriterTest::ConfigurationWriterTest()
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

   ConfigurationWriterTest::~ConfigurationWriterTest()
   {
      // Finalise framework
      FrameworkMacro::finalize();
   }

//   void ConfigurationWriterTest::SetUp()
//   {
//   }

//   void ConfigurationWriterTest::TearDown()
//   {
//   }

   /**
    * @brief Test default 1D setup
    */
   TEST_F(ConfigurationWriterTest, Default1D)
   {
      // Set dimension
      int dim = 1;

      // Set type string
      std::string type = "Test1D";

      // Create configuration writer
      IoConfig::ConfigurationWriter writer(dim, type);

      // Initialise writer
      writer.init();

      // Write configuration
      writer.write();

      // Finalise writer
      writer.finalize();

      ASSERT_TRUE(true);
   }

   /**
    * @brief Test default 2D setup
    */
   TEST_F(ConfigurationWriterTest, Default2D)
   {
      // Set dimension
      int dim = 2;

      // Set type string
      std::string type = "Test2D";

      // Create configuration writer
      IoConfig::ConfigurationWriter writer(dim, type);

      // Initialise writer
      writer.init();

      // Write configuration
      writer.write();

      // Finalise writer
      writer.finalize();

      ASSERT_TRUE(true);
   }

   /**
    * @brief Test default 3D setup
    */
   TEST_F(ConfigurationWriterTest, Default3D)
   {
      // Set dimension
      int dim = 3;

      // Set type string
      std::string type = "Test3D";

      // Create configuration writer
      IoConfig::ConfigurationWriter writer(dim, type);

      // Initialise writer
      writer.init();

      // Write configuration
      writer.write();

      // Finalise writer
      writer.finalize();

      ASSERT_TRUE(true);
   }

   /**
    * @brief Test data setup
    */
   TEST_F(ConfigurationWriterTest, ConfigData)
   {
      // Set dimension
      int dim = 3;

      // Set type string
      std::string type = "TestData";

      // Create configuration writer
      IoConfig::ConfigurationWriter writer(dim, type);

      // Add the physical part
      std::vector<std::string>   names;
      names.push_back("ekman");
      names.push_back("prandtl");
      names.push_back("rayleigh");
      IoConfig::SharedPhysicalPart   spPhys(new IoConfig::PhysicalPart(names));
      writer.addPart(IoConfig::SimulationBlocks::PHYSICAL, spPhys);

      // Add the boundary part
      names.clear();
      names.push_back("temperature");
      names.push_back("velocityz");
      IoConfig::SharedBoundaryPart   spBound(new IoConfig::BoundaryPart(names));
      writer.addPart(IoConfig::SimulationBlocks::BOUNDARY, spBound);

      // Initialise writer
      writer.init();

      // Set truncation information
      writer.rspTruncation()->setIValue("dim1D", 9);
      writer.rspTruncation()->setIValue("dim2D", 11);
      writer.rspTruncation()->setIValue("dim3D", 13);

      // Set parallel information
      writer.rspParallel()->setIValue("cpus", 42);

      // Set timestepping information
      writer.rspTimestepping()->setFValue("time", 1.5);
      writer.rspTimestepping()->setFValue("timestep", 1e-5);

      // Set run information
      writer.rspRun()->setFValue("sim", 3.0);
      writer.rspRun()->setFValue("wall", 5.23);

      // Set IO information
      writer.rspIo()->setFValue("ascii", 100);
      writer.rspIo()->setFValue("hdf5", 501);

      // Set physical information
      writer.rspPhysical()->setFValue("ekman", 1.5e-7);
      writer.rspPhysical()->setFValue("prandtl", 1.1);
      writer.rspPhysical()->setFValue("rayleigh", 1728.23);

      // Set boundary information
      writer.rspBoundary()->setIValue("temperature", 0);
      writer.rspBoundary()->setIValue("velocityz", 2);

      // Write configuration
      writer.write();

      // Finalise writer
      writer.finalize();

      ASSERT_TRUE(true);
   }

}
}

/// Main to execute all test from test case
int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
