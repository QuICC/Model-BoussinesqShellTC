/** 
 * @file StateFileReaderTest.cpp
 * @brief Implementation of test case for StateFileReader
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoVariable/StateFileReader.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "SpatialSchemes/3D/TFTScheme.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StateFileReader implementation
    */
   class StateFileReaderTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         StateFileReaderTest();

         /**
          * @brief Destructor
          */
         virtual ~StateFileReaderTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   StateFileReaderTest::StateFileReaderTest()
   {
   }

   StateFileReaderTest::~StateFileReaderTest()
   {
   }

//   void StateFileReaderTest::SetUp()
//   {
//   }

//   void StateFileReaderTest::TearDown()
//   {
//   }

   /**
    * @brief Test full file
    */
   TEST_F(StateFileReaderTest, FullFile)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 11; dim(1) = 13; dim(2) = 12;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Chect test data
      int i_, j_, k_;
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)));
            }
         }
      }

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
   }

   /**
    * @brief Test truncated file in first dimension
    */
   TEST_F(StateFileReaderTest, Trunc1D)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 8; dim(1) = 13; dim(2) = 12;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Create test data
      int i_, j_, k_;
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)));
            }
         }
      }

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
   }

   /**
    * @brief Test truncated read in second dimension
    */
   TEST_F(StateFileReaderTest, Trunc2D)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 11; dim(1) = 11; dim(2) = 12;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Create test data
      int i_, j_, k_;
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)));
            }
         }
      }

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
   }

   /**
    * @brief Test truncated file read in third dimension
    */
   TEST_F(StateFileReaderTest, Trunc3D)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 11; dim(1) = 13; dim(2) = 7;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Create test data
      int i_, j_, k_;
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)));
            }
         }
      }

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
   }

   /**
    * @brief Test full file with extended first dimension
    */
   TEST_F(StateFileReaderTest, Extended1D)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 16; dim(1) = 13; dim(2) = 12;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Create test data
      int i_, j_, k_;
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               if( i_ < 12)
               {
                  EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)));
               } else
               {
                  EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(0,0));
               }
            }
         }
      }

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
   }

   /**
    * @brief Test full file with extended second dimension
    */
   TEST_F(StateFileReaderTest, Extended2D)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 11; dim(1) = 15; dim(2) = 12;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Create test data
      int i_, j_, k_;
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               if( k_ < 14)
               {
                  EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)));
               } else
               {
                  EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(0,0));
               }
            }
         }
      }

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
   }

   /**
    * @brief Test full file with extended second dimension
    */
   TEST_F(StateFileReaderTest, Extended3D)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      // Set type string
      std::string type = TFTScheme::type();
      bool isRegular = TFTScheme::isRegular();

      // Set spectral and physical dimensions
      ArrayI dim(3); dim(0) = 11; dim(1) = 13; dim(2) = 21;
     
      // Create the load splitter
      Parallel::LoadSplitter splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Initialise the load splitter
      splitter.init<TFTScheme>(dim);

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Store the shared resolution object
      SharedResolution spRes = best.first;

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);

      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Make sure all information was provided
      if(state.isFull())
      {
         // Initialise state file
         state.init();

         // Write state file
         state.read();

         // Finalize state file
         state.finalize();
      } else
      {
         ASSERT_TRUE(false);
      }

      // Create test data
      int i_, j_, k_;
      for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               if( j_ < 13)
               {
                  EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(static_cast<MHDFloat>(-k_),static_cast<MHDFloat>(i_)*0.01+static_cast<MHDFloat>(j_)));
               } else
               {
                  EXPECT_EQ(spScalar->dom(0).perturbation().point(i,j,k), MHDComplex(0,0));
               }
            }
         }
      }

      // Synchronize over CPUs
      FrameworkMacro::synchronize();
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
