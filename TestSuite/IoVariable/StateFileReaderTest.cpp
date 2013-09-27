/** 
 * @file StateFileReaderTest.cpp
 * @brief Implementation of test case for StateFileReader
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "IoVariable/StateFileReader.hpp"
#include "LoadSplitter/LoadSplitter.hpp"
#include "TypeSelectors/SpatialSchemeSelector.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the StateFileReader implementation
    */
   class StateFileReaderTest : public ::testing::Test {
      public:
         /**
          * @brief Setup the shared resolution
          */
         void setupResolution(const int dim1D, const int dim2D, const int dim3D);

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

         /**
          * @brief Acceptable absolute error
          */
         double mError;

         /**
          * @brief Maximum order in 1D
          */
         int mMax1D;

         /**
          * @brief Maximum order in 2D
          */
         int mMax2D;

         /**
          * @brief Maximum order in 3D
          */
         int mMax3D;
         
         /**
          * @brief Shared resolution 
          */
         SharedResolution mspRes;
   };

   StateFileReaderTest::StateFileReaderTest()
      : mError(1e-10), mMax1D(13), mMax2D(15), mMax3D(17)
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

   void StateFileReaderTest::setupResolution(const int dim1D, const int dim2D, const int dim3D)
   {
      // Create load splitter
      Parallel::LoadSplitter   splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Create test resolution
      ArrayI dims(3);
      dims(0) = dim1D; dims(1) = dim2D; dims(2) = dim3D;

      // Initialise the load splitter with spatial scheme
      splitter.init<Schemes::SpatialType>(dims);

      // Get best splitting
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Get resolution
      this->mspRes = best.first;

      // Set additional options on final resolution object
      Schemes::SpatialType::tuneResolution(this->mspRes);
   }

   /**
    * @brief Test full file
    */
   TEST_F(StateFileReaderTest, FullFile)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      this->setupResolution(this->mMax1D, this->mMax2D, this->mMax3D);

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Set type string
      std::string type = Schemes::SpatialType::type();
      bool isRegular = Schemes::SpatialType::isRegular();

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
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
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

      this->setupResolution(2*this->mMax1D/3, this->mMax2D, this->mMax3D);

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Set type string
      std::string type = Schemes::SpatialType::type();
      bool isRegular = Schemes::SpatialType::isRegular();

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
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
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

      this->setupResolution(this->mMax1D, 2*this->mMax2D/3, this->mMax3D);

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Set type string
      std::string type = Schemes::SpatialType::type();
      bool isRegular = Schemes::SpatialType::isRegular();

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
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
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

      this->setupResolution(this->mMax1D, this->mMax2D, 2*this->mMax3D/3);

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Set type string
      std::string type = Schemes::SpatialType::type();
      bool isRegular = Schemes::SpatialType::isRegular();

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
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
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

      this->setupResolution(3*this->mMax1D/2, this->mMax2D, this->mMax3D);

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Set type string
      std::string type = Schemes::SpatialType::type();
      bool isRegular = Schemes::SpatialType::isRegular();

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
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
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

      this->setupResolution(this->mMax1D, 3*this->mMax2D/2, this->mMax3D);

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Set type string
      std::string type = Schemes::SpatialType::type();
      bool isRegular = Schemes::SpatialType::isRegular();

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
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
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

      this->setupResolution(this->mMax1D, this->mMax2D, 3*this->mMax3D/2);

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(this->mspRes));

      // Set type string
      std::string type = Schemes::SpatialType::type();
      bool isRegular = Schemes::SpatialType::isRegular();

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
      for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
      {
         k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
         for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); ++j)
         {
            j_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
            for(int i = 0; i < this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL); ++i)
            {
               i_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
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
