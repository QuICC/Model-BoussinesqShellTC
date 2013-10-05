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
          * @brief Create shared resolution
          */
         SharedResolution createResolution(const int dim1D, const int dim2D, const int dim3D);

         /**
          * @brief Create resolution ratio compared to file resolution
          *
          * @param ratio1D First dimension ratio
          * @param ratio2D Second dimension ratio
          * @param ratio3D Third dimension ratio
          */
         ArrayI resolutionRatio(const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D);

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
         virtual void SetUp();

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};

         /**
          * @brief Shared file resolution
          */
         SharedSimulationResolution mspFileRes;
   };

   StateFileReaderTest::StateFileReaderTest()
   {
   }

   StateFileReaderTest::~StateFileReaderTest()
   {
   }

   void StateFileReaderTest::SetUp()
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      SharedResolution spRes = this->createResolution(10, 10, 10);
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

      // Create state file
      IoVariable::StateFileReader   state("0000", type, isRegular);
      state.expect(PhysicalNames::TEMPERATURE);
      state.addScalar(std::make_pair(PhysicalNames::TEMPERATURE, spScalar));

      // Initialise state file
      state.init();
      this->mspFileRes = state.getFileTruncation();
      state.finalize();
   }

//   void StateFileReaderTest::TearDown()
//   {
//   }

   ArrayI StateFileReaderTest::resolutionRatio(const MHDFloat ratio1D, const MHDFloat ratio2D, const MHDFloat ratio3D)
   {
      // Create test resolution
      ArrayI dims(3);
      dims(0) = ratio1D*(this->mspFileRes->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL)-1);
      dims(1) = ratio2D*(this->mspFileRes->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL)-1);
      dims(2) = ratio3D*(this->mspFileRes->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::SPECTRAL)-1);

      return dims;
   }

   SharedResolution StateFileReaderTest::createResolution(const int dim1D, const int dim2D, const int dim3D)
   {
      // Create load splitter
      Parallel::LoadSplitter   splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Create test resolution
      ArrayI dims(3);
      dims(0) = dim1D; dims(1) = dim2D; dims(2) = dim3D;

      // Initialise the load splitter with spatial scheme
      splitter.init<Schemes::SpatialSelector>(dims);

      // Get best splitting
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Get resolution
      SharedResolution spRes = best.first;

      // Set additional options on final resolution object
      Schemes::SpatialSelector::tuneResolution(spRes);

      return spRes;
   }

   /**
    * @brief Test full file
    */
   TEST_F(StateFileReaderTest, FullFile)
   {
      // Synchronize over CPUs
      FrameworkMacro::synchronize();

      ArrayI dims = this->resolutionRatio(1.0, 1.0, 1.0);
      SharedResolution spRes = this->createResolution(dims(0), dims(1), dims(2));

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

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
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM); ++i)
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

      ArrayI dims = this->resolutionRatio(2.0/3.0, 1.0, 1.0);
      SharedResolution spRes = this->createResolution(dims(0), dims(1), dims(2));

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

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
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM); ++i)
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

      ArrayI dims = this->resolutionRatio(1.0, 2.0/3.0, 1.0);
      SharedResolution spRes = this->createResolution(dims(0), dims(1), dims(2));

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

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
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM); ++i)
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

      ArrayI dims = this->resolutionRatio(1.0, 1.0, 2.0/3.0);
      SharedResolution spRes = this->createResolution(dims(0), dims(1), dims(2));

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

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
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM); ++i)
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

      ArrayI dims = this->resolutionRatio(3.0/2.0, 1.0, 1.0);
      SharedResolution spRes = this->createResolution(dims(0), dims(1), dims(2));

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

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
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               if( i_ < this->mspFileRes->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM))
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

      ArrayI dims = this->resolutionRatio(1.0, 3.0/2.0, 1.0);
      SharedResolution spRes = this->createResolution(dims(0), dims(1), dims(2));

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

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
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               if( k_ < this->mspFileRes->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::TRANSFORM))
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

      ArrayI dims = this->resolutionRatio(1.0, 1.0, 3.0/2.0);
      SharedResolution spRes = this->createResolution(dims(0), dims(1), dims(2));

      // Create scalar variable
      Datatypes::SharedScalarVariableType   spScalar(new Datatypes::ScalarVariableType(spRes));

      // Set type string
      std::string type = Schemes::SpatialSelector::type();
      bool isRegular = Schemes::SpatialSelector::isRegular();

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
            for(int i = 0; i < spRes->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM); ++i)
            {
               i_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATB1D>(i,k);
               if( j_ < this->mspFileRes->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::TRANSFORM))
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
