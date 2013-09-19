/** 
 * @file SpatialSchemeTest.cpp
 * @brief Implementation of test cases for a spatial scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "LoadSplitter/LoadSplitter.hpp"
#include "TypeSelectors/SpatialSchemeSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/ParallelSelector.hpp"
#include "TransformGroupers/IForwardGrouper.hpp"
#include "TransformGroupers/IBackwardGrouper.hpp"
#include "TransformCoordinators/TransformCoordinatorTools.hpp"
#include "Variables/RequirementTools.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the SpatialScheme implementations
    */
   class SpatialSchemeTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         SpatialSchemeTest();

         /**
          * @brief Destructor
          */
         virtual ~SpatialSchemeTest();

         /**
          * @brief Do Set-up work before each test
          */
         virtual void SetUp();

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown();
         
      private:
         /**
          * @brief Shared resolution 
          */
         SharedResolution mspRes;

         /**
          * @brief Transform coordinator 
          */
         Transform::TransformCoordinatorType   mCoord;

         /**
          * @brief Storage for a shared forward transform grouper
          */
         Transform::SharedIForwardGrouper   mspFwdGrouper;

         /**
          * @brief Storage for a shared backward transform grouper
          */
         Transform::SharedIBackwardGrouper   mspBwdGrouper;

         /**
          * @brief Storage for scalar equations
          */
         std::vector<Equations::SharedIScalarEquation> mScalarEquations;

         /**
          * @brief Storage for vector equations
          */
         std::vector<Equations::SharedIVectorEquation> mVectorEquations;

         /**
          * @brief Map between name and pointer for the scalar variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>  mScalarVariables;

         /**
          * @brief Map between name and pointer for the vector variables
          */
         std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>  mVectorVariables;
   };

   SpatialSchemeTest::SpatialSchemeTest()
   {
   }

   SpatialSchemeTest::~SpatialSchemeTest()
   {
   }

   void SpatialSchemeTest::SetUp()
   {
      // Storage for the variable requirements
      VariableRequirement varInfo;

      // Initialise the variables and set general variable requirements
      RequirementTools::initVariables(varInfo, this->mScalarVariables, this->mVectorVariables, this->mScalarEquations, this->mVectorEquations, this->mspRes);

      // Storage for the nonlinear requirement info
      std::set<PhysicalNames::Id>   nonInfo;

      // Map variables to the equations and set nonlinear requirements
      RequirementTools::mapEquationVariables(nonInfo, this->mScalarEquations, this->mVectorEquations, this->mScalarVariables, this->mVectorVariables);

      // Create load splitter
      Parallel::LoadSplitter   splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Create test resolution
      ArrayI dims(3);
      dims(0) = 64; dims(1) = 64; dims(2) = 64;

      // Initialise the load splitter with spatial scheme
      splitter.init<Schemes::SpatialType>(dims);

      // Get best splitted resolution
      this->mspRes = splitter.bestSplitting().first;

      // Set additional options on final resolution object
      Schemes::SpatialType::tuneResolution(this->mspRes);

      // Initialise the transform grouper
      Parallel::setGrouper(splitter.bestSplitting().second, this->mspFwdGrouper, this->mspBwdGrouper);

      // Set the transform options
      std::map<NonDimensional::Id,MHDFloat> runOptions;

      // Initialise the transform coordinator
      Transform::TransformCoordinatorTools::init(this->mCoord, this->mspFwdGrouper, this->mspBwdGrouper, varInfo, nonInfo, this->mspRes, runOptions);
   }

//   void SpatialSchemeTest::TearDown()
//   {
//   }

   /**
    * @brief Forward transform loop with exact solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, FwdLoopExact)
   {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

   /**
    * @brief Forward transform loop with full resolution solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, FwdLoopFull)
   {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

   /**
    * @brief Backward transform loop with exact solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, BwdLoopExact)
   {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

   /**
    * @brief Backward transform loop with full resolution solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, BwdLoopFull)
   {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

   /**
    * @brief Error accumulation through many loops with full resolution solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, Accumulation)
   {
      ASSERT_TRUE(false) << "##########################################" << std::endl << "## Tests have not yet been implemented! ##" << std::endl << "##########################################";
   }

}
}

/**
 * @brief Main function to execute all test cases
 *
 * @param argc Number of arguments
 * @param argv Arguments
 */
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
