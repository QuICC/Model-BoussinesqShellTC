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
#include "Equations/Tests/TestSpatialSchemeForwardScalar.hpp"
#include "Equations/Tests/TestSpatialSchemeBackwardScalar.hpp"

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
         
         /**
          * @brief Initialise transform coordinator
          */
         void initCoord();
         
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

         /**
          * @brief Acceptable absolute error
          */
         double mError;

         /**
          * @brief Number of accumulation loops
          */
         int mNAccumul;

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
   };

   SpatialSchemeTest::SpatialSchemeTest()
      //: mError(1e-12), mNAccumul(10), mMax1D(63), mMax2D(57), mMax3D(64)
      : mError(1e-10), mNAccumul(10), mMax1D(13), mMax2D(15), mMax3D(17)
   {
   }

   SpatialSchemeTest::~SpatialSchemeTest()
   {
   }

   void SpatialSchemeTest::SetUp()
   {
      // Create load splitter
      Parallel::LoadSplitter   splitter(FrameworkMacro::id(), FrameworkMacro::nCpu());

      // Create test resolution
      ArrayI dims(3);
      dims(0) = this->mMax1D; dims(1) = this->mMax2D; dims(2) = this->mMax3D;

      // Initialise the load splitter with spatial scheme
      splitter.init<Schemes::SpatialSelector>(dims);

      // Get best splitting
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

      // Get resolution
      this->mspRes = best.first;

      // Set additional options on final resolution object
      Schemes::SpatialSelector::tuneResolution(this->mspRes);

      // Initialise the transform grouper
      Parallel::setGrouper(best.second, this->mspFwdGrouper, this->mspBwdGrouper);
   }

//   void SpatialSchemeTest::TearDown()
//   {
//   }

   void SpatialSchemeTest::initCoord()
   {
      // Storage for the variable requirements
      VariableRequirement varInfo;

      // Storage for the nonlinear requirement info
      std::set<PhysicalNames::Id>   nonInfo;

      // Initialise the variables and set general variable requirements
      RequirementTools::initVariables(varInfo, this->mScalarVariables, this->mVectorVariables, this->mScalarEquations, this->mVectorEquations, this->mspRes);

      // Map variables to the equations and set nonlinear requirements
      RequirementTools::mapEquationVariables(nonInfo, this->mScalarEquations, this->mVectorEquations, this->mScalarVariables, this->mVectorVariables);

      // Set the transform options
      std::map<NonDimensional::Id,MHDFloat> runOptions;

      // Initialise the transform coordinator
      Transform::TransformCoordinatorTools::init(this->mCoord, this->mspFwdGrouper, this->mspBwdGrouper, varInfo, nonInfo, this->mspRes, runOptions);
   }

   /**
    * @brief Test the resolution object
    *
    * This test will fail in parallel setup, but it is difficult to setup correctly (depends on parallelisation algorithm)
    *
    * @param SpatialSchemeTest   Test fixture ID
    * @param Placeholder         Test ID
    */
   TEST_F(SpatialSchemeTest, Resolution)
   {
      EXPECT_EQ(this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL),this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>());
      EXPECT_EQ(this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::PHYSICAL),this->mspRes->cpu()->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATF1D>());
      EXPECT_EQ(this->mspRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::PHYSICAL),this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATF1D>());

      //EXPECT_EQ(this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM),this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>());
      EXPECT_EQ(this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::TRANSFORM),this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>());
      EXPECT_EQ(this->mspRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::TRANSFORM),this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>());
   }

   /**
    * @brief Forward transform loop with exact zero solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, FwdLoopZero)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeForwardScalar>   spEq(new Equations::TestSpatialSchemeForwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeForwardScalar::ZERO);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      // Synchronize  framework
      FrameworkMacro::synchronize();

      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);

      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
      // Check gradient solution
      EXPECT_LT(spEq->computeGradientError(), this->mError);
   }

   /**
    * @brief Forward transform loop with exact constant solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, FwdLoopConstant)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeForwardScalar>   spEq(new Equations::TestSpatialSchemeForwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeForwardScalar::CONSTANT);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      // Synchronize  framework
      FrameworkMacro::synchronize();

      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);

      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
      // Check gradient solution
      EXPECT_LT(spEq->computeGradientError(), this->mError);
   }

   /**
    * @brief Forward transform loop with exact solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, FwdLoopExact)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeForwardScalar>   spEq(new Equations::TestSpatialSchemeForwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeForwardScalar::EXACT);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      // Synchronize  framework
      FrameworkMacro::synchronize();

      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);

      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
      // Check gradient solution
      EXPECT_LT(spEq->computeGradientError(), this->mError);
   }

   /**
    * @brief Backward transform loop with exact zero solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, BwdLoopZero)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::ZERO);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      // Synchronize  framework
      FrameworkMacro::synchronize();

      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
   }

   /**
    * @brief Backward transform loop with exact constant solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, BwdLoopConstant)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::CONSTANT);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      // Synchronize  framework
      FrameworkMacro::synchronize();

      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
   }

   /**
    * @brief Backward transform loop with exact solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, BwdLoopExact)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::EXACT);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      // Synchronize  framework
      FrameworkMacro::synchronize();

      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
   }

   /**
    * @brief Backward transform loop with exact full resolution solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, BwdLoopFull)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::FULL);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      // Synchronize  framework
      FrameworkMacro::synchronize();

      this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

      this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
   }

   /**
    * @brief Error accumulation through many loops with zero solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, AccumulationZero)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::ZERO);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      for(int i = 0; i < this->mNAccumul; i++)
      {
         // Synchronize  framework
         FrameworkMacro::synchronize();

         this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

         this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);
      }

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
   }

   /**
    * @brief Error accumulation through many loops with constant solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, AccumulationConstant)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::CONSTANT);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      for(int i = 0; i < this->mNAccumul; i++)
      {
         // Synchronize  framework
         FrameworkMacro::synchronize();

         this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

         this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);
      }

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
   }

   /**
    * @brief Error accumulation through many loops with exact solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, AccumulationExact)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::EXACT);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      for(int i = 0; i < this->mNAccumul; i++)
      {
         // Synchronize  framework
         FrameworkMacro::synchronize();

         this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

         this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);
      }

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
   }

   /**
    * @brief Error accumulation through many loops with full resolution solution test
    *
    * @param SpatialSchemeTest Test fixture ID
    * @param Placeholder   Test ID
    */
   TEST_F(SpatialSchemeTest, AccumulationFull)
   {
      // Set correct equation
      SharedPtrMacro<Equations::TestSpatialSchemeBackwardScalar>   spEq(new Equations::TestSpatialSchemeBackwardScalar());
      spEq->setIdentity(PhysicalNames::TEMPERATURE);
      spEq->setSolutionType(Equations::TestSpatialSchemeBackwardScalar::FULL);
      this->mScalarEquations.push_back(spEq);

      // Initialise transform coordinator
      this->initCoord();

      for(int i = 0; i < this->mNAccumul; i++)
      {
         // Synchronize  framework
         FrameworkMacro::synchronize();

         this->mspBwdGrouper->transform(this->mScalarVariables, this->mVectorVariables, this->mCoord);

         this->mspFwdGrouper->transform(this->mScalarEquations, this->mVectorEquations, this->mCoord);
      }

      // Synchronize  framework
      FrameworkMacro::synchronize();

      // Check scalar solution
      EXPECT_LT(spEq->computeScalarError(), this->mError);
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
