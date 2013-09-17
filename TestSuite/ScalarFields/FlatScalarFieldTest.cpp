/** 
 * @file FlatScalarFieldTest.cpp
 * @brief Implementation of test cases for the flat storage scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "gtest/gtest.h"

#include "Framework/FrameworkMacro.h"
#include "ScalarFields/FlatScalarField.hpp"
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

namespace TestSuite {

   /**
    * @brief Test fixture for the FlatScalarField implementation
    */
   class FlatScalarFieldTest : public ::testing::Test {
      public:

      protected:
         /**
          * @brief Constructor
          */
         FlatScalarFieldTest();

         /**
          * @brief Destructor
          */
         virtual ~FlatScalarFieldTest();

         /**
          * @brief Do Set-up work before each test
          */
         //virtual void SetUp() {};

         /**
          * @brief Do tear-down work after each test
          */
         //virtual void TearDown() {};
   };

   FlatScalarFieldTest::FlatScalarFieldTest()
   {
   }

   FlatScalarFieldTest::~FlatScalarFieldTest()
   {
   }

//   void FlatScalarFieldTest::SetUp()
//   {
//   }

//   void FlatScalarFieldTest::TearDown()
//   {
//   }

   /**
    * @brief Test 1D flat storage scalar field 
    *
    * @param FlatScalarFieldTest Test fixture ID
    * @param OneDimensional      Test ID
    */
   TEST_F(FlatScalarFieldTest, OneDimensional)
   {
      // Create 1D ScalarField setup
      int dim1D = 10;
      SharedPtrMacro<Datatypes::ScalarFieldSetup<Dimensions::ONED> >   spSetup(new Datatypes::ScalarFieldSetup<Dimensions::ONED>(dim1D));

      // Create real 1D scalar field
      Datatypes::FlatScalarField<MHDFloat, Dimensions::ONED> dScalar(spSetup);

      // Test size of data storage
      ASSERT_TRUE(dScalar.data().rows() == spSetup->dataRows());
      ASSERT_TRUE(dScalar.data().cols() == spSetup->dataCols());

      // Set test array
      Array dTest = Array::LinSpaced(dim1D, 0, dim1D-1);
      dScalar.setData(dTest);

      // Output scalar array through point
      Array dOut(dim1D);
      for(int i = 0; i < dim1D; ++i)
      {
         dOut(i) = dScalar.point(i);
      }

      // Compare scalar and test array
      ASSERT_TRUE((dTest.array() == dOut.array()).all());

      dScalar.setPoint(-42.0, dim1D/3);

      // Test value set through setPoint
      ASSERT_EQ(dScalar.point(dim1D/3), -42.0);

      //
      // COMPLEX CASE
      //

      // Create complex 1D scalar field
      Datatypes::FlatScalarField<MHDComplex, Dimensions::ONED> zScalar(spSetup);

      // Test size of data storage
      ASSERT_TRUE(zScalar.data().rows() == spSetup->dataRows());
      ASSERT_TRUE(zScalar.data().cols() == spSetup->dataCols());

      // Set test array
      ArrayZ zTest(dim1D);
      zTest.setConstant(MHDComplex(3.0, -2.0));
      zScalar.setData(zTest);

      // Output scalar array through point
      ArrayZ zOut(dim1D);
      for(int i = 0; i < dim1D; ++i)
      {
         zOut(i) = zScalar.point(i);
      }

      // Compare scalar and test array
      ASSERT_TRUE((zTest.array() == zOut.array()).all());

      zScalar.setPoint(-42.0, dim1D/2);

      // Test value set through setPoint
      ASSERT_EQ(zScalar.point(dim1D/2), MHDComplex(-42.0,0.0));
   }

   /**
    * @brief Test 2D flat storage scalar field 
    *
    * @param FlatScalarFieldTest Test fixture ID
    * @param TwoDimensional      Test ID
    */
   TEST_F(FlatScalarFieldTest, TwoDimensional)
   {
      // Create 2D scalar field setup
      int dim2D = 10;
      SharedArrayI   spDim1D(new ArrayI(dim2D));
      spDim1D->setConstant(10);
      SharedPtrMacro<Datatypes::ScalarFieldSetup<Dimensions::TWOD> >   spSetup(new Datatypes::ScalarFieldSetup<Dimensions::TWOD>(spDim1D, 10));

      // Create real 2D scalar field
      Datatypes::FlatScalarField<MHDFloat, Dimensions::TWOD> dScalar(spSetup);

      // Create real 2D scalar field
      // Test size of data storage
      ASSERT_TRUE(dScalar.data().rows() == spSetup->dataRows());
      ASSERT_TRUE(dScalar.data().cols() == spSetup->dataCols());

      // Set test array
      Matrix dTest1 = Matrix::Random(dim2D, dim2D);
      Matrix dTest2 = Matrix::Random(dim2D, dim2D);
      Matrix dTest = dTest1 + dTest2;
      dScalar.setData(dTest1+dTest2);

      Matrix dOut(dim2D, dim2D);
      // Get scalar data
      for(int j = 0; j < dim2D; ++j)
      {
         for(int i = 0; i < (*spDim1D)(j); ++i)
         {
            dOut(i,j) = dScalar.point(i,j);
         }
      }
      // Compare test and output
      ASSERT_TRUE((dTest.array() == dOut.array()).all());

      // Get scalar data
      dOut.setConstant(0.0);
      for(int i = 0; i < dim2D; ++i)
      {
         dOut.col(i) = dScalar.profile(i);
      }
      // Compare test and output
      ASSERT_TRUE((dTest.array() == dOut.array()).all());

      // Test setPoint
      dScalar.setPoint(-42.0, dim2D/3, dim2D/2);
      ASSERT_EQ(dScalar.point(dim2D/3, dim2D/2), -42.0);

      // Test setProfile
      Array dProf = Array::Random(dim2D);
      dScalar.setProfile(dProf, 3);
      ASSERT_TRUE((dScalar.profile(3).array() == dProf.array()).all()); 

      //
      // COMPLEX CASE
      //

      // Create complex 2D scalar field
      Datatypes::FlatScalarField<MHDComplex, Dimensions::TWOD> zScalar(spSetup);

      // Test size of data storage
      ASSERT_TRUE(zScalar.data().rows() == spSetup->dataRows());
      ASSERT_TRUE(zScalar.data().cols() == spSetup->dataCols());

      // Set test array
      MatrixZ zTest1 = MatrixZ::Random(dim2D, dim2D);
      MatrixZ zTest2 = MatrixZ::Random(dim2D, dim2D);
      MatrixZ zTest = zTest1*zTest2;
      zScalar.setData(zTest1*zTest2);

      MatrixZ zOut(dim2D, dim2D);
      // Get scalar data
      for(int j = 0; j < dim2D; ++j)
      {
         for(int i = 0; i < (*spDim1D)(j); ++i)
         {
            zOut(i,j) = zScalar.point(i,j);
         }
      }
      // Compare test and output
      ASSERT_TRUE((zTest.array() == zOut.array()).all());

      // Get scalar data
      zOut.setConstant(0.0);
      for(int i = 0; i < dim2D; ++i)
      {
         zOut.col(i) = zScalar.profile(i);
      }
      // Compare test and output
      ASSERT_TRUE((zTest.array() == zOut.array()).all());

      // Test setPoint
      zScalar.setPoint(MHDComplex(-42.0, 21.0), dim2D/2, dim2D/3);
      ASSERT_EQ(zScalar.point(dim2D/2, dim2D/3), MHDComplex(-42.0,21.0));

      // Test setProfile
      ArrayZ zProf = ArrayZ::Random(dim2D);
      zScalar.setProfile(zProf, dim2D-1);
      ASSERT_TRUE((zScalar.profile(dim2D-1).array() == zProf.array()).all()); 
   }

   /**
    * @brief Test 3D flat storage scalar field
    *
    * @param FlatScalarFieldTest Test fixture ID
    * @param ThreeDimensional    Test ID
    */
   TEST_F(FlatScalarFieldTest, ThreeDimensional)
   {
      // Create 3D scalar field setup
      int dim3D = 10;
      SharedArrayI   spDim1D(new ArrayI(dim3D));
      spDim1D->setConstant(10);
      SharedArrayI   spDim2D(new ArrayI(dim3D));
      spDim2D->setConstant(10);
      SharedPtrMacro<Datatypes::ScalarFieldSetup<Dimensions::THREED> >   spSetup(new Datatypes::ScalarFieldSetup<Dimensions::THREED>(spDim1D, spDim2D, dim3D));

      // Create real 3D scalar field
      Datatypes::FlatScalarField<MHDFloat, Dimensions::THREED> dScalar(spSetup);

      // Test size of data storage
      ASSERT_TRUE(dScalar.data().rows() == spSetup->dataRows());
      ASSERT_TRUE(dScalar.data().cols() == spSetup->dataCols());

      // Set test array
      Matrix dTest1 = Matrix::Random(spSetup->dataRows(), spSetup->dataCols());
      Matrix dTest2 = Matrix::Random(spSetup->dataRows(), spSetup->dataCols());
      Matrix dTest = dTest1 - dTest2;

      dScalar.setData(dTest1-dTest2);
      Matrix dOut(spSetup->dataRows(), spSetup->dataCols());
      // Get scalar data
      for(int k = 0; k < dim3D; ++k)
      {
         for(int j = 0; j < (*spDim2D)(k); ++j)
         {
            for(int i = 0; i < (*spDim1D)(k); ++i)
            {
               dOut(i, spDim2D->head(k).sum() + j) = dScalar.point(i,j, k);
            }
         }
      }
      // Compare test and output
      ASSERT_TRUE((dTest.array() == dOut.array()).all());

      // Get scalar data
      dTest.setRandom();
      dScalar.setData(dTest);
      dOut.setConstant(0.0);
      for(int k = 0; k < dim3D; ++k)
      {
         for(int j = 0; j < (*spDim2D)(k); ++j)
         {
            dOut.col(spDim2D->head(k).sum() + j) = dScalar.profile(j,k);
         }
      }
      // Compare test and output
      ASSERT_TRUE((dTest.array() == dOut.array()).all());

      // Get scalar data
      dTest.setRandom();
      dScalar.setData(dTest);
      dOut.setConstant(0.0);
      for(int k = 0; k < dim3D; ++k)
      {
         dOut.block(0, spDim2D->head(k).sum(), (*spDim1D)(k), (*spDim2D)(k)) = dScalar.slice(k);
      }
      // Compare test and output
      ASSERT_TRUE((dTest.array() == dOut.array()).all());

      // Test setPoint
      dScalar.setPoint(-42.0, (*spDim1D)(dim3D/4)/3, (*spDim2D)(dim3D/4)/2, dim3D/4);
      ASSERT_EQ(dScalar.point((*spDim1D)(dim3D/4)/3, (*spDim2D)(dim3D/4)/2, dim3D/4), -42.0);

      // Test setProfile
      Array dProf = Array::Random((*spDim2D)(dim3D/3));
      dScalar.setProfile(dProf, (*spDim2D)(dim3D/3)/2, dim3D/3);
      ASSERT_TRUE((dScalar.profile((*spDim2D)(dim3D/3)/2, dim3D/3).array() == dProf.array()).all()); 

      // Test setSlice
      Matrix dSlice = Matrix::Random((*spDim1D)(dim3D/2), (*spDim2D)(dim3D/2));
      dScalar.setSlice(dSlice, dim3D/2);
      ASSERT_TRUE((dScalar.slice(dim3D/2).array() == dSlice.array()).all()); 

      //
      // COMPLEX CASE
      //

      // Create complex 3D scalar field
      Datatypes::FlatScalarField<MHDComplex, Dimensions::THREED> zScalar(spSetup);

      // Test size of data storage
      ASSERT_TRUE(zScalar.data().rows() == spSetup->dataRows());
      ASSERT_TRUE(zScalar.data().cols() == spSetup->dataCols());

      // Set test array
      MatrixZ zTest1 = MatrixZ::Random(spSetup->dataRows(), spSetup->dataCols());
      MatrixZ zTest2 = MatrixZ::Random(spSetup->dataRows(), spSetup->dataCols());
      MatrixZ zTest = zTest1 - zTest2;

      zScalar.setData(zTest1-zTest2);
      MatrixZ zOut(spSetup->dataRows(), spSetup->dataCols());
      // Get scalar data
      for(int k = 0; k < dim3D; ++k)
      {
         for(int j = 0; j < (*spDim2D)(k); ++j)
         {
            for(int i = 0; i < (*spDim1D)(k); ++i)
            {
               zOut(i, spDim2D->head(k).sum() + j) = zScalar.point(i,j, k);
            }
         }
      }
      // Compare test and output
      ASSERT_TRUE((zTest.array() == zOut.array()).all());

      // Get scalar data
      zTest.setRandom();
      zScalar.setData(zTest);
      zOut.setConstant(0.0);
      for(int k = 0; k < dim3D; ++k)
      {
         for(int j = 0; j < (*spDim2D)(k); ++j)
         {
            zOut.col(spDim2D->head(k).sum() + j) = zScalar.profile(j,k);
         }
      }
      // Compare test and output
      ASSERT_TRUE((zTest.array() == zOut.array()).all());

      // Get scalar data
      zTest.setRandom();
      zScalar.setData(zTest);
      zOut.setConstant(0.0);
      for(int k = 0; k < dim3D; ++k)
      {
         zOut.block(0, spDim2D->head(k).sum(), (*spDim1D)(k), (*spDim2D)(k)) = zScalar.slice(k);
      }
      // Compare test and output
      ASSERT_TRUE((zTest.array() == zOut.array()).all());

      // Test setPoint
      zScalar.setPoint(MHDComplex(21.2, -42.0), (*spDim1D)(dim3D/4)/3, (*spDim2D)(dim3D/4)/2, dim3D/4);
      ASSERT_EQ(zScalar.point((*spDim1D)(dim3D/4)/3, (*spDim2D)(dim3D/4)/2, dim3D/4), MHDComplex(21.2,-42.0));

      // Test setProfile
      ArrayZ zProf = ArrayZ::Random((*spDim2D)(dim3D/3));
      zScalar.setProfile(zProf, (*spDim2D)(dim3D/3)/2, dim3D/3);
      ASSERT_TRUE((zScalar.profile((*spDim2D)(dim3D/3)/2, dim3D/3).array() == zProf.array()).all()); 

      // Test setSlice
      MatrixZ zSlice = MatrixZ::Random((*spDim1D)(dim3D/2), (*spDim2D)(dim3D/2));
      zScalar.setSlice(zSlice, dim3D/2);
      ASSERT_TRUE((zScalar.slice(dim3D/2).array() == zSlice.array()).all()); 
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
