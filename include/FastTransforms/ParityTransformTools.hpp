/**
 * @file ParityTransformTools.hpp
 * @brief Tools for transform that are aware of parity
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PARITYTRANSFORMTOOLS_HPP
#define PARITYTRANSFORMTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Tools for transform that are aware of parity
    */
   class ParityTransformTools
   {
      public:
         /**
          * @brief Extract parity modes from full field
          */
         static void extractParityModes(Matrix& rSelected, const Matrix& data, const MatrixI& info, const int rows);
         
         /**
          * @brief Set parity modes into full field
          */
         static void setParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows);

         /**
          * @brief Add the parity modes into the whole data
          */
         static void addParityModes(Matrix& rData, const Matrix& selected, const MatrixI& info, const int rows);

         /**
          * @brief Scale the parity modes
          */
         static void scaleParityModes(Matrix& rData, const MatrixI& info, const MHDFloat scale, const int rows);

         /**
          * @brief Extract the parity modes from the whole data
          */
         static void extractParityModes(Matrix& rSelected, const MatrixZ& data, const bool isReal, const MatrixI& info, const int rows);

         /**
          * @brief Set the parity modes into the whole data
          */
         static void setParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows);

         /**
          * @brief Add the parity modes into the whole data
          */
         static void addParityModes(MatrixZ& rData, const Matrix& selected, const bool isReal, const MatrixI& info, const int rows);

         /**
          * @brief Scale the parity modes
          */
         static void scaleParityModes(MatrixZ& rData, const bool isReal, const MatrixI& info, const MHDFloat scale, const int rows);

         /**
          * @brief Apply matrix operator to data with proper parity
          */
         static void applyOperator(Matrix& rData, const SparseMatrix& op, const MatrixI& infoPair, const MHDFloat scale, const int rows);

         /**
          * @brief Apply matrix operator to data with proper parity
          */
         static void applyOperator(MatrixZ& rData, const bool isReal, const SparseMatrix& opPair, const MatrixI& infoPair, const MHDFloat scale, const int rows);
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         ParityTransformTools();

         /**
          * @brief Empty Destructor
          */
         ~ParityTransformTools();

   };

}
}

#endif // PARITYTRANSFORMTOOLS_HPP
