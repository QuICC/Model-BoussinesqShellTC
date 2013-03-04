/** \file IVariableHdf5Reader.hpp
 *  \brief Implementation of a generic variable data file reader
 */

#ifndef IVARIABLEHDF5READER_HPP
#define IVARIABLEHDF5READER_HPP

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/PhysicalNames.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoHdf5/IHdf5Reader.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of a generic variable data file reader
    */
   class IVariableHdf5Reader: public IoHdf5::IHdf5Reader
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name    File name
          * @param ext     File extension
          * @param header  File header
          * @param type    Type string of file
          * @param version Version string of file
          */
         IVariableHdf5Reader(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
          * @brief Destructor
          */
         virtual ~IVariableHdf5Reader();

         /**
          * @brief Add name of expected variable to be added
          */
         void expect(const PhysicalNames::Id id);

         /**
          * @brief Make sure all the expected variables have been added
          */
         bool isFull() const;

         /**
          * @brief Add scalar variable to file
          *
          * @param scalar Scalar variable to add
          */
         void addScalar(const std::pair<PhysicalNames::Id,Datatypes::SharedScalarVariableType>& scalar);

         /**
          * @brief Add vector variable to file
          *
          * @param vector Vector variable to add
          */
         void addVector(const std::pair<PhysicalNames::Id,Datatypes::SharedVectorVariableType>& vector);

         /**
          * @brief Write State to file
          */
         virtual void read() = 0;
         
      protected:
         /// Typedef for the scalar const iterator
         typedef std::map<PhysicalNames::Id,Datatypes::SharedScalarVariableType>::const_iterator  scalar_iterator;

         /// Typedef for the vector const iterator
         typedef std::map<PhysicalNames::Id,Datatypes::SharedVectorVariableType>::const_iterator  vector_iterator;

         /// Typedef for the scalar iterator range
         typedef std::pair<scalar_iterator,scalar_iterator>  scalar_iterator_range;

         /// Typedef for the scalar iterator range
         typedef std::pair<vector_iterator, vector_iterator>  vector_iterator_range;

         /**
          * @brief Resolution information
          */
         SharedResolution mspRes;

         /**
          * @brief Set the read arguments
          */
         void setReadArguments(const Dimensions::Space::Id id);

         /**
          * @brief Read truncation information
          */
         void readTruncation();

         /**
          * @brief Check truncation compatibily between data and file
          */
         void checkTruncation(const Dimensions::Space::Id id);

         /**
          * @brief Get iterator range to scalars
          */
         scalar_iterator_range   scalarRange();

         /**
          * @brief Get iterator range to vectors
          */
         vector_iterator_range   vectorRange();

      private:
         /**
          * @brief Set the resolution and use it for preliminary initialisation
          * 
          * @param spRes      Resolution information
          */
         void setResolution(SharedResolution spRes);

         /**
          * @brief Set the maximum number of IO operations
          */
         void setCollIo();

         /**
          * @brief Simulation resolution based on file parameters
          */
         SharedSimulationResolution mspFileRes;

         /**
          * @brief Storage for the names of the expecte variables
          */
         std::set<PhysicalNames::Id>  mExpected;

         /**
          * @brief Storage for the scalars
          */
         std::map<PhysicalNames::Id,Datatypes::SharedScalarVariableType>   mScalars;

         /**
          * @brief Storage for the vectors
          */
         std::map<PhysicalNames::Id,Datatypes::SharedVectorVariableType>   mVectors;
   };

   /// Typedef for a smart reference counting pointer of a Variable HDF5 numbering reader
   typedef SharedPtrMacro<IVariableHdf5Reader>   SharedIVariableHdf5Reader;

}
}

#endif // IVARIABLEHDF5READER_HPP
