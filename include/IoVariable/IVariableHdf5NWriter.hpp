/** \file IVariableHdf5NWriter.hpp
 *  \brief Implementation of a generic variable to file writer
 */

#ifndef IVARIABLEHDF5WRITER_HPP
#define IVARIABLEHDF5WRITER_HPP

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
#include "IoHdf5/IHdf5NWriter.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of a generic variable to file writer
    */
   class IVariableHdf5NWriter: public IoHdf5::IHdf5NWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name       Filename
          * @param ext        File extension
          * @param header     Header string of file
          * @param type       Type string of file
          * @param version    Version string of file
          * @param id         ID of the dimension space
          * @param isRegular  Is data regular?
          */
         IVariableHdf5NWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const bool isRegular);

         /**
          * @brief Destructor
          */
         virtual ~IVariableHdf5NWriter();

         /**
          * @brief Add name of expected variable to be added
          */
         void expect(const PhysicalNames::Id id);

         /**
          * @brief Get dimension space file is working on
          */
         Dimensions::Space::Id   space() const;

         /**
          * @brier Set the physical parameters of the simulation
          *
          * @param parameters Physical parameters
          */
         void setPhysical(const std::map<std::string,MHDFloat>& parameters);

         /**
          * @brier Set the mesh grid arrays
          *
          * @param mesh    Grid arrays of the mesh
          */
         void setMesh(const std::vector<Array>& mesh);

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
         virtual void write() = 0;

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
          * 
          * @param spRes      Resolution information
          */
         SharedResolution mspRes;

         /**
          * @brief Physical parameters of the simulation
          */
         std::map<std::string,MHDFloat> mPhysical;

         /**
          * @brief Storage for the mesh
          */
         std::vector<Array> mMesh;

         /**
          * @brief Set the size of the dataset
          */
         void setDatasetSize();

         /**
          * @brief Set the offsets of the dataset
          */
         void setDatasetOffsets();

         /**
          * @brief Write truncation information
          */
         void writeTruncation();

         /**
          * @brief Write Physical parameters to file
          */
         void writePhysical();

         /**
          * @brief Get iterator range to scalars
          */
         scalar_iterator_range   scalarRange();

         /**
          * @brief Get iterator range to vectors
          */
         vector_iterator_range   vectorRange();

         /**
          * @brief Is file working on regular data?
          */
         bool mIsRegular;

      private:
         /**
          * @brief Set the resolution and use it for preliminary initialisation
          */
         void setResolution(SharedResolution spRes);

         /**
          * @brief Set the maximum number of IO operations
          *
          * @param id ID of the dimension space (spectral or physical)
          */
         void setCollIo();

         /**
          * @brief The dimension space the file is working on
          */
         Dimensions::Space::Id mSpaceId;

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

   /// Typedef for a smart reference counting pointer of a Variable HDF5 numbering writer
   typedef SharedPtrMacro<IVariableHdf5NWriter>   SharedIVariableHdf5NWriter;

}
}

#endif // IVARIABLEHDF5WRITER_HPP
