/** 
 * @file IVariableAsciiWriter.hpp
 * @brief Implementation of a generic variable to ASCII file writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IVARIABLEASCIIWRITER_HPP
#define IVARIABLEASCIIWRITER_HPP

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoAscii/IAsciiWriter.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/TransformCommSelector.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of a generic variable to ASCII file writer
    */
   class IVariableAsciiWriter: public IoAscii::IAsciiWriter
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
          * @param mode       Write mode of file
          */
         IVariableAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const IAsciiWriter::WriteMode mode = IAsciiWriter::EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~IVariableAsciiWriter();

         /**
          * @brief Add name of expected variable to be added
          *
          * @param id ID of field
          */
         void expect(const PhysicalNames::Id id);

         /**
          * @brief Get dimension space file is working on
          */
         Dimensions::Space::Id   space() const;

         /**
          * @brief Set the physical parameters of the simulation
          *
          * @param parameters Physical parameters
          */
         void setPhysical(const std::map<std::string,MHDFloat>& parameters);

         /**
          * @brief Set the simulation time parameters
          *
          * @param time       Reached simulation time
          * @param timestep   Last timestep size
          */
         void setSimTime(const MHDFloat time, const MHDFloat timestep);

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
          * @brief Perform heavy calculations
          */
         virtual void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Write State to file
          */
         virtual void write() = 0;

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const;

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
          * @brief Get resolution
          */
         const Resolution& res() const;

         /**
          * @brief Get iterator range to scalars
          */
         scalar_iterator_range   scalarRange();

         /**
          * @brief Get iterator range to vectors
          */
         vector_iterator_range   vectorRange();

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
          * @brief Time
          */
         MHDFloat mTime;

         /**
          * @brief Timestep
          */
         MHDFloat mTimestep;

      private:
         /**
          * @brief Set the resolution and use it for preliminary initialisation
          */
         void setResolution(SharedResolution spRes);

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
   typedef SharedPtrMacro<IVariableAsciiWriter>   SharedIVariableAsciiWriter;

   inline bool IVariableAsciiWriter::isHeavy() const
   {
      return false;
   }

}
}

#endif // IVARIABLEASCIIWRITER_HPP
