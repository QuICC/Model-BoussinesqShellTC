/** 
 * @file VisualizationFileWriter.hpp
 * @brief Implementation of the HDF5 visualisation file writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef VISUALIZATIONFILEWRITER_HPP
#define VISUALIZATIONFILEWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableHdf5NWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * \brief Implementation of the HDF5 visualisation file writer
    */
   class VisualizationFileWriter: public IVariableHdf5NWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type    Type of the file (typically scheme name)
          */
         VisualizationFileWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~VisualizationFileWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:
         /**
          * @brief Write the mesh to file
          */
         void writeMesh();

         /**
          * @brief Create group for scalar field
          *
          * @param name    Name of the field
          * @param scalar  Scalar field values
          */
         void writePhysicalScalar(const std::string& name, const Datatypes::PhysicalScalarType& scalar);

         /**
          * @brief Create group for vector field
          *
          * @param name    Name of the field
          * @param vector  Vector of components
          */
         void writePhysicalVector(const std::string& name, const std::vector<Datatypes::PhysicalScalarType>& vector);

      private:

   };

   /// Typedef for a shared pointer of a Hdf5Writer
   typedef SharedPtrMacro<VisualizationFileWriter> SharedVisualizationFileWriter;

}
}

#endif // VISUALIZATIONFILEWRITER_HPP
