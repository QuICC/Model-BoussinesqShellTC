/** 
 * @file SphereAngularMomentumWriter.hpp
 * @brief Implementation of the ASCII sphere angular momentum 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_SPHEREANGULARMOMENTUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHEREANGULARMOMENTUMWRITER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII sphere angular momentum
    */
   class SphereAngularMomentumWriter: public IVariableAsciiWriter
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
         SphereAngularMomentumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereAngularMomentumWriter();

         /**
          * @brief Compute energy for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const;

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:
         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

         /**
          * @brief Angular momentum in X, Y, Z directions
          */
         Array mMomentum;

         /**
          * @brief Angular momentum operator
          */
         MatrixI mAngMomLM;

         /**
          * @brief Angular momentum operator
          */
         Matrix mAngMomOp;

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphereAngularMomentumWriter> SharedSphereAngularMomentumWriter;

   inline bool SphereAngularMomentumWriter::isHeavy() const
   {
      return true;
   }

}
}

#endif // QUICC_IO_VARIABLE_SPHEREANGULARMOMENTUMWRITER_HPP
