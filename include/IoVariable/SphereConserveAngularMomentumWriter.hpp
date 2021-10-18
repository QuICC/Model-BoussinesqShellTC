/** 
 * @file SphereConserveAngularMomentumWriter.hpp
 * @brief Implementation of the ASCII sphere angular momentum conservation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_SPHERECONSERVEANGULARMOMENTUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERECONSERVEANGULARMOMENTUMWRITER_HPP

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
#include "IoVariable/SphereAngularMomentumWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII sphere angular momentum conservation
    */
   class SphereConserveAngularMomentumWriter: public SphereAngularMomentumWriter
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
         SphereConserveAngularMomentumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereConserveAngularMomentumWriter();

         /**
          * @brief Compute energy for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:
         /**
          * @brief Unit angular momentum in n = 0 mode
          */
         Array mUnitMomentum;

      private:
   };

   /// Typedef for a shared pointer of a angular momemtum conservation file
   typedef SharedPtrMacro<SphereConserveAngularMomentumWriter> SharedSphereConserveAngularMomentumWriter;

}
}

#endif // QUICC_IO_VARIABLE_SPHERECONSERVEANGULARMOMENTUMWRITER_HPP
