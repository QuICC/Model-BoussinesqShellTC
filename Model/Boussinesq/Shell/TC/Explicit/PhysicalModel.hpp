/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq thermal convection in a spherical
 * shell (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_TC_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_TC_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Model/Boussinesq/Shell/TC/ITCModel.hpp"
#include "QuICC/SpatialScheme/3D/SLFl.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace TC {

namespace Explicit {

/**
 * @brief Implementation of the Boussinesq thermal convection spherical shell
 * model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */
class PhysicalModel : public ITCModel
{
public:
   /// Typedef for the spatial scheme used
   typedef SpatialScheme::SLFl SchemeType;

   /**
    * @brief Constructor
    */
   PhysicalModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~PhysicalModel() = default;

   /// Python model script module name
   virtual std::string PYMODULE() override;

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

} // namespace Explicit
} // namespace TC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_TC_EXPLICIT_PHYSICALMODEL_HPP
