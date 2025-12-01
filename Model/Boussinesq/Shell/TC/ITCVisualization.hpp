/**
 * @file ITCVisualization.hpp
 * @brief Implementation of the Boussinesq thermal convection in a spherical
 * shell (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SHELL_TC_ITCVISUALIZATION_HPP
#define QUICC_MODEL_BOUSSINESQ_SHELL_TC_ITCVISUALIZATION_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Generator/VisualizationGenerator.hpp"
#include "QuICC/Model/IVisualizationGeneratorBuilder.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace TC {

/**
 * @brief Implementation of the Boussinesq thermal convection spherical shell
 * model (Toroidal/Poloidal formulation)
 */
class ITCVisualization : public IVisualizationGeneratorBuilder<VisualizationGenerator>
{
public:
   /**
    * @brief Constructor
    */
   ITCVisualization() = default;

   /**
    * @brief Destructor
    */
   virtual ~ITCVisualization() = default;

   /// Formulation used for vector fields
   virtual VectorFormulation::Id SchemeFormulation() override;

   /**
    * @brief Version string
    */
   std::string version() const final;

   /**
    * @brief Add the visualization generation equations
    *
    * @param spGen   Shared visualization generator
    */
   virtual void addVisualizers(SharedVisualizationGenerator spVis) override;

protected:
private:
};

} // namespace TC
} // namespace Shell
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SHELL_TC_ITCVISUALIZATION_HPP
