/** \file Equation2DEigenTools.cpp
 *  \brief Source of the tools for schemes with two eigen direction
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Tools/Equation2DEigenTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   void Equation2DEigenTools::boundaryBlock2DEigen(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D)
   {
      // Get 1D dimension
      int n1D = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create equation ID
      SpectralFieldId eqId = std::make_pair(eq.name(), compId);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(n1D);

      // Create spectral boundary operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::BcType bound1D(n1D);

      // Initialise output matrices
      mat.first.resize(n1D,n1D);
      mat.second.resize(n1D,n1D);

      if(c1D != 0)
      {
         // Safety check to abort early
         if(!eq.bcIds().hasEquation(eqId))
         {
            throw Exception("Missing condition(s) for boundary operator!");
         }

         // Set boundary conditions on fieldId
         if(eq.bcIds().hasField(eqId,fieldId))
         {
            // Set 1D boundary conditions
            if(eq.bcIds().bcs(eqId,fieldId).count(Dimensions::Simulation::SIM1D) > 0)
            {
               mat = Spectral::BoundaryConditions::tauMatrix(bound1D, eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second);
               if(mat.first.nonZeros() > 0)
               {
                  if(c1D != 1.0)
                  {
                     mat.first *= c1D;
                  }
               }

               if(mat.second.nonZeros() > 0)
               {
                  if(c1D != 1.0)
                  {
                     mat.second *= c1D;
                  }
               }
            }
         }

         // Prune matrices for safety
         mat.first.prune(1e-32);
         mat.second.prune(1e-32);
      }
   }
}
}
