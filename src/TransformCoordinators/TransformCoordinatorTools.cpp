/** 
 * @file TransformCoordinatorTools.cpp
 * @brief Source of the requirement tools to work with variables and equations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "TransformCoordinators/TransformCoordinatorTools.hpp"

// Project includes
//
#include "Timers/StageTimer.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   TransformCoordinatorTools::TransformCoordinatorTools()
   {
   }

   TransformCoordinatorTools::~TransformCoordinatorTools()
   {
   }

   void TransformCoordinatorTools::init(TransformCoordinatorType& rCoord, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::vector<Transform::IntegratorTree>& integratorTree, const std::vector<Transform::ProjectorTree>& projectorTree, SharedResolution spRes, const std::map<NonDimensional::Id,MHDFloat>& runOptions)
   {
      StageTimer stage;
      stage.start("initializing transforms");

      // Get the list of required options
      std::set<NonDimensional::Id>  requests;
      rCoord.requiredOptions(requests);

      // Storage for the options' values
      std::map<NonDimensional::Id,MHDFloat> options;

      // Extract required options
      std::map<NonDimensional::Id,MHDFloat>::const_iterator it;
      for(it = runOptions.begin(); it != runOptions.end(); ++it)
      {
         if(requests.count(it->first) > 0)
         {
            options.insert(std::make_pair(it->first, it->second));
         }
      }

      // Make sure everything was present
      if(requests.size() != options.size())
      {
         throw Exception("TransformCoordinatorTools: requested options have not been set!");
      }

      // Transfer options to transforms
      rCoord.setOptions(options);

      // Initialise the transform coordinator
      rCoord.initTransforms(spRes, integratorTree, projectorTree);

      stage.done();
      stage.start("initializing communicators");

      // Initialise the communicator
      rCoord.initCommunicator(spRes);

      stage.done();
      stage.start("initializing converters");

      // Get the buffer pack sizes
      ArrayI packs1DFwd = spFwdGrouper->packs1D(integratorTree);
      ArrayI packs2DFwd = spFwdGrouper->packs2D(integratorTree);
      ArrayI packs1DBwd = spBwdGrouper->packs1D(projectorTree);
      ArrayI packs2DBwd = spBwdGrouper->packs2D(projectorTree);

      // Initialise the converters
      rCoord.communicator().initConverter(spRes, packs1DFwd, packs1DBwd, packs2DFwd, packs2DBwd, spFwdGrouper->split);

      stage.done();
   }
}
}
