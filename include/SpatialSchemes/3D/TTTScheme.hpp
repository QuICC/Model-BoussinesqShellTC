/** 
 * @file TTTScheme.hpp
 * @brief Implementation of the Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TTTSCHEME_HPP
#define TTTSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/FftSelector.hpp"
#include "Enums/Splitting.hpp"
#include "Resolutions/Resolution.hpp"
#include "SpatialSchemes/3D/IRegular3DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace QuICC {

namespace Schemes {

   /**
    * @brief Implementation of Chebyshev(FFT) + Chebyshev(FFT) + Chebyshev(FFT) scheme
    */
   class TTTScheme: public IRegular3DScheme
   {
      public:
         /**
          * @brief Get type string for the scheme
          */
         static std::string type();

         /**
          * @brief Tune the shared resolution used by simulation
          */
         static void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr);

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncations 
          */
         explicit TTTScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~TTTScheme(); 

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const;

         /**
          * @brief Add the transform setups to resolution
          */
         virtual void addTransformSetups(SharedResolution spRes) const;
         
      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions();

         /**
          * @brief Set transform costs
          */
         virtual void setCosts();

         /**
          * @brief Set transform scalings
          */
         virtual void setScalings();

         /**
          * @brief Set transform memory footprint
          */
         virtual void setMemoryScore();

      private:
         /**
          * @brief Construct setup object for first transform
          */
         Transform::SharedFftSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::SharedFftSetup  spSetup2D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for third transform
          */
         Transform::SharedFftSetup  spSetup3D(SharedResolution spRes) const; 
   };

}
}

#endif // TTTSCHEME_HPP
