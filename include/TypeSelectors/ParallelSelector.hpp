/** 
 * @file ParallelSelector.hpp
 * @brief Typedefs for the parallelisation algorithms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PARALLELSELECTOR_HPP
#define PARALLELSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "LoadSplitter/Algorithms/SplittingDescription.hpp"
#include "TransformConfigurators/ForwardSerialConfigurator.hpp"
#include "TransformConfigurators/BackwardSerialConfigurator.hpp"
#include "TransformConfigurators/ForwardSingle1DConfigurator.hpp"
#include "TransformConfigurators/BackwardSingle1DConfigurator.hpp"
#include "TransformConfigurators/ForwardSingle2DConfigurator.hpp"
#include "TransformConfigurators/BackwardSingle2DConfigurator.hpp"
#include "TransformConfigurators/ForwardTubularConfigurator.hpp"
#include "TransformConfigurators/BackwardTubularConfigurator.hpp"
#include "TransformGroupers/ForwardEquationGrouper.hpp"
#include "TransformGroupers/BackwardEquationGrouper.hpp"
#include "TransformGroupers/ForwardSingle1DGrouper.hpp"
#include "TransformGroupers/BackwardSingle1DGrouper.hpp"
#include "TransformGroupers/ForwardSingle2DGrouper.hpp"
#include "TransformGroupers/BackwardSingle2DGrouper.hpp"
#include "TransformGroupers/ForwardTransformGrouper.hpp"
#include "TransformGroupers/BackwardTransformGrouper.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      /// Transform configurator selector template
      template <Splitting::Algorithms::Id TAlgo> struct ConfigSelector;

      /// Transform configurator selector specialised for SERIAL case
      template <> struct ConfigSelector<Splitting::Algorithms::SERIAL>
      {
         /// Typedef for forward configurator
         typedef ForwardSerialConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSerialConfigurator  BwdConfigType;
      };

      #ifdef GEOMHDISCC_MPIALGO_SINGLE1D
      /// Transform configurator selector specialised for SINGLE1D case
      template <> struct ConfigSelector<Splitting::Algorithms::SINGLE1D>
      {
         /// Typedef for forward configurator
         typedef ForwardSingle1DConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSingle1DConfigurator  BwdConfigType;
      };
      #endif //GEOMHDISCC_MPIALGO_SINGLE1D

      #ifdef GEOMHDISCC_MPIALGO_SINGLE2D
      /// Transform configurator selector specialised for SINGLE2D case
      template <> struct ConfigSelector<Splitting::Algorithms::SINGLE2D>
      {
         /// Typedef for forward configurator
         typedef ForwardSingle2DConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSingle2DConfigurator  BwdConfigType;
      };
      #endif //GEOMHDISCC_MPIALGO_SINGLE2D

      #ifdef GEOMHDISCC_MPIALGO_TUBULAR
      /// Transform configurator selector specialised for TUBULAR case
      template <> struct ConfigSelector<Splitting::Algorithms::TUBULAR>
      {
         /// Typedef for forward configurator
         typedef ForwardTubularConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardTubularConfigurator  BwdConfigType;
      };
      #endif //GEOMHDISCC_MPIALGO_TUBULAR

      #ifdef GEOMHDISCC_MPIALGO_COUPLED2D
      /// Transform configurator selector specialised for COUPLED2D case
      template <> struct ConfigSelector<Splitting::Algorithms::COUPLED>
      {
         /// Typedef for forward configurator
         typedef ForwardSingle1DConfigurator   FwdConfigType;

         /// Typedef for forward configurator
         typedef BackwardSingle1DConfigurator  BwdConfigType;
      };
      #endif //GEOMHDISCC_MPIALGO_COUPLED2D

      /// Transform grouper selector template
      template <Splitting::Groupers::Id TGrouper,Splitting::Algorithms::Id TAlgo> struct GrouperSelector;

      #ifdef GEOMHDISCC_TRANSGROUPER_EQUATION
      /// Transform grouper selector for EQUATION grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::EQUATION,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardEquationGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardEquationGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };
      #endif //GEOMHDISCC_TRANSGROUPER_EQUATION

      #ifdef GEOMHDISCC_TRANSGROUPER_SINGLE1D
      /// Transform grouper selector for SINGLE1D grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::SINGLE1D,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardSingle1DGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardSingle1DGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };
      #endif //GEOMHDISCC_TRANSGROUPER_SINGLE1D

      #ifdef GEOMHDISCC_TRANSGROUPER_SINGLE2D
      /// Transform grouper selector for SINGLE2D grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::SINGLE2D,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardSingle2DGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardSingle2DGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };
      #endif //GEOMHDISCC_TRANSGROUPER_SINGLE2D

      #ifdef GEOMHDISCC_TRANSGROUPER_TRANSFORM
      /// Transform grouper selector for TRANSFORM grouper
      template <Splitting::Algorithms::Id TAlgo> struct GrouperSelector<Splitting::Groupers::TRANSFORM,TAlgo>: public ConfigSelector<TAlgo>
      {
         /// Typedef for forward grouper
         typedef ForwardTransformGrouper<typename ConfigSelector<TAlgo>::FwdConfigType>   FwdGrouperType;

         /// Typedef for backward grouper
         typedef BackwardTransformGrouper<typename ConfigSelector<TAlgo>::BwdConfigType>  BwdGrouperType;
      };
      #endif //GEOMHDISCC_TRANSGROUPER_TRANSFORM
   }

   namespace Parallel
   {
      /**
       * @brief Set the transform grouper depending on setup
       */
      void setGrouper(const SplittingDescription& descr, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper);

      /**
       * @brief Templated grouper selector 
       */
      template <Splitting::Groupers::Id TGroup, Splitting::Algorithms::Id TAlgo> void setGrouper(Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper);

      /**
       * @brief Templated grouper selector 
       */
      template <Splitting::Groupers::Id TGroup> void setGrouper(const Splitting::Algorithms::Id algo, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper);

      template <Splitting::Groupers::Id TGroup, Splitting::Algorithms::Id TAlgo> void setGrouper(Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper)
      {
         SharedPtrMacro<typename Transform::GrouperSelector<TGroup,TAlgo>::FwdGrouperType> spFwd(new typename Transform::GrouperSelector<TGroup,TAlgo>::FwdGrouperType);
         SharedPtrMacro<typename Transform::GrouperSelector<TGroup,TAlgo>::BwdGrouperType> spBwd(new typename Transform::GrouperSelector<TGroup,TAlgo>::BwdGrouperType);

         spFwdGrouper = spFwd;
         spBwdGrouper = spBwd;
      }
   
      template <Splitting::Groupers::Id TGroup> void setGrouper(const Splitting::Algorithms::Id algo, Transform::SharedIForwardGrouper& spFwdGrouper, Transform::SharedIBackwardGrouper& spBwdGrouper)
      {
         if(algo == Splitting::Algorithms::SERIAL)
         {
            setGrouper<TGroup,Splitting::Algorithms::SERIAL>(spFwdGrouper, spBwdGrouper);

      #ifdef GEOMHDISCC_MPI
         #ifdef GEOMHDISCC_MPIALGO_SINGLE1D
         } else if(algo == Splitting::Algorithms::SINGLE1D)
         {
            setGrouper<TGroup,Splitting::Algorithms::SINGLE1D>(spFwdGrouper, spBwdGrouper);
         #endif //GEOMHDISCC_MPIALGO_SINGLE1D
         #ifdef GEOMHDISCC_MPIALGO_SINGLE2D
         } else if(algo == Splitting::Algorithms::SINGLE2D)
         {
            setGrouper<TGroup,Splitting::Algorithms::SINGLE2D>(spFwdGrouper, spBwdGrouper);
         #endif //GEOMHDISCC_MPIALGO_SINGLE2D
         #ifdef GEOMHDISCC_MPIALGO_TUBULAR
         } else if(algo == Splitting::Algorithms::TUBULAR)
         {
            setGrouper<TGroup,Splitting::Algorithms::TUBULAR>(spFwdGrouper, spBwdGrouper);
         #endif //GEOMHDISCC_MPIALGO_TUBULAR
         #ifdef GEOMHDISCC_MPIALGO_COUPLED2D
         } else if(algo == Splitting::Algorithms::COUPLED2D)
         {
            setGrouper<TGroup,Splitting::Algorithms::COUPLED2D>(spFwdGrouper, spBwdGrouper);
         #endif //GEOMHDISCC_MPIALGO_COUPLED2D
      #endif //GEOMHDISCC_MPI
         } else
         {
            throw Exception("Unknown algorithm for transform grouper setup");
         }
      }
   }

}

#endif // PARALLELSELECTOR_HPP
