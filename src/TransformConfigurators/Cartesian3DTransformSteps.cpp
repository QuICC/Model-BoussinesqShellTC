/** 
 * @file Cartesian3DTransformSteps.cpp
 * @brief Source of the implementation of the 3D cartesian transform steps
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/Cartesian3DTransformSteps.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace TransformSteps {

   std::vector<TransformPath>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      if(isNL)
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::ADD);
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::ADD);
      }

      return transform;
   }

   #if defined QUICC_SPATIALSCHEME_TFF_TORPOL

   std::vector<TransformPath>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 2);
      std::vector<TransformPath> transform;
      FieldComponents::Spectral::Id curlId = components.at(0).first;
      int curlFlag = components.at(0).second;
      FieldComponents::Spectral::Id curlcurlId = components.at(1).first;
      int curlcurlFlag = components.at(1).second;

      if(isNL)
      {
         // Integrate for standard second order equation
         if(curlFlag == 0)
         {
            // Compute curl component and mean in X
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFFM);
            transform.back().addEdge(Integrator2DType::INTGM);
            transform.back().addEdge(Integrator1DType::INTGI2, curlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGDIFF);
            transform.back().addEdge(Integrator1DType::INTGI2, curlId, Arithmetics::SUB);

         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         // Integrate for standard fourth order double curl equation
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl with Dz component and mean in Y
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGDIFF);
            transform.back().addEdge(Integrator1DType::INTGI4D1, curlcurlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFFM);
            transform.back().addEdge(Integrator2DType::INTGM);
            transform.back().addEdge(Integrator1DType::INTGI4D1MI2, curlcurlId, Arithmetics::ADD);

            // Compute curlcurl without Dz component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGLAPLH);
            transform.back().addEdge(Integrator1DType::INTGI4, curlcurlId, Arithmetics::SUB);

         // Integrate for second order double curl equation
         } else if(curlcurlFlag == 1)
         {
            // Compute curlcurl with Dz component
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGDIFF);
            transform.back().addEdge(Integrator1DType::INTGI2D1, curlcurlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFFM);
            transform.back().addEdge(Integrator2DType::INTGM);
            transform.back().addEdge(Integrator1DType::INTGI2D1MI2, curlcurlId, Arithmetics::ADD);

            // Compute curlcurl without Dz component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGLAPLH);
            transform.back().addEdge(Integrator1DType::INTGI2, curlcurlId, Arithmetics::SUB);

         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

      // The following assumes the physical values are obtained from a Toroidal/Poloidal decomposition
      } else
      {
         if(curlFlag == 0 && curlcurlFlag == 0)
         {
            // Extract Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFF);
            transform.back().addEdge(Integrator2DType::INTGINVLAPLH);
            transform.back().addEdge(Integrator1DType::INTG, curlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGDIFFFINVLAPLH);
            transform.back().addEdge(Integrator1DType::INTG, curlId, Arithmetics::SUB);

            // Extract Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGINVLAPLH);
            transform.back().addEdge(Integrator1DType::INTG, curlcurlId, Arithmetics::SUB);

            // Extract X mean into Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::X, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGMEANONLY);
            transform.back().addEdge(Integrator1DType::INTG, curlId, Arithmetics::ADD);

            // Extract Y mean into Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::Y, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGMEANONLY);
            transform.back().addEdge(Integrator1DType::INTG, curlcurlId, Arithmetics::ADD);
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }
      }

      return transform;
   }

   #else

   std::vector<TransformPath>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 3);
      std::vector<TransformPath> transform;

      if(isNL)
      {
         transform.push_back(TransformPath(FieldComponents::Physical::ONE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::ONE, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::TWO, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::TWO, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::THREE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::THREE, Arithmetics::ADD);
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::ONE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::ONE, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::TWO, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::TWO, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::THREE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   #endif //defined QUICC_SPATIALSCHEME_TFF_TORPOL

   std::vector<TransformPath>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(Projector2DType::PROJ);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      return transform;
   }

   std::vector<TransformPath>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req)
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::ONE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::DIFF2);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF2);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::DIFF, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::THREE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF2, pairId, Arithmetics::ADD);
      }

      return transform;
   }

   #if defined QUICC_SPATIALSCHEME_TFF_TORPOL

   std::vector<TransformPath>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::X)->second)
      {
         // Toroidal contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::X, Arithmetics::ADD);

         // Mean contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJMEANONLY);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::X, Arithmetics::ADD);

         // Poloidal contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::X, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Y)->second)
      {
         // Toroidal contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Y, Arithmetics::SUB);

         // Poloidal contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::Y, Arithmetics::ADD);

         // Mean contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJMEANONLY);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Y, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         // Poloidal contributions
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::LAPLH);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Z, Arithmetics::SUB);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      throw Exception("Vector gradient for Toroidal/Poloidal field not implemented");
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::X)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::X, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Y)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Y, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::Z, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::X)->second)
      {
         // Toroidal contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::X, Arithmetics::ADD);

         // Poloidal contributions
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF2);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::X, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFFSLAPLH);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::X, Arithmetics::SUB);

         // Mean contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJMEANONLY);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::X, Arithmetics::SUB);
      }

      if(req.find(FieldComponents::Physical::Y)->second)
      {
         // Toroidal contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::Y, Arithmetics::ADD);

         // Mean contribution
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJMEANONLY);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Y, Arithmetics::ADD);

         // Poloidal contributions
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF2);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Y, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFFFLAPLH);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Y, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::LAPLH);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Z, Arithmetics::SUB);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardDivergence()
   {
      // The divergence is zero be construction in this case!
      throw Exception("Divergence should not be used in Toroidal/Poloidal expansion");

      std::vector<TransformPath> transform;

      return transform;
   }

   #else

   std::vector<TransformPath>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::ONE, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::TWO, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::SUB);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THREE, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardDivergence()
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIFF);
      transform.back().addEdge(Projector2DType::PROJ);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(Projector2DType::DIFF);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(Projector2DType::PROJ);
      transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      return transform;
   }

   #endif //defined QUICC_SPATIALSCHEME_TFF_TORPOL

}
}
}
