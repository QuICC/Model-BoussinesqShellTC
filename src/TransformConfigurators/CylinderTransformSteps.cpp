/** 
 * @file CylinderTransformSteps.cpp
 * @brief Source of the implementation of the physical <-> spectral transform steps in a whole cylinder
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
#include "TransformConfigurators/CylinderTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

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

   #if defined GEOMHDISCC_SPATIALSCHEME_WFT_TORPOL

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
         // Integrate for curl equation, second order
         if(curlFlag == 0)
         {
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI2);
            transform.back().addEdge(Integrator2DType::INTGDIFF);
            transform.back().addEdge(Integrator1DType::INTGI4DIVR, curlId, Arithmetics::SUB);

            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI2);
            transform.back().addEdge(Integrator2DType::INTG);
            transform.back().addEdge(Integrator1DType::INTGI4DIVRDIFFR, curlId, Arithmetics::ADD);
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         // Integrate for double curl equation, fourth order
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI4D1);
            transform.back().addEdge(Integrator2DType::INTG);
            transform.back().addEdge(Integrator1DType::INTGI6DIVRDIFFR, curlcurlId, Arithmetics::ADD);

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI4D1);
            transform.back().addEdge(Integrator2DType::INTGDIFF);
            transform.back().addEdge(Integrator1DType::INTGI6DIVR, curlcurlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI4);
            transform.back().addEdge(Integrator2DType::INTG);
            transform.back().addEdge(Integrator1DType::INTGI6LAPLH, curlcurlId, Arithmetics::SUB);

         // Integrate for double curl equation, second order (typically induction equation)
         } else if(curlcurlFlag == 1)
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI4D1);
            transform.back().addEdge(Integrator2DType::INTG);
            transform.back().addEdge(Integrator1DType::INTGI6DIVRDIFFR, curlcurlId, Arithmetics::ADD);

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI4D1);
            transform.back().addEdge(Integrator2DType::INTGDIFF);
            transform.back().addEdge(Integrator1DType::INTGI6DIVR, curlcurlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGI4);
            transform.back().addEdge(Integrator2DType::INTG);
            transform.back().addEdge(Integrator1DType::INTGI6LAPLH, curlcurlId, Arithmetics::SUB);
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

      // The following assumes the physical values are obtained from a Toroidal/Poloidal decomposition
      } else
      {
         if(curlFlag == 0 && curlcurlFlag == 0)
         {
            // Compute Toroidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTGDIFF);
            transform.back().addEdge(Integrator1DType::INTGINVLAPLH, curlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTG);
            transform.back().addEdge(Integrator1DType::INTGINVLAPLH, curlId, Arithmetics::SUB);

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG);
            transform.back().addEdge(Integrator2DType::INTG);
            transform.back().addEdge(Integrator1DType::INTGINVLAPLH, curlcurlId, Arithmetics::SUB);
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
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
      transform.back().addEdge(IntegratorNDType::INTG);
      transform.back().addEdge(Integrator2DType::INTG);
      transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::R, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
      transform.back().addEdge(IntegratorNDType::INTG);
      transform.back().addEdge(Integrator2DType::INTG);
      transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::THETA, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Physical::Z, FieldType::VECTOR));
      transform.back().addEdge(IntegratorNDType::INTG);
      transform.back().addEdge(Integrator2DType::INTG);
      transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::Z, Arithmetics::ADD);

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_WFT_TORPOL

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

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::ADD);
      } 

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::Z, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req)
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      throw Exception("Second derivative is not implementated yet!");

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::ONE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::THREE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_WFT_TORPOL

   std::vector<TransformPath>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::R, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      { 
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THETA, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::LAPLH);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Z, Arithmetics::SUB);
      }

      return transform;
   } 

   std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::ADD);
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

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::R, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVRLAPLH);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::DIFF2, FieldComponents::Physical::R, Arithmetics::SUB);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         // Toroidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THETA, Arithmetics::ADD);

         // Poloidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFFLAPLH);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::ADD);

         // Poloidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF2, FieldComponents::Physical::THETA, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         // Toroidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::LAPLH);
         transform.back().addEdge(Projector2DType::PROJ);
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

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Z, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::ADD);
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

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::R, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THETA, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SUB);
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Z, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVRDIFFR);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::Z, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardDivergence()
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIVRDIFFR);
      transform.back().addEdge(Projector2DType::PROJ);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIVR);
      transform.back().addEdge(Projector2DType::DIFF);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::Z, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(Projector2DType::PROJ);
      transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_WFT_TORPOL

}
}
}
