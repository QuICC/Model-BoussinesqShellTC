/** 
 * @file ShellTransformSteps.cpp
 * @brief Source of the implementation of the physical <-> spectral transform steps in a spherical shell
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
#include "TransformConfigurators/ShellTransformSteps.hpp"

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
         transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
         transform.back().addEdge(Integrator2DType::INTG, -1, Arithmetics::SET);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::SET);
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
         transform.back().addEdge(Integrator2DType::INTG, -1, Arithmetics::SET);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::SET);
      }

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

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
            // Compute curl component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIVSIN, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGT, curlId, Arithmetics::SET);

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGT, curlId, Arithmetics::SUB);
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         // Integrate for standard fourth order spherical equation
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGLL, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGQ4, curlcurlId, Arithmetics::SETNEG);

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGS4, curlcurlId, Arithmetics::ADD);

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIVSIN, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGS4, curlcurlId, Arithmetics::ADD);

         // Integrate for second order spherical equation
         } else if(curlcurlFlag == 1)
         {
            // Compute curlcurl Q component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGLL, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGQ2, curlcurlId, Arithmetics::SET);

            // Compute curlcurl S component
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGS2, curlcurlId, Arithmetics::SUB);

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIVSIN, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGS2, curlcurlId, Arithmetics::SUB);
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
            transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTGDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIVLLDIVSIN, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTG, curlId, Arithmetics::SET);

            transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIVLLDIFF, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTG, curlId, Arithmetics::SUB);

            // Compute Poloidal component
            transform.push_back(TransformPath(FieldComponents::Physical::R, FieldType::VECTOR));
            transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator2DType::INTGDIVLL, -1, Arithmetics::SET);
            transform.back().addEdge(Integrator1DType::INTGR, curlcurlId, Arithmetics::SET);
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
      transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
      transform.back().addEdge(Integrator2DType::INTG, -1, Arithmetics::SET);
      transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::R, Arithmetics::SET);

      transform.push_back(TransformPath(FieldComponents::Physical::THETA, FieldType::VECTOR));
      transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
      transform.back().addEdge(Integrator2DType::INTGDIVLL, -1, Arithmetics::SET);
      transform.back().addEdge(Integrator1DType::INTGR, FieldComponents::Spectral::THETA, Arithmetics::SET);

      transform.push_back(TransformPath(FieldComponents::Physical::PHI, FieldType::VECTOR));
      transform.back().addEdge(IntegratorNDType::INTG, -1, Arithmetics::SET);
      transform.back().addEdge(Integrator2DType::INTG, -1, Arithmetics::SET);
      transform.back().addEdge(Integrator1DType::INTGR  FieldComponents::Spectral::PHI, Arithmetics::SET);

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

   std::vector<TransformPath>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
      transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::SET);

      return transform;
   }

   std::vector<TransformPath>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SET);
      } 

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::PHI, Arithmetics::SET);
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
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::SET);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::SET);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::SET);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::SET);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::SET);
      }

      pairId = std::make_pair(FieldComponents::Physical::THREE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::SET);
      }

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

   std::vector<TransformPath>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJLL, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      { 
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THETA, Arithmetics::SET);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIVRDIFFR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::PHI, Arithmetics::SETNEG);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::DIVRDIFFR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::PHI, Arithmetics::ADD);
      }

      return transform;
   } 

   std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::PHI, Arithmetics::SET);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJLL, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         // Toroidal part
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVRDIFFR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SET);

         // Poloidal part 1
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::RADLAPL, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THETA, Arithmetics::SUB);

         // Poloidal part 2
         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR2, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSINLL, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THETA, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TOR, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVRDIFFR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::PHI, Arithmetics::SET);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::RADLAPL, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::PHI, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Spectral::POL, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR2, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFFLL, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::PHI, Arithmetics::SUB);
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
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::PHI, Arithmetics::SET);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SET);
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::PHI, Arithmetics::SET);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::R, Arithmetics::SETNEG);

         transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSINDIFFSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::R, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THETA, Arithmetics::SET);

         transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVRDIFFR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THETA, Arithmetics::SUB);
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::DIFF, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::PHI, Arithmetics::SETNEG);

         transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIVRDIFFR, -1, Arithmetics::SET);
         transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::PHI, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardDivergence()
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::R, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIFF, -1, Arithmetics::SET);
      transform.back().addEdge(Projector2DType::PROJ, -1, Arithmetics::SET);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::SET);

      transform.push_back(TransformPath(FieldComponents::Spectral::THETA, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
      transform.back().addEdge(Projector2DType::DIVSINDIFFSIN, -1, Arithmetics::SET);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::PHI, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIVR, -1, Arithmetics::SET);
      transform.back().addEdge(Projector2DType::DIVSIN, -1, Arithmetics::SET);
      transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

}
}
}
