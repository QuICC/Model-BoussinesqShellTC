/** 
 * @file SphereTransformSteps.cpp
 * @brief Source of the implementation of the physical <-> spectral transform steps in a whole sphere
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
#include "TransformConfigurators/SphereTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<IntegratorBranch>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<IntegratorBranch> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      if(isNL)
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, scalId, FieldType::SCALAR));
      } else
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, scalId, FieldType::SCALAR));
      }

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

   std::vector<IntegratorBranch>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 2);
      std::vector<IntegratorBranch> transform;
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
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntgPhysType::INTGDIFF, IntgPartType::INTGDIVSIN, IntgSpecType::INTGT, curlId, FieldType::VECTOR, Arithmetics::SET));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntgPhysType::INTG, IntgPartType::INTGDIFF, IntgSpecType::INTGT, curlId, FieldType::VECTOR, Arithmetics::SUB));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         // Integrate for standard fourth order spherical equation
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl Q component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntgPhysType::INTG, IntgPartType::INTGLL, IntgSpecType::INTGQ4, curlcurlId, FieldType::VECTOR, Arithmetics::SETNEG));

            // Compute curlcurl S component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntgPhysType::INTG, IntgPartType::INTGDIFF, IntgSpecType::INTGS4, curlcurlId, FieldType::VECTOR, Arithmetics::ADD));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntgPhysType::INTGDIFF, IntgPartType::INTGDIVSIN, IntgSpecType::INTGS4, curlcurlId, FieldType::VECTOR, Arithmetics::ADD));

         // Integrate for second order spherical equation
         } else if(curlcurlFlag == 1)
         {
            // Compute curlcurl Q component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntgPhysType::INTG, IntgPartType::INTGLL, IntgSpecType::INTGQ2, curlcurlId, FieldType::VECTOR, Arithmetics::SET));

            // Compute curlcurl S component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntgPhysType::INTG, IntgPartType::INTGDIFF, IntgSpecType::INTGS2, curlcurlId, FieldType::VECTOR, Arithmetics::SUB));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntgPhysType::INTGDIFF, IntgPartType::INTGDIVSIN, IntgSpecType::INTGS2, curlcurlId, FieldType::VECTOR, Arithmetics::SUB));
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
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntgPhysType::INTGDIFF, IntgPartType::INTGDIVLLDIVSIN, IntgSpecType::INTG, curlId, FieldType::VECTOR, Arithmetics::SET));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntgPhysType::INTG, IntgPartType::INTGDIVLLDIFF, IntgSpecType::INTG, curlId, FieldType::VECTOR, Arithmetics::SUB));

            // Compute Poloidal component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntgPhysType::INTG, IntgPartType::INTGDIVLL, IntgSpecType::INTGR, curlcurlId, FieldType::VECTOR, Arithmetics::SET));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }
      }

      return transform;
   }

   #else

   std::vector<IntegratorBranch>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      std::vector<IntegratorBranch> transform;

      transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::R, FieldType::VECTOR));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::THETA, FieldType::VECTOR));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::PHI, FieldType::VECTOR));

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

   std::vector<ProjectorBranch>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR));

      return transform;
   }

   std::vector<ProjectorBranch>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjSpecType::DIVR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      } 

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjSpecType::DIVR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::PHI, FieldType::GRADIENT));
      }

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::DIVR, ProjPartType::PROJLL, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjSpecType::PROJ, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::THETA, FieldType::VECTOR, Arithmetics::SET));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::DIVRDIFFR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::VECTOR, Arithmetics::ADD));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjSpecType::PROJ, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::PHI, FieldType::VECTOR, Arithmetics::SETNEG));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::DIVRDIFFR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::PHI, FieldType::VECTOR, Arithmetics::ADD));
      }

      return transform;
   } 

   std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjSpecType::DIVR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjSpecType::DIVR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::PHI, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjSpecType::DIVR, ProjPartType::PROJLL, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::CURL, Arithmetics::SET));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjSpecType::DIVRDIFFR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::SET));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::RADLAPL, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::SUB));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::DIVR2, ProjPartType::DIVSINLL, ProjPhysType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::ADD));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjSpecType::DIVRDIFFR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::PHI, FieldType::CURL, Arithmetics::SET));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::RADLAPL, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL, Arithmetics::ADD));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::DIVR2, ProjPartType::DIFFLL, ProjPhysType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL, Arithmetics::SUB));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardDivergence()
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::DIVR, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjSpecType::DIVR, ProjPartType::DIVSINDIFFSIN, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }
   
   #else

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::PHI, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjSpecType::DIVR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjSpecType::DIVR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::PHI, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjSpecType::DIVR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::R, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjSpecType::DIVR, ProjPartType::DIVSINDIFFSIN, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjSpecType::DIVR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjSpecType::DIVRDIFFR, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjSpecType::DIVR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjSpecType::DIVRDIFFR, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardDivergence()
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjSpecType::DIVR, ProjPartType::DIVSINDIFFSIN, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjSpecType::DIVR, ProjPartType::DIVSIN, ProjPhysType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

}
}
}
