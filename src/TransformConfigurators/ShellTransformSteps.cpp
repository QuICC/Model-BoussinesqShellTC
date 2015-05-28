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

   std::vector<IntegratorBranch>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<IntegratorBranch> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      if(isNL)
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, scalId, FieldType::SCALAR));
      } else
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, scalId, FieldType::SCALAR));
      }

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

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
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTGDIFF, IntegratorBranch::Intg2DType::INTGDIVSIN, IntegratorBranch::Intg1DType::INTGT, curlId, FieldType::VECTOR, Arithmetics::SET));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGDIFF, IntegratorBranch::Intg1DType::INTGT, curlId, FieldType::VECTOR, Arithmetics::SUB));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         // Integrate for standard fourth order spherical equation
         if(curlcurlFlag == 0)
         {
            // Compute curlcurl Q component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGLL, IntegratorBranch::Intg1DType::INTGQ4, curlcurlId, FieldType::VECTOR, Arithmetics::SETNEG));

            // Compute curlcurl S component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGDIFF, IntegratorBranch::Intg1DType::INTGS4, curlcurlId, FieldType::VECTOR, Arithmetics::ADD));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTGDIFF, IntegratorBranch::Intg2DType::INTGDIVSIN, IntegratorBranch::Intg1DType::INTGS4, curlcurlId, FieldType::VECTOR, Arithmetics::ADD));

         // Integrate for second order spherical equation
         } else if(curlcurlFlag == 1)
         {
            // Compute curlcurl Q component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGLL, IntegratorBranch::Intg1DType::INTGQ2, curlcurlId, FieldType::VECTOR, Arithmetics::SET));

            // Compute curlcurl S component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGDIFF, IntegratorBranch::Intg1DType::INTGS2, curlcurlId, FieldType::VECTOR, Arithmetics::SUB));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTGDIFF, IntegratorBranch::Intg2DType::INTGDIVSIN, IntegratorBranch::Intg1DType::INTGS2, curlcurlId, FieldType::VECTOR, Arithmetics::SUB));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

      // The following assumes the physical values are obtained froma Toroidal/Poloidal decomposition
      } else
      {
         if(curlFlag == 0 && curlcurlFlag == 0)
         {
            // Compute Toroidal component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTGDIFF, IntegratorBranch::Intg2DType::INTGDIVLLDIVSIN, IntegratorBranch::Intg1DType::INTG, curlId, FieldType::VECTOR, Arithmetics::SET));

            transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGDIVLLDIFF, IntegratorBranch::Intg1DType::INTG, curlId, FieldType::VECTOR, Arithmetics::SUB));

            // Compute Poloidal component
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGDIVLL, IntegratorBranch::Intg1DType::INTGR, curlcurlId, FieldType::VECTOR, Arithmetics::SET));
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

      transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::R, FieldType::VECTOR));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::THETA, FieldType::VECTOR));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::PHI, FieldType::VECTOR));

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

   std::vector<ProjectorBranch>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR));

      return transform;
   }

   std::vector<ProjectorBranch>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      } 

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::PHI, FieldType::GRADIENT));
      }

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::PROJLL, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THETA, FieldType::VECTOR, Arithmetics::SET));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::VECTOR, Arithmetics::ADD));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::PHI, FieldType::VECTOR, Arithmetics::SETNEG));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::PHI, FieldType::VECTOR, Arithmetics::ADD));
      }

      return transform;
   } 

   std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::PHI, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::PROJLL, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::CURL, Arithmetics::SET));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::SET));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::RADLAPL, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::SUB));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVR2, ProjectorBranch::Proj2DType::DIVSINLL, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::ADD));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::PHI, FieldType::CURL, Arithmetics::SET));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::RADLAPL, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL, Arithmetics::ADD));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVR2, ProjectorBranch::Proj2DType::DIFFLL, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL, Arithmetics::SUB));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardDivergence()
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSINDIFFSIN, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }
   
   #else

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::PHI, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::PHI, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::R, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSINDIFFSIN, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::PHI, FieldType::CURL));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardDivergence()
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSINDIFFSIN, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::PHI, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL

}
}
}
