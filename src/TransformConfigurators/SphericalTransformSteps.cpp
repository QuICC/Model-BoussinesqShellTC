/** 
 * @file SphericalTransformSteps.cpp
 * @brief Source of the implementation of the spherical transform steps
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
#include "TransformConfigurators/SphericalTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<IntegratorBranch>  forwardScalar()
   {
      std::vector<IntegratorBranch> transform;

      transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::SCALAR, FieldType::SCALAR));

      return transform;
   }

   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLF_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

   std::vector<IntegratorBranch>  forwardVector()
   {
      std::vector<IntegratorBranch> transform;

      transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGDIFF, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::TOR, FieldType::VECTOR, Arithmetics::SET));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTGDIFF, IntegratorBranch::Intg2DType::INTGDIVSIN, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::TOR, FieldType::VECTOR, Arithmetics::SUB));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTGDIVR, FieldComponents::Spectral::POL, FieldType::VECTOR, Arithmetics::SET));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTGDIFF, IntegratorBranch::Intg1DType::INTGDIVRDIFFR, FieldComponents::Spectral::POL, FieldType::VECTOR, Arithmetics::SET));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTGDIFF, IntegratorBranch::Intg2DType::INTGDIVSIN, IntegratorBranch::Intg1DType::INTGDIVRDIFFR, FieldComponents::Spectral::POL, FieldType::VECTOR, Arithmetics::SUB));

      return transform;
   }

   #else

   std::vector<IntegratorBranch>  forwardVector()
   {
      std::vector<IntegratorBranch> transform;

      transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::R, FieldType::VECTOR));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::THETA, FieldType::VECTOR));

      transform.push_back(IntegratorBranch(FieldComponents::Physical::PHI, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::PHI, FieldType::VECTOR));

      return transform;
   }

   #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLF_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

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

   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLF_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::POL, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::VECTOR, Arithmetics::SET));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THETA, FieldType::VECTOR, Arithmetics::SUB));

      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::PHI, FieldType::VECTOR, Arithmetics::SET));

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
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::PHI)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TOR, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::DIVSIN, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::PHI, FieldType::CURL));
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

   #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_SLFM_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_BLF_TORPOL || defined GEOMHDISCC_SPATIALSCHEME_WLF_TORPOL

}
}
}
