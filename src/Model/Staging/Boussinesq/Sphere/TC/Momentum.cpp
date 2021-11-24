/** 
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in the Boussinesq thermal convection in a sphere model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Sphere/TC/Momentum.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/Cross.hpp"
#include "PolynomialTransforms/WorlandOperators.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace TC {

   Momentum::Momentum(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams), mConserveAngMom(false)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Momentum::~Momentum()
   {
   }

   void Momentum::setCoupling()
   {
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         int start = 1;
         bool hasMOrdering = false;
      #elif defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         int start = 0;
         bool hasMOrdering = true;
      #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      this->defineCoupling(FieldComponents::Spectral::TOR, CouplingInformation::PROGNOSTIC, start, true, false);

      this->defineCoupling(FieldComponents::Spectral::POL, CouplingInformation::PROGNOSTIC, start, true, false);

      this->mConserveAngMom = (this->bcIds().bcId(this->name()) == 1);

      if(this->mConserveAngMom)
      {
         // Initialize angular momentum operator
         this->mAngMomLM = MatrixI::Constant(2, 2, -1);

         SharedResolution spRes = this->unknown().dom(0).spRes();

         if(hasMOrdering)
         {
            // Loop over harmonic order m
            for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
               for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  int l_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

                  if(l_ == 1 && m_ == 0)
                  {
                     this->mAngMomLM(0,0) = k;
                     this->mAngMomLM(1,0) = j;
                  } else if(l_ == 1 && m_ == 1)
                  {
                     this->mAngMomLM(0,1) = k;
                     this->mAngMomLM(1,1) = j;
                  }
               }
            }
         } else
         {
            // Loop over harmonic degree l
            for(int k = 0; k < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               int l_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
               for(int j = 0; j < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
                  if(l_ == 1 && m_ == 0)
                  {
                     this->mAngMomLM(0,0) = k;
                     this->mAngMomLM(1,0) = j;
                  } else if(l_ == 1 && m_ == 1)
                  {
                     this->mAngMomLM(0,1) = k;
                     this->mAngMomLM(1,1) = j;
                  }
               }
            }
         }

         // Compute operator if required
         if(this->mAngMomLM.sum() > -4)
         {
            int nN = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            Polynomial::WorlandOperators::integralR3(this->mAngMomOp, nN, 1);
            assert(this->mAngMomOp.rows() == nN && this->mAngMomOp.cols() == 1);
            this->mAngMomOp.col(0).bottomRows(this->mAngMomOp.rows()-1) /= -this->mAngMomOp(0,0);
            this->mAngMomOp(0,0) = 0;
         }
      }
   }

   void Momentum::setNLComponents()
   {
      this->addNLComponent(FieldComponents::Spectral::TOR, 0);

      this->addNLComponent(FieldComponents::Spectral::POL, 0);
   }

   void Momentum::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      ///
      /// Compute \f$\vec u\wedge\left(\nabla\wedge\vec u\right)\f$
      ///
      switch(id)
      {
         case(FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         case(FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, this->unknown().dom(0).curl(), this->unknown().dom(0).phys(), 1.0);
            break;
         default:
            assert(false);
            break;
      }
   }

   void Momentum::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?(, need curl?)
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, true, false, true));
   }

   void Momentum::tuneSolution(const FieldComponents::Spectral::Id compId, const int k)
   {
      if(this->mConserveAngMom)
      {
         if(compId == FieldComponents::Spectral::TOR)
         {
            ArrayZ mom;
            if(this->mAngMomLM(0,1) == k)
            {
               mom = (this->mAngMomOp.transpose()*this->unknown().dom(0).total().comp(compId).profile(this->mAngMomLM(1,1),this->mAngMomLM(0,1)));
               this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(mom(0), 0, this->mAngMomLM(1,1),this->mAngMomLM(0,1));
            }

            if(this->mAngMomLM(0,0) == k)
            {
               mom = (this->mAngMomOp.transpose()*this->unknown().dom(0).total().comp(compId).profile(this->mAngMomLM(1,0),this->mAngMomLM(0,0)));
               this->rUnknown().rDom(0).rPerturbation().rComp(compId).setPoint(mom(0), 0, this->mAngMomLM(1,0),this->mAngMomLM(0,0));
            }
         }
      }
   }

}
}
}
}
}
