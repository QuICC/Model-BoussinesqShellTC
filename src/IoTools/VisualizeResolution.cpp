/** \file VisualizeResolution.cpp
 *  \brief Source of the implementation of a few formatting functions
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "IoTools/VisualizeResolution.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"

namespace GeoMHDiSCC {

namespace IoTools {

   void VisualizeResolution::show(std::ostream& stream, SharedResolution spRes)
   {
      if(spRes->cpu(0)->nDim() == 1)
      {
         VisualizeResolution::show1D(stream, spRes);

      } else if(spRes->cpu(0)->nDim() == 2)
      {
         VisualizeResolution::show2D(stream, spRes);

      } else if(spRes->cpu(0)->nDim() == 3)
      {
         VisualizeResolution::show3D(stream, spRes);
      }
   }

   void VisualizeResolution::show1D(std::ostream& stream, SharedResolution spRes)
   {
      for(int j = 0; j < spRes->cpu(0)->nDim(); j++)
      {
         // Loop over CPUs
         for(int i = 0; i < spRes->nCpu(); i++)
         {
            for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATF1D>(); m++)
            {
               stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATF1D>(m) << std::endl;
            }
         }
         stream << "- -" << std::endl;

         for(int i = 0; i < spRes->nCpu(); i++)
         {
            for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATB1D>(); m++)
            {
               stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATB1D>(m) << std::endl;
            }
         }
         stream << "- -" << std::endl;
      }
   }

   void VisualizeResolution::show2D(std::ostream& stream, SharedResolution spRes)
   {
      for(int j = 0; j < spRes->cpu(0)->nDim(); j++)
      {
         // Loop over CPUs
         for(int i = 0; i < spRes->nCpu(); i++)
         {
            for(int l = 0; l < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT2D>(); l++)
            {
               for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATF1D>(l); m++)
               {
                  stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATF1D>(m,l) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT2D>(l) << std::endl;
               }
            }
         }
         stream << "- - -" << std::endl;

         for(int i = 0; i < spRes->nCpu(); i++)
         {
            for(int l = 0; l < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT2D>(); l++)
            {
               for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATB1D>(l); m++)
               {
                  stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATB1D>(m,l) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT2D>(l) << std::endl;
               }
            }
         }
         stream << "- - -" << std::endl;
      }
   }

   void VisualizeResolution::show3D(std::ostream& stream, SharedResolution spRes)
   {
      int bwdDim = 0;

      for(int j = 0; j < spRes->cpu(0)->nDim(); j++)
      {
         // Loop over CPUs
         for(int i = 0; i < spRes->nCpu(); i++)
         {
            for(int k = 0; k < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT3D>(); k++)
            {
               for(int l = 0; l < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT2D>(k); l++)
               {
                  for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATF1D>(k); m++)
                  {
                     stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATF1D>(m,k) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT2D>(l,k) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT3D>(k) << std::endl;
                  }
               }
            }
         }
         stream << "- - - -" << std::endl;

         for(int i = 0; i < spRes->nCpu(); i++)
         {
            for(int k = 0; k < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT3D>(); k++)
            {
               for(int l = 0; l < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT2D>(k); l++)
               {
                  if(static_cast<Dimensions::Transform::Id>(j) != Dimensions::Transform::TRA3D)
                  {
                     bwdDim = spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATB1D>(k);
                  } else
                  {
                     bwdDim = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
                  }
                  for(int m = 0; m < bwdDim; m++)
                  {
                     stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATB1D>(m,k) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT2D>(l,k) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT3D>(k) << std::endl;
                  }
               }
            }
         }
         stream << "- - - -" << std::endl;
      }
   }
}
}
