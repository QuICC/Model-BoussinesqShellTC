/** \file VisualizeResolution.cpp
 *  \brief Source of the implementation of a few formatting functions
 */

// System includes
//
#include <iostream>
#include <vector>
#include <set>
#include <tr1/tuple>

// External includes
//

// Class include
//
#include "IoTools/VisualizeResolution.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/DimensionTools.hpp"

namespace GeoMHDiSCC {

namespace IoTools {

   void VisualizeResolution::show(std::ostream& stream, std::string graphName, SharedResolution spRes)
   {
      if(spRes->cpu(0)->nDim() == 1)
      {
         VisualizeResolution::show1D(stream, graphName, spRes);

      } else if(spRes->cpu(0)->nDim() == 2)
      {
         VisualizeResolution::show2D(stream, graphName, spRes);

      } else if(spRes->cpu(0)->nDim() == 3)
      {
         VisualizeResolution::show3D(stream, graphName, spRes);
      }
   }

   void VisualizeResolution::show1D(std::ostream& stream, std::string graphName, SharedResolution spRes)
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

   void VisualizeResolution::show2D(std::ostream& stream, std::string graphName, SharedResolution spRes)
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

   void VisualizeResolution::show3D(std::ostream& stream, std::string graphName, SharedResolution spRes)
   {
      // Open graph
      stream << "graph " + graphName + " {" << std::endl;

      // Set subraph string identifier
      std::vector<std::string> cName;
      cName.push_back("A");
      cName.push_back("B");

      // Set subgraph colors
      std::vector<std::string> color;
      color.push_back("blue");
      color.push_back("green");

      // Loop over the two transposes
      for(int i = 0; i < 2; i++)
      {
         // Cast integer to dimension ID
         Dimensions::Transform::Id fwdDim = static_cast<Dimensions::Transform::Id>(i);

         // Create forward coordinates list
         std::vector<std::set<std::tr1::tuple<int,int,int> > > coordFwd;

         // Loop over all CPUs
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // Initialise coordinate list for current CPU
            coordFwd.push_back(std::set<std::tr1::tuple<int,int,int> >());

            // Loop over slow data dimension
            int i_, j_, k_;
            for(int k = 0; k < spRes->cpu(cpu)->dim(fwdDim)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               // Get global slow index
               k_ = spRes->cpu(cpu)->dim(fwdDim)->idx<Dimensions::Data::DAT3D>(k);

               // Loop over middle data dimension
               for(int j = 0; j < spRes->cpu(cpu)->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(k); ++j)
               {
                  // Get global middle index
                  j_ = spRes->cpu(cpu)->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j,k);

                  // Loop over middle Forward data dimension
                  for(int i = 0; i < spRes->cpu(cpu)->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(k); ++i)
                  {
                     // Get global fast index
                     i_ = spRes->cpu(cpu)->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,k);

                     // Add coordinate to list
                     coordFwd.back().insert(std::tr1::make_tuple(i_,k_,j_));
                  }
               }
            }
         }

         // Create backward coordinate list
         std::vector<std::set<std::tr1::tuple<int,int,int> > > coordBwd;

         // Loop over all CPUs
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // Initialise coordinate list for current CPU
            coordBwd.push_back(std::set<std::tr1::tuple<int,int,int> >());

            // Loop over slow data dimension
            int i_, j_, k_;
            for(int k = 0; k < spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               // Get global slow index
               k_ = spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT3D>(k);

               // Loop over middle data dimension
               for(int j = 0; j < spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(k); ++j)
               {
                  // Get global middle index
                  j_ = spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j,k);

                  // Loop over middle backward data dimension
                  for(int i = 0; i < spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(k); ++i)
                  {
                     // Get global fast index
                     i_ = spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,k);

                     // Add coordinate to list
                     coordBwd.back().insert(std::tr1::make_tuple(j_,i_,k_));
                  }
               }
            }
         }

         // Open subgraph
         stream << "\t" << "subgraph cluster_" + cName.at(i) + " {" << std::endl;
         // Set color of subgraph
         stream << "\t\t" << "color=" + color.at(i) + ";" << std::endl;

         // Set label of subgraph
         if(fwdDim == Dimensions::Transform::TRA1D)
         {
            stream << "\t\t" << "label=\"Transpose 1D/2D\";" << std::endl;
         } else
         {
            stream << "\t\t" << "label=\"Transpose 2D/3D\";" << std::endl;
         }

         // Loop over the CPUs
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // Create labels for CPU nodes
            stream << "\t\t" + cName.at(i) << cpu << " [label=\"CPU " << cpu << "\"];" << std::endl;
         }

         // Create set for shared coordinates
         std::set<std::tr1::tuple<int,int,int> > shared;

         // Loop over forward cpus
         for(int cpu1D = 0; cpu1D < spRes->nCpu(); cpu1D++)
         {
            // Loop over backward cpus
            for(int cpu2D = 0; cpu2D < spRes->nCpu(); cpu2D++)
            {
               // Compute coordinate intersection
               std::set_intersection(coordFwd.at(cpu1D).begin(), coordFwd.at(cpu1D).end(), coordBwd.at(cpu2D).begin(), coordBwd.at(cpu2D).end(), std::inserter(shared, shared.begin()));

               // Create link if at least one coordinate is shared
               if(shared.size() > 0)
               {
                  stream << "\t\t" + cName.at(i) << cpu1D << " -- " + cName.at(i) << cpu2D << " [color=" + color.at(i) + "];" << std::endl;
               }

               // Clear shared list
               shared.clear();
            }
         }

         // Close subgraph
         stream << "\t" << "}" << std::endl;
      }

      // Close graph
      stream << "}" << std::endl;
   }
}
}
