/** 
 * @file SplittingAlgorithm.cpp
 * @brief Source of the base of the implementation of the load splitting algorithms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <algorithm>
#include <set>
#include <map>
#include <tr1/tuple>

// External includes
//

// Class include
//
#include "LoadSplitter/Algorithms/SplittingAlgorithm.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "LoadSplitter/Algorithms/SplittingTools.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   SplittingAlgorithm::SplittingAlgorithm(const int id, const int nCpu, const ArrayI& dim, const Splitting::Algorithms::Id algo)
      : mAlgo(algo), mGrouper(Splitting::Groupers::EQUATION), mId(id), mNCpu(nCpu), mDims(dim.size()), mSimDim(dim)
   {
   }

   SplittingAlgorithm::~SplittingAlgorithm()
   {
   }

   int SplittingAlgorithm::id() const
   {
      return this->mId;
   }

   int SplittingAlgorithm::nCpu() const
   {
      return this->mNCpu;
   }

   int SplittingAlgorithm::dims() const
   {
      return this->mDims;
   }

   int SplittingAlgorithm::factor(const int i) const
   {
      // Assert on index of requested factor
      assert(i < this->mFactors.size());

      return this->mFactors(i);
   }

   const ArrayI& SplittingAlgorithm::factors() const
   {
      return this->mFactors;
   }

   int SplittingAlgorithm::maxFactor() const
   {
      return this->mFactors.maxCoeff();
   }

   std::pair<int, std::pair<SharedResolution, SplittingDescription> > SplittingAlgorithm::scoreSplitting()
   {
      // Storage for all the shared core resolutions
      std::vector<SharedCoreResolution>  coreRes;

      // Storage for all the shared transform resolutions
      std::vector<SharedTransformResolution>  transformRes;

      // Loop over all CPUs
      for(int id = 0; id < this->nCpu(); id++)
      {
         // Clear content of the transform resolutions
         transformRes.clear();

         // Loop over all dimensions
         for(int j = 0; j < this->dims(); j++)
         {
            // Add new shared transform resolution
            transformRes.push_back(this->splitDimension(static_cast<Dimensions::Transform::Id>(j), id));
         }

         // Create new shared core resolution
         coreRes.push_back(SharedCoreResolution(new CoreResolution(transformRes)));
      }

      // Create shared resolution
      ArrayI transDim = this->mspScheme->getTransformSpace();
      SharedResolution  spRes(new Resolution(coreRes, this->mSimDim, transDim));

      // Add the transform setups to the resolution
      this->mspScheme->addTransformSetups(spRes);

      // Compute the score of the obtained resolution
      int score = this->computeScore(spRes);

      // Create splitting description
      SplittingDescription descr(this->mAlgo, this->mGrouper, this->mDims, this->mFactors, score);

      // Return combination of score and shared resolution/description
      return std::make_pair(score, std::make_pair(spRes,descr));
   }

   void SplittingAlgorithm::initFactors(const int nFactors)
   {
      // Initialise the storage for the factors
      this->mFactors.resize(nFactors);

      // Compute the factors
      SplittingTools::factorizeNCpu(this->mNCpuFactors, nFactors, this->nCpu());

      // Filter the factors
      SplittingTools::filterFactors(this->mNCpuFactors, nFactors, this->nCpu());
   }

   bool SplittingAlgorithm::useNextFactors()
   {
      // Get iterator through known factors
      std::list<int>::iterator  it = this->mNCpuFactors.begin();

      // Check if list is empty
      if(it == this->mNCpuFactors.end())
      {
         return false;

      // Setup system with next factors
      } else
      {
         // Extract the next factors to try and remove them from list
         for(int i = 0; i < this->mFactors.size(); i++)
         {
            this->mFactors(i) = *it;
            it = this->mNCpuFactors.erase(it);
         }

         return true;
      }
   }

   double SplittingAlgorithm::communicationScore(SharedResolution spRes, ArrayI& details)
   {
      // The worst possible value is obtained for an all-to-all communication
      // at each (possible) communication step
      int worst = (this->dims()-1)*this->nCpu();

      // Initialise current structure score
      details.resize(spRes->cpu(0)->nDim()-1);
      details.setConstant(worst);

      // Storage for the communication structure
      std::vector<std::multimap<int,int> >  structure;

      // Handle 1D resolution
      if(spRes->cpu(0)->nDim() == 1)
      {
         throw Exception("Requested computation of communication structure score for 1D resolution!");

      // Handle 2D resolution
      } else if(spRes->cpu(0)->nDim() == 2)
      {
         // Create storage for structure
         structure.push_back(std::multimap<int,int>());

         // Storage for the communication structure
         std::map<std::tr1::tuple<int,int>, int> bwdMap;

         // Storage for a tuple object
         std::tr1::tuple<int,int> point;

         // Position iterators for insert positions
         std::map<std::tr1::tuple<int,int>, int>::iterator   mapPos;

         // Loop over CPUs
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // initialise the position hint for inserts
            mapPos = bwdMap.begin();

            // Loop over second dimension
            for(int j = 0; j < spRes->cpu(cpu)->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DAT2D>(); j++)
            {
               // Loop over backward dimension
               for(int k = 0; k < spRes->cpu(cpu)->dim(Dimensions::Transform::TRA2D)->dim<Dimensions::Data::DATB1D>(j); k++)
               {
                  // Generate point information
                  point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DAT2D>(j), spRes->cpu(cpu)->dim(Dimensions::Transform::TRA2D)->idx<Dimensions::Data::DATB1D>(k,j));

                  // Get insertion position to use as next starting point to speed up insertion
                  mapPos = bwdMap.insert(mapPos, std::make_pair(point, cpu));
               }
            }
         }

         // Make sure the content is the same also
         std::set<std::pair<int,int> > filter;
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // Loop over second dimension
            for(int j = 0; j < spRes->cpu(cpu)->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(); j++)
            {
               // Loop over forward dimension
               for(int k = 0; k < spRes->cpu(cpu)->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>(j); k++)
               {
                  // Generate point information
                  point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DATF1D>(k,j), spRes->cpu(cpu)->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j));

                  // Look for same key in backward list
                  mapPos = bwdMap.find(point);

                  // Check that both position are the same
                  if(mapPos == bwdMap.end())
                  {
                     throw Exception("The computed index sets don't match!");
                  }

                  // Add corresponding communication edge to filter
                  filter.insert(std::make_pair(cpu, mapPos->second));
               }
            }
         }

         // Store obtained minimized structure
         std::set<std::pair<int,int> >::iterator filIt;
         for(filIt = filter.begin(); filIt != filter.end(); filIt++)
         {
            structure.at(static_cast<int>(Dimensions::Transform::TRA1D)).insert(*filIt);
         }

         // Clear all the data
         bwdMap.clear();

      // Handle 3D resolution
      } else if(spRes->cpu(0)->nDim() == 3)
      {
         // Extract communication structure from resolution object
         std::map<std::tr1::tuple<int,int,int>, int> bwdMap;

         // Storage for a tuple object
         std::tr1::tuple<int,int,int> point;

         // Position iterator for insert calls
         std::map<std::tr1::tuple<int,int,int>, int>::iterator   mapPos;

         // Loop over possible data exchanges
         for(int ex = 0; ex < spRes->cpu(0)->nDim()-1; ex++)
         {
            // Create storage for structure
            structure.push_back(std::multimap<int,int>());

            // Loop over CPUs
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               // initialise the position hint for inserts
               mapPos = bwdMap.begin();

               // Loop over third dimension
               for(int i = 0; i < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex+1))->dim<Dimensions::Data::DAT3D>(); i++)
               {
                  // Loop over second dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex+1))->dim<Dimensions::Data::DAT2D>(i); j++)
                  {
                     // Loop over backward dimension
                     for(int k = 0; k < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex+1))->dim<Dimensions::Data::DATB1D>(i); k++)
                     {
                        // Generate point information
                        point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex+1))->idx<Dimensions::Data::DAT2D>(j, i), spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex+1))->idx<Dimensions::Data::DAT3D>(i), spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex+1))->idx<Dimensions::Data::DATB1D>(k,i));

                        // Get insertion position to use as next starting point to speed up insertion
                        mapPos = bwdMap.insert(mapPos,std::make_pair(point, cpu));
                     }
                  }
               }
            }

            // Make sure both maps contain the same and exctract communication structure
            std::set<std::pair<int,int> > filter;
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               // Loop over third dimension
               for(int i = 0; i < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex))->dim<Dimensions::Data::DAT3D>(); i++)
               {
                  // Loop over second dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex))->dim<Dimensions::Data::DAT2D>(i); j++)
                  {
                     // Loop over forward dimension
                     for(int k = 0; k < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex))->dim<Dimensions::Data::DATF1D>(i); k++)
                     {
                        // Generate point information
                        point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex))->idx<Dimensions::Data::DATF1D>(k,i), spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex))->idx<Dimensions::Data::DAT2D>(j,i), spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(ex))->idx<Dimensions::Data::DAT3D>(i));

                        // Look for same key in backward list
                        mapPos = bwdMap.find(point);

                        // Check that both position are the same
                        if(mapPos == bwdMap.end())
                        {
                           throw Exception("The computed index sets don't match!");
                        }

                        // Add corresponding communication edge to filter
                        filter.insert(std::make_pair(cpu, mapPos->second));
                     }
                  }
               }
            }

            // Store obtained minimized structure
            std::set<std::pair<int,int> >::iterator filIt;
            for(filIt = filter.begin(); filIt != filter.end(); filIt++)
            {
               structure.at(ex).insert(*filIt);
            }

            // Clear all the data for next loop
            bwdMap.clear();
         }
      }

      // Loop over possible data exchanges
      for(int ex = 0; ex < spRes->cpu(0)->nDim()-1; ex++)
      {
         // Storage for the group sizes
         std::set<int>  groups;

         // Loop over CPUs
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            groups.insert(structure.at(ex).count(cpu));
         }

         // Store the group size (use the last one in set, to take into account case with unequal sizes. End score will be bad due to inbalance)
         details(ex) = *groups.rbegin();
      }

      // Return ratio of both structures (higher is better)
      return static_cast<double>(worst)/static_cast<double>(details.sum());
   }

   double SplittingAlgorithm::balancingScore(SharedResolution spRes, Array& balance)
   {
      // Storage for the per CPU loads for each dimension
      std::vector<std::map<int, double> >   loads;

      // Handle 1D resolution
      if(spRes->cpu(0)->nDim() == 1)
      {
         throw Exception("Requested computation of load balancing score for 1D resolution!");

      // Handle 2D resolution
      } else if(spRes->cpu(0)->nDim() == 2)
      {
         // Loop over dimensions
         for(int d = 0; d < spRes->cpu(0)->nDim(); d++)
         {
            // Create storage
            loads.push_back(std::map<int, double>());

            // Loop over CPUs
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               // Initialise CPU load to zero
               loads.at(d)[cpu] = 0.0;

               // Loop over second dimension
               for(int j = 0; j < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(d))->dim<Dimensions::Data::DAT2D>(); j++)
               {
                  // Increment load by 1
                  loads.at(d).find(cpu)->second += 1.0;
               }
            }
         }

      // Handle 3D resolution
      } else if(spRes->cpu(0)->nDim() == 3)
      {
         // Loop over dimensions
         for(int d = 0; d < spRes->cpu(0)->nDim(); d++)
         {
            // Create storage
            loads.push_back(std::map<int, double>());

            // Loop over CPUs
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               // Initialise CPU fload to zero
               loads.at(d)[cpu] = 0.0;

               // Loop over third dimension
               for(int i = 0; i < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(d))->dim<Dimensions::Data::DAT3D>(); i++)
               {
                  // Loop over second dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(static_cast<Dimensions::Transform::Id>(d))->dim<Dimensions::Data::DAT2D>(i); j++)
                  {
                     // Increment load by 1
                     loads.at(d).find(cpu)->second += 1.0;
                  }
               }
            }
         }
      }

      // Get total load
      double optimal = 0.0;
      Array perCpu(spRes->nCpu());

      // Loop over dimensions
      std::map<int, double>::const_iterator  it;
      for(int d = 0; d < spRes->cpu(0)->nDim(); d++)
      {
         // Reset loads
         optimal = 0.0;
         perCpu.setConstant(0.0);

         for(it = loads.at(d).begin(); it != loads.at(d).end(); it++)
         {
            perCpu(it->first) += it->second;
            optimal += it->second;
         }

         // Convert total load to optimal load per CPU
         optimal = optimal/spRes->nCpu();

         // Load balance
         if(perCpu.minCoeff() > optimal)
         {
            balance(d) *= optimal/perCpu.maxCoeff();
         } else if(perCpu.maxCoeff() < optimal)
         {
            balance(d) *= perCpu.minCoeff()/optimal;
         } else
         {
            balance(d) *= std::min(perCpu.minCoeff()/optimal, optimal/perCpu.maxCoeff());
         }
      }

      // Compute score
      double score = 1.0;

      for(int i = 0; i < balance.size(); i++)
      {
         score *= balance(i);
      }

      return score;
   }

}
}
