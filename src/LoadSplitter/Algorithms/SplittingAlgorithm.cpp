/** \file SplittingAlgorithm.cpp
 *  \brief Source of the base of the implementation of the load splitting algorithms
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

namespace GeoMHDiSCC {

   SplittingAlgorithm::SplittingAlgorithm(const int id, const int nCpu, const ArrayI& dim, const Splittings::Algorithms::Id algo)
      :  mId(id), mNCpu(nCpu), mDims(dim.size()), mAlgo(algo), mSimDim(dim), mGrouper(Splittings::Groupers::EQUATION)
   {
   }

   SplittingAlgorithm::~SplittingAlgorithm()
   {
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
            transformRes.push_back(this->splitDimension(j, id));
         }

         // Create new shared core resolution
         coreRes.push_back(SharedCoreResolution(new CoreResolution(transformRes)));
      }

      // Create shared resolution
      SharedResolution  res(new Resolution(coreRes, this->mSimDim));

      // Compute the score of the obtained resolution
      int score = this->computeScore(res);

      // Create splitting description
      SplittingDescription descr(this->mAlgo, this->mGrouper, this->mDims, this->mFactors, score);

      // Return combination of score and shared resolution/description
      return std::make_pair(score, std::make_pair(res,descr));
   }

   void SplittingAlgorithm::factoriseNCpu(const int factors)
   {
      // Initialise the storage
      this->mFactors.resize(factors);

      // Select factorisation algorithm depending on number of factors
      if(factors == 1)
      {
         this->mFactors(0) = this->nCpu();

         // Add factor
         this->mNCpuFactors.push_back(this->nCpu());

      // Factorise CPUs into two groups
      } else if(factors == 2)
      {
         // Get the maximum factor
         int factor = static_cast<int>(std::sqrt(this->nCpu()));

         // Compute smaller factors
         while(factor > 0)
         {
            if(this->nCpu() % factor == 0)
            {
               // Add factor
               this->mNCpuFactors.push_back(factor);

               // Add nCpu / factor
               this->mNCpuFactors.push_back(this->nCpu()/factor);

               // Add reversed splitting order
               if(factor != this->nCpu()/factor)
               {
                  // Add factor
                  this->mNCpuFactors.push_back(this->nCpu()/factor);

                  // Add nCpu / factor
                  this->mNCpuFactors.push_back(factor);
               }
            }
            --factor;
         }
      } else
      {
         throw Exception("SplittingAlgorithm::factoriseNCpu", "No factorisation algorithm available for requested factors!");
      }
   }

   void SplittingAlgorithm::filterFactors()
   {
      // Get iterator through known factors
      std::list<int>::iterator  it = this->mNCpuFactors.begin();
      std::list<int>::iterator  itF;

      bool suitable;

      // Loop over all known factors
      while(it != this->mNCpuFactors.end())
      {
         // Extract factors to test
         itF = it;
         for(int i = 0; i < this->mFactors.size(); i++)
         {
            this->mFactors(i) = *itF;
            itF++;
         }

         // Test if factors are usable splitting factors
         suitable = this->confirmFactors();

         // Remove the unusable from list
         if(not suitable)
         {
            // Erase the unusable factors
            for(int i = 0; i < this->mFactors.size(); i++)
            {
               it = this->mNCpuFactors.erase(it);
            }

         // Move to the next set of factors
         } else
         {
            std::advance(it, this->mFactors.size());
         }
      }
   }

   bool SplittingAlgorithm::confirmFactors() const
   {
      // Loop over all factors
      for(int i =0; i < this->mFactors.size(); i++)
      {
         // We don't want the extrem cases (no splitting in one direction)
         if(this->factor(i) == 1 && this->nCpu() > 1)
         {
            return false;
         }
      }

      // In all other cases accept the splitting
      return true;
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

   void SplittingAlgorithm::balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id) const
   {
      // Avoid splitting with zero elements
      if(tot < parts)
      {
         throw Exception("SplittingAlgorithm::balancedSplit", "Number of parts is bigger than total!");
      }

      // Compute part assigned to id
      if(parts > 1)
      {
         nN = 0;
         n0 = 0;
         for(int i = 0; i < tot; i++)
         {
            if(i % parts == id)
            {
               nN++;
            }
            else if(i % parts < id)
            {
               n0++;
            }
         }

      // Single part, use total
      } else if(parts == 1)
      {
         n0 = 0;
         nN = tot;

      // Can't split into less than 1 part
      } else
      {
         throw Exception("SplittingAlgorithm::balancedSplit", "Number of parts < 1!");
      }
   }

   void SplittingAlgorithm::splitMapped(const std::multimap<int, int>& mapped, ArrayI &rIdx, const int id) const
   {
      // Create typedef to simply notation
      typedef  std::multimap<int, int>::const_iterator MapIt;

      // Create iterator on map
      MapIt it;

      // Create pair of iterators for range
      std::pair<MapIt, MapIt> range;

      // resize the array of indexes
      rIdx.resize(mapped.count(id));

      // Get range of indexes for given id
      range = mapped.equal_range(id);

      // Put indexes into a set to be sure to get them sorted
      std::set<int>  sorter;
      for(it = range.first; it != range.second; ++it)
      {
         sorter.insert(it->second);
      }

      // Extract the ordered indexes from set and store in output array
      std::set<int>::iterator setIt;
      int i = 0;
      for(setIt = sorter.begin(); setIt != sorter.end(); ++setIt, ++i)
      {
         rIdx(i) = *setIt;
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
         throw Exception("SplittingAlgorithm::communicationScore", "Tried computation of communication structure score for 1D resolution!");

      // Handle 2D resolution
      } else if(spRes->cpu(0)->nDim() == 2)
      {
         // Create storage for structure
         structure.push_back(std::multimap<int,int>());

         // Storage for the communication structure
         std::map<std::tr1::tuple<int,int>, int> fwdMap, bwdMap;

         // Storage for a tuple object
         std::tr1::tuple<int,int> point;

         // Position iterators for insert positions
         std::map<std::tr1::tuple<int,int>, int>::iterator   mapPos;

         // Loop over CPUs
         int ex = 0;
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // Initialise the position hint for inserts
            mapPos = fwdMap.begin();

            // Loop over second dimension
            for(int j = 0; j < spRes->cpu(cpu)->dim(ex)->dim2D(); j++)
            {
               // Loop over forward dimension
               for(int k = 0; k < spRes->cpu(cpu)->dim(ex)->dimFwd(j); k++)
               {
                  // Generate point information
                  point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(ex)->idxFwd(k,j), spRes->cpu(cpu)->dim(ex)->idx2D(j));

                  // Get insertion position to use as next starting point to speed up insertion
                  mapPos = fwdMap.insert(mapPos, std::make_pair(point, cpu));
               }
            }

            // initialise the position hint for inserts
            mapPos = bwdMap.begin();

            // Loop over second dimension
            for(int j = 0; j < spRes->cpu(cpu)->dim(ex+1)->dim2D(); j++)
            {
               // Loop over backward dimension
               for(int k = 0; k < spRes->cpu(cpu)->dim(ex+1)->dimBwd(j); k++)
               {
                  // Generate point information
                  point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(ex+1)->idx2D(j), spRes->cpu(cpu)->dim(ex+1)->idxBwd(k,j));

                  // Get insertion position to use as next starting point to speed up insertion
                  mapPos = bwdMap.insert(mapPos, std::make_pair(point, cpu));
               }
            }
         }

         // Check that both sets have the same size
         if(fwdMap.size() != bwdMap.size())
         {
            throw Exception("SplittingAlgorithm::communicationScore", "The size of the computed index sets don't match!");
         }

         // Make sure the content is the same also
         std::map<std::tr1::tuple<int,int>, int>::const_iterator it;
         std::set<std::pair<int,int> > filter;
         mapPos = bwdMap.begin();
         for(it = fwdMap.begin(); it != fwdMap.end(); it++)
         {
            // Check that both position are the same
            if(it->first != mapPos->first)
            {
               throw Exception("SplittingAlgorithm::communicationScore", "The computed index sets don't match!");
            }

            // Add corresponding communication edge to filter
            filter.insert(std::make_pair(it->second, mapPos->second));

            // Increment second iterator
            mapPos++;
         }

         // Store obtained minimized structure
         std::set<std::pair<int,int> >::iterator filIt;
         for(filIt = filter.begin(); filIt != filter.end(); filIt++)
         {
            structure.at(ex).insert(*filIt);
         }

         // Clear all the data
         fwdMap.clear();
         bwdMap.clear();

      // Handle 3D resolution
      } else if(spRes->cpu(0)->nDim() == 3)
      {
         // Extract communication structure from resolution object
         std::map<std::tr1::tuple<int,int,int>, int> fwdMap, bwdMap;

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
               mapPos = fwdMap.begin();

               // Loop over third dimension
               for(int i = 0; i < spRes->cpu(cpu)->dim(ex)->dim3D(); i++)
               {
                  // Loop over second dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(ex)->dim2D(i); j++)
                  {
                     // Loop over forward dimension
                     for(int k = 0; k < spRes->cpu(cpu)->dim(ex)->dimFwd(j,i); k++)
                     {
                        // Generate point information
                        point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(ex)->idxFwd(k,j,i), spRes->cpu(cpu)->dim(ex)->idx2D(j,i), spRes->cpu(cpu)->dim(ex)->idx3D(i));

                        // Get insertion position to use as next starting point to speed up insertion
                        mapPos = fwdMap.insert(mapPos, std::make_pair(point, cpu));
                     }
                  }
               }

               // initialise the position hint for inserts
               mapPos = bwdMap.begin();

               // Loop over third dimension
               for(int i = 0; i < spRes->cpu(cpu)->dim(ex+1)->dim3D(); i++)
               {
                  // Loop over second dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(ex+1)->dim2D(i); j++)
                  {
                     // Loop over backward dimension
                     for(int k = 0; k < spRes->cpu(cpu)->dim(ex+1)->dimBwd(j,i); k++)
                     {
                        // Generate point information
                        point = std::tr1::make_tuple(spRes->cpu(cpu)->dim(ex+1)->idx2D(j, i), spRes->cpu(cpu)->dim(ex+1)->idx3D(i), spRes->cpu(cpu)->dim(ex+1)->idxBwd(k,j,i));

                        // Get insertion position to use as next starting point to speed up insertion
                        mapPos = bwdMap.insert(mapPos,std::make_pair(point, cpu));
                     }
                  }
               }
            }

            // Both sets should have the same set. Padding and unaliased runs unfortunately invalidate this check
            // As a workaround let's only check for correct order. In addition with the next check we should be on 
            // the safe side
            if(fwdMap.size() > bwdMap.size())
            {
               throw Exception("SplittingAlgorithm::communicationScore", "The size of the computed index sets don't match!");
            }

            // Make sure both maps contain the same and exctract communication structure
            std::map<std::tr1::tuple<int,int,int>, int>::const_iterator it;
            std::set<std::pair<int,int> > filter;
            for(it = fwdMap.begin(); it != fwdMap.end(); it++)
            {
               // Look for same key in backward list
               mapPos = bwdMap.find(it->first);

               // Check that both position are the same
               if(mapPos == bwdMap.end())
               {
                  throw Exception("SplittingAlgorithm::communicationScore", "The computed index sets don't match!");
               }

               // Add corresponding communication edge to filter
               filter.insert(std::make_pair(it->second, mapPos->second));
            }

            // Store obtained minimized structure
            std::set<std::pair<int,int> >::iterator filIt;
            for(filIt = filter.begin(); filIt != filter.end(); filIt++)
            {
               structure.at(ex).insert(*filIt);
            }

            // Clear all the data for next loop
            fwdMap.clear();
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
         throw Exception("SplittingAlgorithm::balancingScore", "Tried computation of load balancing score for 1D resolution!");

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
               for(int j = 0; j < spRes->cpu(cpu)->dim(d)->dim2D(); j++)
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
               for(int i = 0; i < spRes->cpu(cpu)->dim(d)->dim3D(); i++)
               {
                  // Loop over second dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(d)->dim2D(i); j++)
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
