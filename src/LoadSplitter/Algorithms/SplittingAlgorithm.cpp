/** 
 * @file SplittingAlgorithm.cpp
 * @brief Source of the base of the implementation of the load splitting algorithms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

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

      // Initialise description
      SplittingDescription descr;

      // Load splitting might fail
      int status = 0;

      // Loop over all CPUs
      for(int id = 0; id < this->nCpu(); id++)
      {
         // Clear content of the transform resolutions
         transformRes.clear();

         // Loop over all dimensions
         for(int j = 0; j < this->dims(); j++)
         {
            SharedTransformResolution  spTRes = this->splitDimension(static_cast<Dimensions::Transform::Id>(j), id, status);

            #ifdef GEOMHDISCC_MPI
               MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            #endif //GEOMHDISCC_MPI

            // Splitting fail, abort
            if(status != 0)
            {
               break;
            }

            #ifdef GEOMHDISCC_DEBUG
               descr.vtpFiles.at(j)->representResolution(spTRes, id);
            #endif //GEOMHDISCC_DEBUG

            // Clear unused indexes for remote resolutions
            if(id != FrameworkMacro::id())
            {
               spTRes->clearIndexes();
            }

            transformRes.push_back(spTRes);
         }

         // Splitting fail, abort
         if(status != 0)
         {
            break;
         }

         // Create new shared core resolution
         coreRes.push_back(SharedCoreResolution(new CoreResolution(transformRes)));
      }

      SharedResolution  spRes;

      // Splitting was successful
      Array score = Array::Constant(4,1.0);
      if(status == 0)
      {
         // Create shared resolution
         ArrayI transDim = this->mspScheme->getTransformSpace();
         spRes = SharedResolution(new Resolution(coreRes, this->mSimDim, transDim));

         // Add the transform setups to the resolution
         this->mspScheme->addTransformSetups(spRes);

         // Add index counter to resolution
         this->mspScheme->addIndexCounter(spRes);

         // Compute the score of the obtained resolution
         score = this->computeScore(spRes);
      } else
      {
         // Set large negative score (splitting is unusable)
         score(0) = -9999;

         // Clear communication structure
         std::vector<std::multimap<int,int> >().swap(this->mCommStructure);
      }

      // Create splitting description
      descr.algorithm = this->mAlgo;
      descr.grouper = this->mGrouper;
      descr.dims = this->mDims;
      descr.factors = this->mFactors;
      descr.score = score;
      descr.structure = this->mCommStructure;

      // Return combination of score and shared resolution/description
      return std::make_pair(static_cast<int>(score.prod()), std::make_pair(spRes,descr));
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

      // Clear the communication structure
      std::vector<std::multimap<int,int> >().swap(this->mCommStructure);

      Dimensions::Transform::Id dimId;
      int i_;
      int j_;
      int k_;

      // Handle 1D resolution
      if(spRes->cpu(0)->nDim() == 1)
      {
         throw Exception("Requested computation of communication structure score for 1D resolution!");

      // Handle 2D resolution
      } else if(spRes->cpu(0)->nDim() == 2)
      {
         // Create storage for structure
         this->mCommStructure.push_back(std::multimap<int,int>());

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
            this->mCommStructure.at(static_cast<int>(Dimensions::Transform::TRA1D)).insert(*filIt);
         }

         // Clear all the data
         bwdMap.clear();

      // Handle 3D resolution
      } else if(spRes->cpu(0)->nDim() == 3)
      {
         // Simplify syntax
         typedef std::tr1::tuple<int,int,int>   Coordinate;

         // Extract communication structure from resolution object
         std::set<Coordinate> bwdMap;
         std::set<Coordinate> fwdMap;

         // Storage for a coordinate
         Coordinate point;

         // Position iterator for insert calls
         std::set<Coordinate>::iterator   mapPos;

         // Loop over possible data exchanges 
         for(int ex = 0; ex < spRes->cpu(0)->nDim()-1; ex++)
         {
            // Create storage for structure
            this->mCommStructure.push_back(std::multimap<int,int>());

            // initialise the position hint for inserts
            mapPos = bwdMap.begin();

            dimId = static_cast<Dimensions::Transform::Id>(ex+1);
            // Loop over third dimension
            for(int k = 0; k < spRes->cpu()->dim(dimId)->dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = spRes->cpu()->dim(dimId)->idx<Dimensions::Data::DAT3D>(k);

               // Loop over second dimension
               for(int j = 0; j < spRes->cpu()->dim(dimId)->dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = spRes->cpu()->dim(dimId)->idx<Dimensions::Data::DAT2D>(j,k);

                  // Loop over backward dimension
                  for(int i = 0; i < spRes->cpu()->dim(dimId)->dim<Dimensions::Data::DATB1D>(k); i++)
                  {
                     i_ = spRes->cpu()->dim(dimId)->idx<Dimensions::Data::DATB1D>(i,k);

                     // Generate point information
                     point = spRes->counter()->makeKey(dimId, i_, j_, k_);

                     // Get insertion position to use as next starting point to speed up insertion
                     mapPos = bwdMap.insert(mapPos,point);
                  }
               }
            }

            // Loop over CPUs
            MatrixI  matRemote;
            int matched = 0;
            int toMatch = -1;
            std::set<std::pair<int,int> > filter;
            dimId = static_cast<Dimensions::Transform::Id>(ex);
            for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
            {
               matched = 0;

               // Local CPU
               if(cpu == FrameworkMacro::id())
               {
                  // Loop over third dimension
                  for(int k = 0; k < spRes->cpu()->dim(dimId)->dim<Dimensions::Data::DAT3D>(); k++)
                  {
                     k_ = spRes->cpu()->dim(dimId)->idx<Dimensions::Data::DAT3D>(k);

                     // Loop over second dimension
                     for(int j = 0; j < spRes->cpu()->dim(dimId)->dim<Dimensions::Data::DAT2D>(k); j++)
                     {
                        j_ = spRes->cpu()->dim(dimId)->idx<Dimensions::Data::DAT2D>(j,k);

                        // Loop over forward dimension
                        for(int i = 0; i < spRes->cpu()->dim(dimId)->dim<Dimensions::Data::DATF1D>(k); i++)
                        {
                           i_ = spRes->cpu()->dim(dimId)->idx<Dimensions::Data::DATF1D>(i,k);

                           // Generate point information
                           point = spRes->counter()->makeKey(dimId, i_, j_, k_);

                           // Look for same key in backward list
                           mapPos = bwdMap.find(point);

                           // Key was present, drop enntry and extend filter
                           if(mapPos != bwdMap.end())
                           {
                              // Add corresponding communication edge to filter
                              filter.insert(std::make_pair(cpu, FrameworkMacro::id()));

                              // Delete found coordinate
                              bwdMap.erase(mapPos);
                           } else
                           {
                              fwdMap.insert(point);
                           }
                        }
                     }
                  }

                  // Store size of forward coordinates
                  toMatch = fwdMap.size();

               #ifdef GEOMHDISCC_MPI
                  // Convert coordinates set to matrix to send through MPI
                  matRemote.resize(3, fwdMap.size());
                  int i =0;
                  for(std::set<Coordinate>::iterator it = fwdMap.begin(); it != fwdMap.end(); ++it)
                  {
                     matRemote(0,i) = std::tr1::get<0>(*it);
                     matRemote(1,i) = std::tr1::get<1>(*it);
                     matRemote(2,i) = std::tr1::get<2>(*it);
                     i++;
                  }

                  // Broadcast size
                  MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);

                  // Broadcast data
                  MPI_Bcast(matRemote.data(), matRemote.cols()*matRemote.rows(), MPI_INT, cpu, MPI_COMM_WORLD); 

               // Remote CPU   
               } else
               {
                  // Get size
                  MPI_Bcast(&toMatch, 1, MPI_INT, cpu, MPI_COMM_WORLD);

                  // Get remote keys as matrix
                  matRemote.resize(3, toMatch);
                  MPI_Bcast(matRemote.data(), matRemote.cols()*matRemote.rows(), MPI_INT, cpu, MPI_COMM_WORLD); 

                  // Compare received data to stored indexes
                  for(int i = 0; i < toMatch; i++)
                  {
                     point = std::tr1::make_tuple(matRemote(0,i), matRemote(1,i), matRemote(2,i));

                     mapPos = bwdMap.find(point);

                     // Check if point is in backward map
                     if(mapPos != bwdMap.end())
                     {
                        // Add corresponding communication edge to filter
                        filter.insert(std::make_pair(cpu, FrameworkMacro::id()));

                        // Delete found entry
                        bwdMap.erase(mapPos);

                        // Increase matched counter
                        matched++;
                     }
                  }
               }

               // Synchronize
               FrameworkMacro::synchronize();
               #else
               }
               #endif // GEOMHDISCC_MPI
            }

            #ifdef GEOMHDISCC_MPI
               // Gather total number of match entries
               MPI_Allreduce(MPI_IN_PLACE, &matched, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            #endif // GEOMHDISCC_MPI

            // Check that everything matched
            if(toMatch != matched)
            {
               throw Exception("The computed index sets don't match!");
            }

            #ifdef GEOMHDISCC_MPI
               // Store current filter
               MatrixI locFilter(2, filter.size());
               int i = 0;
               for(std::set<std::pair<int,int> >::iterator it = filter.begin(); it != filter.end(); ++it)
               {
                  locFilter(0, i) = it->first;
                  locFilter(1, i) = it->second;
                  i++;
               }

               // Gather full communication structure
               for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
               {
                  int filterSize = 0;
                  if(cpu == FrameworkMacro::id())
                  {
                     filterSize = locFilter.cols();

                     // Get size
                     MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);

                     // Get remote keys as matrix
                     MPI_Bcast(locFilter.data(), locFilter.cols()*locFilter.rows(), MPI_INT, cpu, MPI_COMM_WORLD); 

                  } else
                  {
                     // Get size
                     MPI_Bcast(&filterSize, 1, MPI_INT, cpu, MPI_COMM_WORLD);

                     // Get remote keys as matrix
                     matRemote.resize(2, filterSize);
                     MPI_Bcast(matRemote.data(), matRemote.cols()*matRemote.rows(), MPI_INT, cpu, MPI_COMM_WORLD); 

                     for(int i = 0; i < filterSize; ++i)
                     {
                        filter.insert(std::make_pair(matRemote(0,i), matRemote(1,i)));
                     }
                  }

                  // Synchronize
                  FrameworkMacro::synchronize();
               }
            #endif // GEOMHDISCC_MPI

            // Store obtained minimized structure
            std::set<std::pair<int,int> >::iterator filIt;
            for(filIt = filter.begin(); filIt != filter.end(); filIt++)
            {
               this->mCommStructure.at(ex).insert(*filIt);
            }

            // Clear all the data for next loop
            bwdMap.clear();
            fwdMap.clear();
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
            groups.insert(this->mCommStructure.at(ex).count(cpu));
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
