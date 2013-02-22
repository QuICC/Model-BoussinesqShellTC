/** \file LoadSplitter.cpp
 *  \brief Source of the workload splitter
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "Base/LoadSplitter/LoadSplitter.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/Enums/Splittings.hpp"
#include "Base/IO/Ascii/FormatToolbox.hpp"
// Splitting algorithms
#include "LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "LoadSplitter/Algorithms/SerialSplitting.hpp"
#include "LoadSplitter/Algorithms/SingleSplitting.hpp"
#include "LoadSplitter/Algorithms/TubularSplitting.hpp"

namespace GeoMHDiSCC {

   LoadSplitter::LoadSplitter(const int id, const int nCpu)
      : mSkip(0), mId(id), mNCpu(nCpu)
   {
   }

   int LoadSplitter::nCpu() const
   {
      return this->mNCpu;
   }

   int LoadSplitter::id() const
   {
      return this->mId;
   }

   bool LoadSplitter::moreSplittings() const
   {
      return ((this->mScores.size() - this->mSkip) > 0);
   }

   void LoadSplitter::initAlgorithms(const ArrayI& dim)
   {
      // Check for serial version of code request
      if(this->nCpu() == 1)
      {
         #ifdef GEOMHDISCC_MPI
            // We don't deal with a single CPU in MPI mode
            throw Exception("LoadSplitter::initAlgorithms", "Tried to initialise (MPI) parallelisation with a single CPU");
         #else
            this->mAlgorithms.push_back(SharedSplittingAlgorithm(new SerialSplitting(this->id(), this->nCpu(), dim)));
         #endif //GEOMHDISCC_MPI

      // Setup the parallel version (initialise all algorithms and then choose the best one)
      } else
      {
         #ifdef GEOMHDISCC_MPI
            // Check that problem is at least 2D
            if(dim.size() > 1)
            {
               #ifdef GEOMHDISCC_MPIALGO_SINGLE1D
                  // Add the single splitting algorithm for first data exchange
                  this->mAlgorithms.push_back(SharedSplittingAlgorithm(new SingleSplitting(this->id(), this->nCpu(), dim, Splittings::Locations::FIRST)));
               #endif //GEOMHDISCC_MPIALGO_SINGLE1D

               // Check if problem is 3D
               if(dim.size() == 3)
               {
                  #ifdef GEOMHDISCC_MPIALGO_SINGLE2D
                     // Add the single splitting algorithm for second data exchange
                     this->mAlgorithms.push_back(SharedSplittingAlgorithm(new SingleSplitting(this->id(), this->nCpu(), dim, Splittings::Locations::SECOND)));
                  #endif //GEOMHDISCC_MPIALGO_SINGLE2D

                  #ifdef GEOMHDISCC_MPIALGO_TUBULAR
                     // Add the tubular splitting algorithm
                     this->mAlgorithms.push_back(SharedSplittingAlgorithm(new TubularSplitting(this->id(), this->nCpu(), dim)));
                  #endif //GEOMHDISCC_MPIALGO_TUBULAR
               }

            // 1D problems can't be parallelised with the current algorithms
            } else
            {
               throw Exception("LoadSplitter::initAlgorithms", "There is no parallelisation algorithm for 1D problems!");
            }
         #else
            // nCpu can't be bigger than one without MPI
            throw Exception("LoadSplitter::initAlgorithms", "Initialisation of nCPU > 1 without MPI is not possible!");
         #endif //GEOMHDISCC_MPI
      }

      // Safety check to make sure at least one algorithm has been initialised
      if(this->mAlgorithms.size() == 0)
      {
         throw Exception("LoadSplitter::initAlgorithms", "No algorithm has been initialised!");
      }
   }

   void LoadSplitter::initScores()
   {
      // Get iterator over splitting algorithms
      std::vector<SharedSplittingAlgorithm>::iterator it;

      // Loop over all initialised algorithms
      for(it = this->mAlgorithms.begin(); it != this->mAlgorithms.end(); it++)
      {
         // Loop over possible factorisations of nCPU
         while((*it)->useNextFactors())
         {
            // Check if obtained factorisation is applicable to problem
            if((*it)->applicable())
            {
               // Get scored resolution object
               std::pair<int, std::pair<SharedResolution, SplittingDescription> > scored = (*it)->scoreSplitting();

               // Only scores bigger than 0 are considered
               if(scored.first > 0)
               {
                  this->mScores.insert(scored);
               }
            }
         }
      }
   }

   std::pair<SharedResolution,SplittingDescription> LoadSplitter::bestSplitting() const
   {
      // Make sure there is at least one successful splitting (score > 0)
      if(this->mScores.size() > 0)
      {
         // Describe the splitting with the highest score
         this->describeSplitting(this->mScores.rbegin()->second.second);

         // Return the splitting with the highest score
         return this->mScores.rbegin()->second;
      } else
      {
         throw Exception("LoadSplitter::bestSplitting", "No usable splitting has been found!");
      }
   }

   void LoadSplitter::describeSplitting(const SplittingDescription& descr) const
   {
      // Output a short description of the selected splitting. Make it look nice ;)
      if(this->id() == 0)
      {
         // Print load splitting header
         FormatToolbox::printLine('-');
         FormatToolbox::printCentered("Load distribution", '*');
         FormatToolbox::printLine('-');

         std::string tmpStr;

         // Print grouper information
         switch(descr.grouper)
         {
            case(Splittings::Groupers::EQUATION):
               tmpStr = "Equation";
               break;
            case(Splittings::Groupers::SINGLE1D):
               tmpStr = "Single 1D";
               break;
            case(Splittings::Groupers::SINGLE2D):
               tmpStr = "Single 2D";
               break;
            case(Splittings::Groupers::TRANSFORM):
               tmpStr = "Transform";
               break;
         }
         FormatToolbox::printCentered("Grouper: " + tmpStr);

         // Print Algorithm information
         switch(descr.algorithm)
         {
            case(Splittings::Algorithms::SERIAL):
               tmpStr = "Serial";
               break;
            case(Splittings::Algorithms::SINGLE1D):
               tmpStr = "Single 1D";
               break;
            case(Splittings::Algorithms::SINGLE2D):
               tmpStr = "Single 2D";
               break;
            case(Splittings::Algorithms::TUBULAR):
               tmpStr = "Tubular";
               break;
         }
         FormatToolbox::printCentered("Algorithm: " + tmpStr);

         // Print factorisation information
         tmpStr = "";
         std::stringstream oss;
         for(int i = 0; i < descr.factors.size();i++)
         {
            oss << descr.factors(i);
            tmpStr += oss.str();
            oss.str("");
            if(i < descr.factors.size() - 1)
            {
               tmpStr += " x ";
            }
         }
         if(descr.algorithm == Splittings::Algorithms::SINGLE1D)
         {
            tmpStr += " x 1";
         } else if(descr.algorithm == Splittings::Algorithms::SINGLE2D)
         {
            tmpStr = "1 x " + tmpStr;
         }
         FormatToolbox::printCentered("Factorisation: " + tmpStr);

         oss.str("");
         oss << descr.score;
         FormatToolbox::printCentered("Score: " + oss.str());

         FormatToolbox::printNewline();
      }
   }

   std::pair<SharedResolution,SplittingDescription> LoadSplitter::nextSplitting()
   {
      // Make sure there is at least an additional splitting
      if(this->mScores.size() > this->mSkip)
      {
         // Create reverse iterator
         std::multimap<int, std::pair<SharedResolution,SplittingDescription> >::reverse_iterator   rit;

         // Set start iterator
         rit = this->mScores.rbegin();

         // Advance iterator
         std::advance(rit, this->mSkip);

         // Increment the skipped splitting counter
         this->mSkip++;

         // Describe the obtained splitting
         this->describeSplitting(rit->second.second);

         // Return the splitting
         return rit->second;
      } else
      {
         throw Exception("LoadSplitter::nextSplitting", "No more usable splitting are available!");
      }
   }

}
