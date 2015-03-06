/** 
 * @file LoadSplitter.cpp
 * @brief Source of the workload splitter
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Debug includes
//
#include "Exceptions/Exception.hpp"

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
#include "LoadSplitter/LoadSplitter.hpp"

// Project includes
//
#include "Enums/Splitting.hpp"
#include "IoTools/Formatter.hpp"

// Splitting algorithms
#include "LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "LoadSplitter/Algorithms/SerialSplitting.hpp"
#include "LoadSplitter/Algorithms/SingleSplitting.hpp"
#include "LoadSplitter/Algorithms/TubularSplitting.hpp"
#include "LoadSplitter/Algorithms/Coupled2DSplitting.hpp"
#include "IoXml/GxlWriter.hpp"
#include "IoXml/VtpWriter.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   LoadSplitter::LoadSplitter(const int id, const int nCpu)
      : mId(id), mNCpu(nCpu), mMaxNScores(1)
   {
   }

   LoadSplitter::~LoadSplitter()
   {
   }

   void LoadSplitter::initAlgorithms(const ArrayI& dim)
   {
      // Check for serial version of code request
      if(this->mNCpu == 1)
      {
         #ifdef GEOMHDISCC_MPI
            // We don't deal with a single CPU in MPI mode
            throw Exception("Tried to initialise (MPI) parallelisation with a single CPU");
         #else
            this->mAlgorithms.push_back(SharedSplittingAlgorithm(new SerialSplitting(this->mId, this->mNCpu, dim)));
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
                  this->mAlgorithms.push_back(SharedSplittingAlgorithm(new SingleSplitting(this->mId, this->mNCpu, dim, Splitting::Locations::FIRST)));
               #endif //GEOMHDISCC_MPIALGO_SINGLE1D

               #ifdef GEOMHDISCC_MPIALGO_COUPLED2D
                  // Add the single splitting algorithm for first data exchange
                  this->mAlgorithms.push_back(SharedSplittingAlgorithm(new Coupled2DSplitting(this->mId, this->mNCpu, dim)));
               #endif //GEOMHDISCC_MPIALGO_COUPLED2D


               // Check if problem is 3D
               if(dim.size() == 3)
               {
                  #ifdef GEOMHDISCC_MPIALGO_SINGLE2D
                     // Add the single splitting algorithm for second data exchange
                     this->mAlgorithms.push_back(SharedSplittingAlgorithm(new SingleSplitting(this->mId, this->mNCpu, dim, Splitting::Locations::SECOND)));
                  #endif //GEOMHDISCC_MPIALGO_SINGLE2D

                  #ifdef GEOMHDISCC_MPIALGO_TUBULAR
                     // Add the tubular splitting algorithm
                     this->mAlgorithms.push_back(SharedSplittingAlgorithm(new TubularSplitting(this->mId, this->mNCpu, dim)));
                  #endif //GEOMHDISCC_MPIALGO_TUBULAR
               }

            // 1D problems can't be parallelised with the current algorithms
            } else
            {
               throw Exception("There is no parallelisation algorithm for 1D problems!");
            }
         #else
            // nCpu can't be bigger than one without MPI
            throw Exception("Initialisation of nCPU > 1 without MPI is not possible!");
         #endif //GEOMHDISCC_MPI
      }

      // Safety check to make sure at least one algorithm has been initialised
      if(this->mAlgorithms.size() == 0)
      {
         throw Exception("No algorithm has been initialised!");
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
            // Check if obtained factorisation is applicable to splitting algorithm
            if((*it)->applicable())
            {
               // Get scored resolution object
               std::pair<int, std::pair<SharedResolution, SplittingDescription> > scored = (*it)->scoreSplitting();

               // Only scores bigger than 0 are considered
               if(scored.first > 0)
               {
                  // If not at minimal size yet, add to scores
                  if(this->mScores.size() < static_cast<unsigned int>(this->mMaxNScores))
                  {
                     this->mScores.insert(scored);
                  } else
                  {
                     if(scored.first > this->mScores.begin()->first)
                     {
                        // Remove lowest score
                        this->mScores.erase(this->mScores.begin());
                        // Add new score
                        this->mScores.insert(scored);
                     }
                  }

                  #ifdef GEOMHDISCC_DEBUG
                  // Show description of tested splittings
                  this->describeSplitting(scored.second.second, true);
                  #endif // GEOMHDISCC_DEBUG
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

         #ifdef GEOMHDISCC_DEBUG
            for(std::vector<IoXml::SharedVtpWriter>::const_iterator it = this->mScores.rbegin()->second.second.vtpFiles.begin(); it != this->mScores.rbegin()->second.second.vtpFiles.end(); ++it)
            {
               (*it)->write();
               (*it)->finalize();
            }
         #endif //GEOMHDISCC_DEBUG

         // Return the splitting with the highest score
         return this->mScores.rbegin()->second;
      } else
      {
         throw Exception("No usable splitting has been found!");
      }
   }

   void LoadSplitter::describeSplitting(const SplittingDescription& descr, const bool isTest) const
   {
      // Output a short description of the selected splitting. Make it look nice ;)
      if(this->mId == 0)
      {
         if(isTest)
         {
            // Print load splitting test header
            IoTools::Formatter::printLine(std::cout, '~');
            IoTools::Formatter::printCentered(std::cout, "Tested load distribution", '~');
            IoTools::Formatter::printLine(std::cout, '~');
         } else
         {
            // Print load splitting header
            IoTools::Formatter::printLine(std::cout, '-');
            IoTools::Formatter::printCentered(std::cout, "Load distribution", '*');
            IoTools::Formatter::printLine(std::cout, '-');
         }

         std::string tmpStr;

         // Print grouper information
         switch(descr.grouper)
         {
            case(Splitting::Groupers::EQUATION):
               tmpStr = "Equation";
               break;
            case(Splitting::Groupers::SINGLE1D):
               tmpStr = "Single 1D";
               break;
            case(Splitting::Groupers::SINGLE2D):
               tmpStr = "Single 2D";
               break;
            case(Splitting::Groupers::TRANSFORM):
               tmpStr = "Transform";
               break;
         }
         IoTools::Formatter::printCentered(std::cout, "Grouper: " + tmpStr);

         // Print Algorithm information
         switch(descr.algorithm)
         {
            case(Splitting::Algorithms::SERIAL):
               tmpStr = "Serial";
               break;
            case(Splitting::Algorithms::SINGLE1D):
               tmpStr = "Single 1D";
               break;
            case(Splitting::Algorithms::SINGLE2D):
               tmpStr = "Single 2D";
               break;
            case(Splitting::Algorithms::TUBULAR):
               tmpStr = "Tubular";
               break;
            case(Splitting::Algorithms::COUPLED2D):
               tmpStr = "Coupled 2D";
               break;
         }
         IoTools::Formatter::printCentered(std::cout, "Algorithm: " + tmpStr);

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
         if(descr.algorithm == Splitting::Algorithms::SINGLE1D || descr.algorithm == Splitting::Algorithms::COUPLED2D)
         {
            tmpStr += " x 1";
         } else if(descr.algorithm == Splitting::Algorithms::SINGLE2D)
         {
            tmpStr = "1 x " + tmpStr;
         }
         IoTools::Formatter::printCentered(std::cout, "Factorisation: " + tmpStr);

         // Create GXL graph format file
         tmpStr.erase(std::remove(tmpStr.begin(),tmpStr.end(),' '),tmpStr.end());
         IoXml::GxlWriter gxl("Communication_graph_" + tmpStr);
         gxl.init();
         gxl.graphCommunication(descr.structure);
         gxl.write();
         gxl.finalize();

         oss.str("");
         oss << static_cast<int>(descr.score.prod());
         IoTools::Formatter::printCentered(std::cout, "Score: " + oss.str());
         oss.str("");
         oss << " (" << descr.score(0) << ", " << std::setprecision(2)<< descr.score(1) << ", " << descr.score(2) << ", " << descr.score(3) << ")";
         IoTools::Formatter::printCentered(std::cout, oss.str());

         IoTools::Formatter::printNewline(std::cout);
      }
   }

   void LoadSplitter::showSplittings(const int n) const
   {
      // Get maximum between number of scores and n
      int maxN = std::min(static_cast<int>(this->mScores.size()), n);

      // Create reverse iterator
      std::multimap<int, std::pair<SharedResolution,SplittingDescription> >::const_reverse_iterator rit;

      // Set start iterator
      rit = this->mScores.rbegin();

      // Loop over scores
      for(int i = 0; i < maxN; ++i)
      {
         // Describe the obtained splitting
         this->describeSplitting(rit->second.second);

         // Increment iterator
         rit++;
      }
   }

}
}
