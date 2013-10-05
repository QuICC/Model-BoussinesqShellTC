/** 
 * @file EquationTools.cpp
 * @brief Source of the requirement tools to work with variables and equations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>
#include <algorithm>

// External includes
//

// Class include
//
#include "Equations/Tools/EquationTools.hpp"

// Project includes
//
#include "Debug/DebuggerMacro.h"

namespace GeoMHDiSCC {

namespace Equations {

namespace Tools {

   void identifySolver(const std::pair<std::vector<SharedIScalarEquation>::iterator,std::vector<SharedIScalarEquation>::iterator>& scalarRange, const std::pair<std::vector<SharedIVectorEquation>::iterator,std::vector<SharedIVectorEquation>::iterator>& vectorRange)
   {
      // Iterators for scalar equations
      std::vector<SharedIScalarEquation>::iterator  scalEqIt;
      std::vector<SharedIScalarEquation>::iterator  doneSEqIt;

      // Current solver indexes for real and complex solvers
      int dIdx = 0;
      int zIdx = 0;

      // Coupling flag
      bool coupled = false;

      // Field identification
      FieldComponents::Spectral::Id compId = FieldComponents::Spectral::SCALAR;
      SpectralFieldId   fieldId;
      
      // Loop over the scalar equations
      for(scalEqIt = scalarRange.first; scalEqIt != scalarRange.second; ++scalEqIt)
      {
         // Build field identity
         fieldId = std::make_pair((*scalEqIt)->name(), compId);

         // Loop over already identified equations
         for(doneSEqIt = scalarRange.first; doneSEqIt != scalEqIt; ++doneSEqIt)
         {
            // loop over the implicit range
            CouplingInformation::FieldId_iterator fIt;
            CouplingInformation::FieldId_range fRange = (*doneSEqIt)->couplingInfo(compId).implicitRange();
            for(fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               // Check if field is in implicit range
               if(*fIt == fieldId)
               {
                  // Set the solver index
                  (*scalEqIt)->setSolverIndex(compId, (*doneSEqIt)->couplingInfo(compId).solverIndex());
                  
                  // Debug statements
                  DebuggerMacro_showValue("Identified coupled scalar solver: ", 2, (*scalEqIt)->name());
                  DebuggerMacro_showValue("---> solver index: ", 2, (*scalEqIt)->couplingInfo(compId).solverIndex());
                  DebuggerMacro_showValue("---> is complex? ", 2, (*scalEqIt)->couplingInfo(compId).isComplex());
                  
                  // Set coupling flag and break out
                  coupled = true;
                  break;
               }
            }

            // Break out of loop if already identified
            if(coupled)
            {
               break;
            }
         }

         // All checked and equation is not coupled
         if(!coupled)
         {
            // Set solver index for complex equation
            if((*scalEqIt)->couplingInfo(compId).isComplex())
            {
               // Set complex solver index
               (*scalEqIt)->setSolverIndex(compId, zIdx);

               // Increment complex solver index
               zIdx++;

            // Set solver index for real equation
            } else
            {
               // Set real solver index
               (*scalEqIt)->setSolverIndex(compId, dIdx);
               
               // Increment real solver index
               dIdx++;
            }

            // Debug statements
            DebuggerMacro_showValue("Identified first scalar solver: ", 2, (*scalEqIt)->name());
            DebuggerMacro_showValue("---> solver index: ", 2, (*scalEqIt)->couplingInfo(compId).solverIndex());
            DebuggerMacro_showValue("---> is complex? ", 2, (*scalEqIt)->couplingInfo(compId).isComplex());
         }

         // Reset coupling flag
         coupled = false;
      }

      // Iterators for vector equations
      std::vector<SharedIVectorEquation>::iterator  vectEqIt;
      std::vector<SharedIVectorEquation>::iterator  doneVEqIt;
      
      // Loop over the vector equations
      for(vectEqIt = vectorRange.first; vectEqIt != vectorRange.second; ++vectEqIt)
      {
         // Get coupled counter for each component
         ArrayI counter((*vectEqIt)->nSpectral());
         counter.setConstant(0);

         // Loop over the (identified) scalar equations
         for(doneSEqIt = scalarRange.first; doneSEqIt != scalarRange.second; ++doneSEqIt)
         {
            // loop over the implicit range
            CouplingInformation::FieldId_iterator fIt;
            CouplingInformation::FieldId_range fRange = (*doneSEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
            for(fIt = fRange.first; fIt != fRange.second; ++fIt)
            {
               IVectorEquation::SpectralComponent_iterator compIt;
               IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
               int i = 0;
               for(compIt = compRange.first; compIt != compRange.second; ++compIt, ++i)
               {
                  // Check if field's first component is in implicit range
                  if(counter(i) == 0 && *fIt == std::make_pair((*vectEqIt)->name(), *compIt))
                  {
                     // Set the solver index
                     (*vectEqIt)->setSolverIndex(*compIt, (*doneSEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).solverIndex());

                     // Debug statements
                     DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                     DebuggerMacro_showValue("---> component: ", 2, *compIt);
                     DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(*compIt).solverIndex());
                     DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(*compIt).isComplex());

                     // Set coupling flag and break out
                     counter(i) = 1;
                  }
               }

               // Break out of loop if already identified
               if(counter.sum() == counter.size())
               {
                  break;
               }
            }

            // Break out of loop if already identified
            if(counter.sum() == counter.size())
            {
               break;
            }
         }

         if(counter.sum() < counter.size())
         {
            // Loop over the (identified) vector equations
            for(doneVEqIt = vectorRange.first; doneVEqIt != vectEqIt; ++doneVEqIt)
            {
               IVectorEquation::SpectralComponent_iterator doneIt;
               IVectorEquation::SpectralComponent_range  doneRange = (*doneVEqIt)->spectralRange();
               for(doneIt = doneRange.first; doneIt != doneRange.second; ++doneIt)
               {
                  // loop over the implicit range of identified components
                  CouplingInformation::FieldId_iterator fIt;
                  CouplingInformation::FieldId_range fRange = (*doneVEqIt)->couplingInfo(*doneIt).implicitRange();
                  for(fIt = fRange.first; fIt != fRange.second; ++fIt)
                  {
                     IVectorEquation::SpectralComponent_iterator compIt;
                     IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
                     int i = 0;
                     for(compIt = compRange.first; compIt != compRange.second; ++compIt, ++i)
                     {
                        // Check if field's first component is in implicit range
                        if(counter(i) == 0 && *fIt == std::make_pair((*vectEqIt)->name(), *compIt))
                        {
                           // Set the solver index
                           (*vectEqIt)->setSolverIndex(*compIt, (*doneVEqIt)->couplingInfo(*doneIt).solverIndex());

                           // Debug statements
                           DebuggerMacro_showValue("Identified coupled vector solver: ", 2, (*vectEqIt)->name());
                           DebuggerMacro_showValue("---> component: ", 2, *compIt);
                           DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(*compIt).solverIndex());
                           DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(*compIt).isComplex());

                           // Set coupling flag and break out
                           counter(i) = 1;
                        }
                     }


                     // Break out of loop if already identified
                     if(counter.sum() == counter.size())
                     {
                        break;
                     }
                  }

                  // Break out of loop if already identified
                  if(counter.sum() == counter.size())
                  {
                     break;
                  }
               }
            }
         }

         // All checked and equation is not coupled
         if(counter.sum() < counter.size())
         {
            IVectorEquation::SpectralComponent_iterator compIt;
            IVectorEquation::SpectralComponent_range  compRange = (*vectEqIt)->spectralRange();
            int i = 0;
            for(compIt = compRange.first; compIt != compRange.second; ++compIt, ++i)
            {
               if(counter(i) == 0)
               {
                  // Set solver index for complex equation
                  if((*vectEqIt)->couplingInfo(*compIt).isComplex())
                  {
                     // Set complex solver index
                     (*vectEqIt)->setSolverIndex(*compIt, zIdx);

                     // Increment complex solver index
                     zIdx++;

                     // Set solver index for real equation
                  } else
                  {
                     // Set real solver index
                     (*vectEqIt)->setSolverIndex(*compIt, dIdx);

                     // Increment real solver index
                     dIdx++;
                  }

                  // Debug statements
                  DebuggerMacro_showValue("Identified first vector solver: ", 2, (*vectEqIt)->name());
                  DebuggerMacro_showValue("---> component: ", 2, *compIt);
                  DebuggerMacro_showValue("---> solver index: ", 2, (*vectEqIt)->couplingInfo(*compIt).solverIndex());
                  DebuggerMacro_showValue("---> is complex? ", 2, (*vectEqIt)->couplingInfo(*compIt).isComplex());
               }
            }
         }

         // Reset coupling flag
         counter.setConstant(0);
      }
   }

   SparseMatrix makeBlockMatrix(const int nBlocks, const int row, const int col)
   {
      SparseMatrix  mat(nBlocks, nBlocks);

      mat.insert(row, col) = 1;
      mat.makeCompressed();
      
      return mat;
   }
}
}
}
