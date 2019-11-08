/** 
 * @file ShellTorPolDissipationWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics dissipation calculation for vectorial field in a spherical shell
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>
#include <cmath>
#include <iostream>

// External includes
//

// Class include
//
#include "IoVariable/ShellTorPolDissipationWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/DissipationTags.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace QuICC {

   namespace IoVariable {

      ShellTorPolDissipationWriter::ShellTorPolDissipationWriter(const std::string &prefix, const std::string &type)
              : IVariableAsciiEWriter(prefix + DissipationTags::BASENAME, DissipationTags::EXTENSION,
                                      prefix + DissipationTags::HEADER, type, DissipationTags::VERSION,
                                      Dimensions::Space::SPECTRAL), mTorDiss(-1.0), mPolDiss(-1.0) {
      }

      ShellTorPolDissipationWriter::~ShellTorPolDissipationWriter() {
      }

      void ShellTorPolDissipationWriter::init() {
         // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
         MHDFloat ro = this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second;
         MHDFloat ri = ro * this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second;
         this->mVolume = (4.0 / 3.0) * Math::PI * (std::pow(ro, 3) - std::pow(ri, 3));

         // Initialise python wrapper
         PythonWrapper::init();
         PythonWrapper::import("quicc.geometry.spherical.shell_radius");

         // Prepare arguments
         PyObject *pArgs, *pValue;
         pArgs = PyTuple_New(4);

         // ... compute a, b factors
         PyObject *pTmp = PyTuple_New(2);
         PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(
                 this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RO))->second));
         PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(
                 this->mPhysical.find(IoTools::IdToHuman::toTag(NonDimensional::RRATIO))->second));
         PythonWrapper::setFunction("linear_r2x");
         pValue = PythonWrapper::callFunction(pTmp);
         PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue, 0));
         PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue, 1));
         // ... create boundray condition (none)
         pValue = PyDict_New();
         PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(0));
         PyTuple_SetItem(pArgs, 3, pValue);

         // Get resolution
         int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>() + 2;
         pValue = PyLong_FromLong(cols);
         PyTuple_SetItem(pArgs, 0, pValue);

         // Call r^2
         PythonWrapper::setFunction("r2");
         pValue = PythonWrapper::callFunction(pArgs);
         // Fill matrix and cleanup
         SparseMatrix tmpR2(cols, cols);
         PythonWrapper::fillMatrix(tmpR2, pValue);
         Py_DECREF(pValue);

         pTmp = PyTuple_GetSlice(pArgs, 0, 3);
         // Call avg
         PythonWrapper::setFunction("integral");
         pValue = PythonWrapper::callFunction(pTmp);
         // Fill matrix and cleanup
         SparseMatrix tmpAvg(cols, cols);
         PythonWrapper::fillMatrix(tmpAvg, pValue);
         Py_DECREF(pValue);
         PythonWrapper::finalize();

         // Store integral
         this->mIntgOp = tmpAvg.leftCols(cols - 2);

         // Store spherical integral (include r^2 factor)
         this->mSphIntgOp = tmpAvg * tmpR2.leftCols(cols - 2);

         IVariableAsciiEWriter::init();
      }

      void ShellTorPolDissipationWriter::compute(Transform::TransformCoordinatorType &coord) {
         // get iterator to field
         vector_iterator vIt;
         vector_iterator_range vRange = this->vectorRange();
         assert(std::distance(vRange.first, vRange.second) == 1);
         assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
         assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

         // Compute integral over Chebyshev expansion and sum harmonics
         this->mTorDiss = 0.0;
         this->mPolDiss = 0.0;
         this->mCentroAntysymDiss = 0.0;
         this->mCentroSymDiss = 0.0;
         this->mEquaAntysymDiss = 0.0;
         this->mEquaSymDiss = 0.0;

         // Dealias poloidal variable data
         coord.communicator().dealiasSpectral(
                 vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolRadLapl = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPolRadLapl2 = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension
         coord.transform1D().project(rOutVarPolRadLapl.rData(), rInVarPol.data(),
                                     Transform::TransformCoordinatorType::Transform1DType::ProjectorType::RADLAPL);

         // Compute |f|^2
         rOutVarPolRadLapl2.rData() =
                 rOutVarPolRadLapl.rData().array() * rOutVarPolRadLapl.rData().conjugate().array();

         // Compute projection transform for first dimension
         coord.transform1D().integrate_full(rInVarPol.rData(), rOutVarPolRadLapl2.data(),
                                            Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

         MHDFloat lfactor = 0.0;
         MHDFloat factor = 0.0;

#ifdef QUICC_SPATIALSCHEME_SLFM
         // Loop over harmonic order m
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if( m == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = l*(l+1.0);

                MHDFloat ModeDiss = factor*lfactor*(this->mSphIntgOp*rInVarPol.slice(k).col(j).real()).sum();
                this->mPolDiss += ModeDiss;


                    // assign the centro symmetry energies
                    if( (l % 2) == 0 )
                    {
                        this->mCentroSymDiss += ModeDiss;
                    } else {

                        this->mCentroAntysymDiss += ModeDiss;
                    }

                    // assign the equatorial symmetry energies
                    if( ((l+m) % 2) == 0 )
                    {
                        this->mEquaSymDiss += ModeDiss;
                    } else {

                        this->mEquaAntysymDiss += ModeDiss;
                    }

            }
         }
#endif //defined QUICC_SPATIALSCHEME_SLFM
#ifdef QUICC_SPATIALSCHEME_SLFL
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = l*(l+1.0);
            int start = 0;

            int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
            // m = 0, no factor of two
            if( firstM == 0)
            {
               MHDFloat ModeDiss = lfactor*(this->mSphIntgOp*rInVarPol.slice(k).col(0).real()).sum();
               this->mPolDiss += ModeDiss;
               start = 1;
               firstM = 1;


               if ( (l % 2) == 0){

                  this->mCentroSymDiss += ModeDiss;
                  //this->mEquaSymDiss += ModeDiss;
               } else {

                  this->mCentroAntysymDiss += ModeDiss;
                  //this->mEquaAntysymDiss += ModeDiss;
               }

            }

            Matrix MatrixModes = 2.0*lfactor*(this->mSphIntgOp*rInVarPol.slice(k).rightCols(rInVarPol.slice(k).cols()-start).real());
            this->mPolDiss += MatrixModes.sum();

            if ( (l % 2) == 0){

              this->mCentroSymDiss += MatrixModes.sum();
            } else {

              this->mCentroAntysymDiss += MatrixModes.sum();
            }

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){

               int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               if(m==0){
                  factor = 1.0;
               } else {
                  factor = 2.0;
               }

               MHDFloat temp = factor*lfactor*(this->mSphIntgOp*rInVarPol.slice(k).col(j).real()).sum();
               if( ((l+m) % 2) == 0)
				   {
                  this->mEquaSymDiss += temp;

				   } else {
                  this->mEquaAntysymDiss += temp;
				   }

            }



         }
#endif //QUICC_SPATIALSCHEME_SLFL

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPol);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolRadLapl2);


         // Dealias poloidal variable data
         coord.communicator().dealiasSpectral(
                 vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

         // Recover dealiased BWD data
         rInVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension
         coord.transform1D().project(rOutVarPol.rData(), rInVarPol.data(),
                                     Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

         // Compute |f|^2
         rOutVarPol.rData() = rOutVarPolRadLapl.rData().array() * rOutVarPol.rData().conjugate().array();
         //
         //rOutVarPol.rData() = rOutVarPolRadLapl.rData().array()*rOutVarPol.rData().conjugate().array() + rOutVarPolRadLapl.rData().conjugate().array()*rOutVarPol.rData().array();

         // Compute projection transform for first dimension
         coord.transform1D().integrate_full(rInVarPol.rData(), rOutVarPol.data(),
                                            Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

         lfactor = 0.0;

#ifdef QUICC_SPATIALSCHEME_SLFM
           // Loop over harmonic order m
           for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
           {
              // m = 0, no factor of two
                  int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
              if( m == 0)
              {
                 factor = 1.0;
              } else
              {
                 factor = 2.0;
              }

              for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
              {
                 int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
                 lfactor = std::pow(l*(l+1.0),2);

                  MHDFloat ModeDiss = -2.*factor*lfactor*(this->mIntgOp*rInVarPol.slice(k).col(j).real()).sum();
                 this->mPolDiss += ModeDiss;


                      // assign the centro symmetry energies
                      if( (l % 2) == 0 )
                      {
                          this->mCentroSymDiss += ModeDiss;
                      } else {

                          this->mCentroAntysymDiss += ModeDiss;
                      }

                      // assign the equatorial symmetry energies
                      if( ((l+m) % 2) == 0 )
                      {
                          this->mEquaSymDiss += ModeDiss;
                      } else {

                          this->mEquaAntysymDiss += ModeDiss;
                      }

              }
           }
#endif //defined QUICC_SPATIALSCHEME_SLFM
#ifdef QUICC_SPATIALSCHEME_SLFL
         // Loop over harmonic degree l
           for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
           {
              int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
              lfactor = std::pow(l*(l+1.0),2);
              int start = 0;

              int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
              // m = 0, no factor of two
              if( firstM == 0)
              {
                 MHDFloat ModeDiss = -2.0*lfactor*(this->mIntgOp*rInVarPol.slice(k).col(0).real()).sum();
                 this->mPolDiss += ModeDiss;
                 start = 1;
                 firstM = 1;


                 if ( (l % 2) == 0){

                    this->mCentroSymDiss += ModeDiss;
                    //this->mEquaSymDiss += ModeDiss;
                 } else {

                    this->mCentroAntysymDiss += ModeDiss;
                    //this->mEquaAntysymDiss += ModeDiss;
                 }

              }

              Matrix MatrixModes = -2.0*2.0*lfactor*(this->mIntgOp*rInVarPol.slice(k).rightCols(rInVarPol.slice(k).cols()-start).real());
              this->mPolDiss += MatrixModes.sum();

              if ( (l % 2) == 0){

                      this->mCentroSymDiss += MatrixModes.sum();
                  } else {

                      this->mCentroAntysymDiss += MatrixModes.sum();
                  }


               for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){

                  int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
                  if(m==0){
                     factor = 1.0;
                  } else {
                     factor = 2.0;
                  }

                  MHDFloat temp = -2.0*factor*lfactor*(this->mIntgOp*rInVarPol.slice(k).col(j).real()).sum();
                  if( ((l+m) % 2) == 0)
                  {
                     this->mEquaSymDiss += temp;

                  } else {
                     this->mEquaAntysymDiss += temp;
                  }

               }



           }
#endif //QUICC_SPATIALSCHEME_SLFL

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPol);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPolRadLapl);
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPol);


         // Dealias poloidal variable data
         coord.communicator().dealiasSpectral(
                 vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

         // Recover dealiased BWD data
         rInVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         rOutVarPol = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension
         coord.transform1D().project(rOutVarPol.rData(), rInVarPol.data(),
                                     Transform::TransformCoordinatorType::Transform1DType::ProjectorType::DIVR);

         // Compute |f|^2
         rOutVarPol.rData() = rOutVarPol.rData().array() * rOutVarPol.rData().conjugate().array();

         // Compute projection transform for first dimension
         coord.transform1D().integrate_full(rInVarPol.rData(), rOutVarPol.data(),
                                            Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

         lfactor = 0.0;

#ifdef QUICC_SPATIALSCHEME_SLFM
           // Loop over harmonic order m
           for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
           {
              // m = 0, no factor of two
                  int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
              if( m == 0)
              {
                 factor = 1.0;
              } else
              {
                 factor = 2.0;
              }

              for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
              {
                 int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
                 lfactor = std::pow(l*(l+1.0),3);

                  MHDFloat ModeDiss = factor*lfactor*(this->mIntgOp*rInVarPol.slice(k).col(j).real()).sum();
                 this->mPolDiss += ModeDiss;


                      // assign the centro symmetry energies
                      if( (l % 2) == 0 )
                      {
                          this->mCentroSymDiss += ModeDiss;
                      } else {

                          this->mCentroAntysymDiss += ModeDiss;
                      }

                      // assign the equatorial symmetry energies
                      if( ((l+m) % 2) == 0 )
                      {
                          this->mEquaSymDiss += ModeDiss;
                      } else {

                          this->mEquaAntysymDiss += ModeDiss;
                      }

              }
           }
#endif //defined QUICC_SPATIALSCHEME_SLFM
#ifdef QUICC_SPATIALSCHEME_SLFL
         // Loop over harmonic degree l
           for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
           {
              int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
              lfactor = std::pow(l*(l+1.0),3);
              int start = 0;

              int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
              // m = 0, no factor of two
              if( firstM == 0)
              {
                 MHDFloat ModeDiss = lfactor*(this->mIntgOp*rInVarPol.slice(k).col(0).real()).sum();
                 this->mPolDiss += ModeDiss;
                 start = 1;
                 firstM = 1;

                 if ( (l % 2) == 0){

                    this->mCentroSymDiss += ModeDiss;
                    //this->mEquaSymDiss += ModeDiss;
                 } else {

                    this->mCentroAntysymDiss += ModeDiss;
                    //this->mEquaAntysymDiss += ModeDiss;
                 }

              }

              Matrix MatrixModes = 2.0*lfactor*(this->mIntgOp*rInVarPol.slice(k).rightCols(rInVarPol.slice(k).cols()-start).real());
              this->mPolDiss += MatrixModes.sum();

            if ( (l % 2) == 0){

                this->mCentroSymDiss += MatrixModes.sum();
            } else {

                this->mCentroAntysymDiss += MatrixModes.sum();
            }

               for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){

                  int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
                  if(m==0){
                     factor = 1.0;
                  } else {
                     factor = 2.0;
                  }

                  MHDFloat temp = factor*lfactor*(this->mIntgOp*rInVarPol.slice(k).col(j).real()).sum();
                  if( ((l+m) % 2) == 0)
                  {
                     this->mEquaSymDiss += temp;

                  } else {
                     this->mEquaAntysymDiss += temp;
                  }

               }



           }
#endif //QUICC_SPATIALSCHEME_SLFL

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPol);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarPol);

         // Normalize by sphere volume: 4/3*pi*(r_o^3 - r_i^3)
         //this->mPolDiss /= 2 * this->mVolume;

         // Dealias toroidal variable data for Q component
         coord.communicator().dealiasSpectral(
                 vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTorQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarTorQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension
         coord.transform1D().project(rOutVarTorQ.rData(), rInVarTorQ.data(),
                                     Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

         // Compute |f|^2
         rOutVarTorQ.rData() = rOutVarTorQ.rData().array() * rOutVarTorQ.rData().conjugate().array();

         // Compute projection transform for first dimension
         coord.transform1D().integrate_full(rInVarTorQ.rData(), rOutVarTorQ.data(),
                                            Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

         // Compute Diss in Q component of QST decomposition (toroidal component)
#ifdef QUICC_SPATIALSCHEME_SLFM
         // Loop over harmonic order m
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
                int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if( m == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = std::pow(l*(l+1.0),2);

               MHDFloat ModeDiss = factor*lfactor*(this->mIntgOp*rInVarTorQ.slice(k).col(j).real()).sum();
               this->mTorDiss += ModeDiss;


               // assign the centro symmetry energies
               if( (l % 2) == 1 )
               {

                 this->mCentroSymDiss += ModeDiss;
               } else {

                 this->mCentroAntysymDiss += ModeDiss;
               }

               // assign the equatorial symmetry energies
               if( ((l+m) % 2) == 1 )
               {
                 this->mEquaSymDiss += ModeDiss;
               } else {

                 this->mEquaAntysymDiss += ModeDiss;
               }

            }
         }
#endif //defined QUICC_SPATIALSCHEME_SLFM
#ifdef QUICC_SPATIALSCHEME_SLFL
         lfactor = 0.0;
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = std::pow(l*(l+1.0),2);
            int start = 0;
            // m = 0, no factor of two

            int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
            if( firstM == 0)
            {
               MHDFloat ModeDiss = lfactor*(this->mIntgOp*rInVarTorQ.slice(k).col(0).real()).sum();
               this->mTorDiss += ModeDiss;
               start = 1;
               firstM = 1;

               if ( (l % 2) == 1)
               {
                   this->mCentroSymDiss += ModeDiss;
                   //this->mEquaSymDiss += ModeDiss;
               } else {

                    this->mCentroAntysymDiss += ModeDiss;
                    //this->mEquaAntysymDiss += ModeDiss;
               }
            }
            Matrix MatrixModes = 2.0*lfactor*(this->mIntgOp*rInVarTorQ.slice(k).rightCols(rInVarTorQ.slice(k).cols()-start).real());
            this->mTorDiss += MatrixModes.sum();

            if ( (l % 2) == 1 )
            {
                this->mCentroSymDiss += MatrixModes.sum();
            } else {
                this->mCentroAntysymDiss += MatrixModes.sum();
            }

               for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){

                  int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
                  if(m==0){
                     factor = 1.0;
                  } else {
                     factor = 2.0;
                  }

                  MHDFloat temp = factor*lfactor*(this->mIntgOp*rInVarTorQ.slice(k).col(j).real()).sum();
                  if( ((l+m) % 2) == 1)
                  {
                     this->mEquaSymDiss += temp;

                  } else {
                     this->mEquaAntysymDiss += temp;
                  }

               }



         }
#endif //QUICC_SPATIALSCHEME_SLFL

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTorQ);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTorQ);

         // Dealias toroidal variable data for S component
         coord.communicator().dealiasSpectral(
                 vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTorS = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarTorS = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension
         coord.transform1D().project(rOutVarTorS.rData(), rInVarTorS.data(),
                                     Transform::TransformCoordinatorType::Transform1DType::ProjectorType::DIFFR);

         // Compute |f|^2
         rOutVarTorS.rData() = rOutVarTorS.rData().array() * rOutVarTorS.rData().conjugate().array();

         // Compute projection transform for first dimension
         coord.transform1D().integrate_full(rInVarTorS.rData(), rOutVarTorS.data(),
                                            Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

         // Compute Diss in S component of QST decomposition
#ifdef QUICC_SPATIALSCHEME_SLFM
         // Loop over harmonic order m
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
                int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if(m == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l= this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = l*(l+1.0);

               MHDFloat ModeDiss = factor*lfactor*(this->mIntgOp*rInVarTorS.slice(k).col(j).real()).sum();
               this->mTorDiss += ModeDiss;

               // assign the centro symmetry energies
               if( (l % 2) == 1 )
               {
                   this->mCentroSymDiss += ModeDiss;
               } else {

                   this->mCentroAntysymDiss += ModeDiss;
               }


               // assign the equatorial symmetry energies
               if( ((l+m) % 2) == 1 )
               {

               this->mEquaSymDiss += ModeDiss;
               } else {

               this->mEquaAntysymDiss += ModeDiss;
               }


            }
         }
#endif //defined QUICC_SPATIALSCHEME_SLFM
#ifdef QUICC_SPATIALSCHEME_SLFL
         lfactor = 0.0;
         // Loop over harmonic degree l
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = l*(l+1.0);
            int start = 0;
            // m = 0, no factor of two
            int firstM = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k);
            if(firstM == 0)
            {
               MHDFloat ModeDiss =  lfactor*(this->mIntgOp*rInVarTorS.slice(k).col(0).real()).sum();
               this->mTorDiss += ModeDiss;
               start = 1;
               firstM = 1;


               if ( (l % 2) == 1 ){

                  this->mCentroSymDiss += ModeDiss;
                  //this->mEquaSymDiss += ModeDiss;
               } else {

                  this->mCentroAntysymDiss += ModeDiss;
                  //this->mEquaAntysymDiss += ModeDiss;
               }

            }



            Matrix MatrixModes = 2.0*lfactor*(this->mIntgOp*rInVarTorS.slice(k).rightCols(rInVarTorS.slice(k).cols()-start).real());
            this->mPolDiss += MatrixModes.sum();

                if ( (l % 2) == 1 ){

                    this->mCentroSymDiss += MatrixModes.sum();
                } else {

                    this->mCentroAntysymDiss += MatrixModes.sum();
                }

               for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){

                  int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
                  if(m==0){
                     factor = 1.0;
                  } else {
                     factor = 2.0;
                  }

                  MHDFloat temp = factor*lfactor*(this->mIntgOp*rInVarTorS.slice(k).col(j).real()).sum();
                  if( ((l+m) % 2) == 1)
                  {
                     this->mEquaSymDiss += temp;

                  } else {
                     this->mEquaAntysymDiss += temp;
                  }

               }

         }
#endif //QUICC_SPATIALSCHEME_SLFL

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTorS);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTorS);

         // Normalize by the volume
         this->mPolDiss /= 2 * this->mVolume;
         this->mTorDiss /= 2 * this->mVolume;

         // Normalize the remainder of the energies by volume
         this->mEquaSymDiss /= 2 * this->mVolume;
         this->mEquaAntysymDiss /= 2 * this->mVolume;
         this->mCentroSymDiss /= 2 * this->mVolume;
         this->mCentroAntysymDiss /= 2 * this->mVolume;

      }

      void ShellTorPolDissipationWriter::write() {
         // Create file
         this->preWrite();

         // Get the "global" Kinetic Diss from MPI code
#ifdef QUICC_MPI
         Array Diss(6);

         Diss(0) = this->mTorDiss;
         Diss(1) = this->mPolDiss;
         Diss(2) = this->mCentroSymDiss;
         Diss(3) = this->mCentroAntysymDiss;
         Diss(4) = this->mEquaSymDiss;
         Diss(5) = this->mCentroAntysymDiss;


         MPI_Allreduce(MPI_IN_PLACE, Diss.data(), 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorDiss = Diss(0);
         this->mPolDiss = Diss(1);
           this->mCentroSymDiss = Diss(2);
         this->mCentroAntysymDiss = Diss(3);
         this->mEquaSymDiss = Diss(4);
         this->mEquaAntysymDiss = Diss(5);
#endif //QUICC_MPI

         // Check if the workflow allows IO to be performed
         if (FrameworkMacro::allowsIO()) {
            this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorDiss + this->mPolDiss << "\t"
                        << this->mTorDiss << "\t" << this->mPolDiss << "\t"
                        << this->mCentroSymDiss << "\t" << this->mCentroAntysymDiss << "\t" << this->mEquaSymDiss
                        << "\t" << this->mEquaAntysymDiss << std::endl;
         }

         // Close file
         this->postWrite();

         // Abort if kinetic Diss is NaN
         if (std::isnan(this->mTorDiss) || std::isnan(this->mPolDiss)) {
            FrameworkMacro::abort(99);

            throw Exception("Toroidal/Poloidal Diss is NaN!");
         }
      }

   }
};