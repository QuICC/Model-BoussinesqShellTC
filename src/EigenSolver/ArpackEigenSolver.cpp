/** \file ArpackEigenSolver.cpp
 *  \brief Implementation of an eigensolver based on ARPACK
 */

// System includes
//
#include <cassert>
#include <algorithm>

// External includes
//

// Class include
//
#include "EigenSolver/ArpackEigenSolver.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "../External/Interfaces/ARPACK_Interface.h"

#include <iostream>
namespace GeoMHDiSCC {

   ArpackEigenSolver::ArpackEigenSolver()
      : mSigma(0.0,0.0), mTol(0.0), mMaxIter(300), mNcv(20), mWhich("LM"), mIpntr(14), mIparam(11)
   {
   }

   ArpackEigenSolver::~ArpackEigenSolver()
   {
   }

   void ArpackEigenSolver::compute(const SparseMatrixZ& lhs, const SparseMatrixZ& rhs)
   {
      // Safety asserts for matrix sizes
      assert(lhs.rows() == lhs.cols());
      assert(rhs.rows() == rhs.cols());
      assert(lhs.rows() == rhs.rows());

      // Create LHS and RHS operators for the eigen solver
      this->mLhs = lhs - this->mSigma*rhs;
      this->mRhs = rhs;

      // Compute factorisation of LHS operator
      this->mSolver.compute(this->mLhs);

      // Check for successful factorisation
      if(this->mSolver.info() != Eigen::Success)
      {
         throw Exception("Factorisation of eigen solver operator failed!");
      }

      // Allocate work memory
      this->allocateMemory();
   }

   void ArpackEigenSolver::allocateMemory()
   {
      this->mResid.resize(this->mLhs.rows());
      this->mV.resize(this->mLhs.rows(), this->mNcv);
      this->mWorkd.resize(3*this->mLhs.rows());
      this->mWorkl.resize(std::pow(3*this->mNcv,2) + 5*this->mNcv);
      this->mRwork.resize(this->mNcv);

      this->mSelect.resize(this->mNcv);
      this->mD.resize(11);
      this->mZ.resize(this->mLhs.rows(),11);
      this->mWorkev.resize(2*this->mNcv);
   }

   void ArpackEigenSolver::solve(ArrayZ& eigenValues, MatrixZ&  eigenVectors)
   {
      // Safety assert for eigenValues and eigenVectors size
      assert(eigenValues.size() > 0);
      assert(eigenVectors.size() == 0 || eigenValues.size() == eigenVectors.cols());
      assert(eigenVectors.size() == 0 || eigenVectors.rows() == this->mLhs.rows());

      // Set Matrix size
      int n = this->mLhs.rows();

      // Extract the requested number of eigenvalues
      int nev = eigenValues.size();

      std::string bmat = "I";
      int lworkl = this->mWorkl.size();
      int ido = 0; 
      int info = 0;
      int ldv = n;

      // Set this->mIparam
      this->mIparam(0) = 1; // ishfts
      this->mIparam(2) = this->mMaxIter; // Maximum number of Arnoldi iterations
      this->mIparam(6) = 1; // mode

      // Loop until finished
      while(ido != 99)
      {
         // ARPACK's ZNAUPD
         znaupd_(&ido, bmat.c_str(), &n, this->mWhich.c_str(), &nev, &this->mTol, this->mResid.data(), &this->mNcv, this->mV.data(), &ldv, this->mIparam.data(), this->mIpntr.data(), this->mWorkd.data(), this->mWorkl.data(), &lworkl, this->mRwork.data(), &info);

         // ZNAUPD's IDO == 1
         if(ido == 1)
         {
            // WARNING: index from ARPACK are given starting from 1!
            this->mWorkd.segment(this->mIpntr(2)-1,n) = this->mWorkd.segment(this->mIpntr(0)-1, n);
            this->mWorkd.segment(this->mIpntr(1)-1,n) = this->mRhs*this->mWorkd.segment(this->mIpntr(2)-1, n);
            this->mWorkd.segment(this->mIpntr(1)-1,n) = this->mSolver.solve(this->mWorkd.segment(this->mIpntr(1)-1,n));

         // ZNAUPD's IDO == 99 
         } else if(ido == 99)
         {
         } else
         {
            throw Exception("ARPACK requested an unknown IDO flag!");
         }
      }

      if(info < 0)
      {
         throw Exception("ARPACK's ZNAUPD failed!");
      } else
      {
         int rvec = (eigenVectors.size() > 0);
         std::string howmny = "A";
         int ldz = n;

         if(nev > 10)
         {
            this->mD.resize(nev+1);
            this->mZ.resize(n, nev);
         }

         // ARPACK's ZNEUPD
         zneupd_(&rvec, howmny.c_str(), this->mSelect.data(), this->mD.data(), this->mZ.data(), &ldz, &this->mSigma, this->mWorkev.data(), bmat.c_str(), &n, this->mWhich.c_str(), &nev, &this->mTol, this->mResid.data(), &this->mNcv, this->mV.data(), &ldv, this->mIparam.data(), this->mIpntr.data(), this->mWorkd.data(), this->mWorkl.data(), &lworkl, this->mRwork.data(), &info);

         if(info < 0)
         {
            std::cerr << info << std::endl;
            throw Exception("ARPACK's ZNEUPD failed!");
         } else
         {
            eigenValues = (1/this->mD.array() + this->mSigma).segment(0,eigenValues.size());
            std::sort(eigenValues.data(), eigenValues.data()+eigenValues.size(), &sortEigenValues);
         }
      }
   }

   void ArpackEigenSolver::setSigma(const MHDComplex sigma)
   {
      this->mSigma = sigma;
   }

   void ArpackEigenSolver::setTol(const MHDFloat tol)
   {
      this->mTol = tol;
   }

   void ArpackEigenSolver::setMaxIter(const int maxIter)
   {
      this->mMaxIter = maxIter;
   }

   void ArpackEigenSolver::setNcv(const int ncv)
   {
      this->mNcv = ncv;
   }

   void ArpackEigenSolver::setWhich(const std::string& which)
   {
      this->mWhich = which;
   }

}
