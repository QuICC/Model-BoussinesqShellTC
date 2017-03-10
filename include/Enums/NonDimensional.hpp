/**
 * @file NonDimensional.hpp
 * @brief Definition of some useful enums for nondimensional paramters 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef NONDIMENSIONAL_HPP
#define NONDIMENSIONAL_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

      /**
       * @brief Simple struct holding the mapping for the nondimensional parameters
       */
      struct NonDimensional {
         /**
          * @brief Enums of the different nondimensional factors
          */
         enum Id {
            //
            // Nondimensional numbers
            /// Eady number
            EADY = 0,
            /// Ekman number
            EKMAN,
            /// Magnetic Ekman number
            MAGEKMAN,
            /// Magnetic Prandtl number
            MAGPRANDTL,
			   // Nicolò Lardelli
			   /// Modified Elsasser number
			   MODELSASSER,
            /// Prandtl number
            PRANDTL,
            /// Poincare number
            POINCARE,
            /// Rayleigh number
            RAYLEIGH,
            /// Roberts number
            ROBERTS,
            /// Rossby number
            ROSSBY,
            /// Taylor number
            TAYLOR,

            // 
            // Geometrical numbers
            /// Outer radius R_o
            RO,
            /// Radii ratio R_i/R_o
            RRATIO,

            //
            // Flags
            /// Flag to switch between interal heating, differential heating, etc
            HEATING,
            
            //
            // Axis scaling factors
            /// 1D axis scale
            SCALE1D,
            /// 2D axis scale
            SCALE2D,
            /// 3D axis scale
            SCALE3D,

            //
            // Greek alphabet
            /// Alpha
            ALPHA,
            /// Beta
            BETA,
            /// Gamma
            GAMMA,
            /// Delta
            DELTA,
            /// Epsilon
            EPSILON,
            /// Zeta
            ZETA,
            /// Eta
            ETA,
            /// Theta
            THETA,
            /// Iota
            IOTA,
            /// Kappa
            KAPPA,
            /// Lambda
            LAMBDA,
            /// Mu
            MU,
            /// Nu
            NU,
            /// Xi
            XI,
            /// Omicron
            OMICRON,
            /// Pi
            PI,
            /// Rho
            RHO,
            /// Sigma
            SIGMA,
            /// Tau
            TAU,
            /// Upsilon
            UPSILON,
            /// Phi
            PHI,
            /// Chi
            CHI,
            /// Psi
            PSI,
            /// Omega
            OMEGA,

            //
            // Special flags
            /// Flag to filter elevator modes
            ELEVATOR,
            /// Flag to use fast mean equation
            FAST_MEAN,
            /// Flag to use rescaled equation
            RESCALED,
         };
      };
}

#endif // NONDIMENSIONAL_HPP
