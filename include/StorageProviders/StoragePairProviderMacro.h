/** \file StoragePairProviderMacro.h
 *  \brief Preprocessor macros used to setup the shared pointer implementation depending on CMake setup.
 *
 *  \mhdBug Needs test
 */

#ifndef STORAGEPAIRPROVIDERMACRO_H
#define STORAGEPAIRPROVIDERMACRO_H

//   // Include the right header
//   #include <StorageProviders/FixedPairProvider.hpp>
//
//   // Fixed size storage pair provider
//   #define StoragePairProviderMacro GeoMHDiSCC::FixedPairProvider

   // Include the right header
   #include <StorageProviders/DynamicPairProvider.hpp>

   // Dynamic size storage pair provider
   #define StoragePairProviderMacro GeoMHDiSCC::DynamicPairProvider

#endif // STORAGEPAIRPROVIDERMACRO_H
