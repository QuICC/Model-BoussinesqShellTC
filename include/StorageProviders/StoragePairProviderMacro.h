/** 
 * @file StoragePairProviderMacro.h
 * @brief Preprocessor macros used to setup the storage pair provider.
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

   /// Macro to select the type of storage pair provider
   #define StoragePairProviderMacro GeoMHDiSCC::DynamicPairProvider

#endif // STORAGEPAIRPROVIDERMACRO_H
