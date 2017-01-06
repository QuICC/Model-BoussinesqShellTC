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
//   #define StoragePairProviderMacro QuICC::FixedPairProvider

   // Include the right header
   #include <StorageProviders/DynamicPairProvider.hpp>

   /// Macro to select the type of storage pair provider
   #define StoragePairProviderMacro QuICC::DynamicPairProvider

#endif // STORAGEPAIRPROVIDERMACRO_H
