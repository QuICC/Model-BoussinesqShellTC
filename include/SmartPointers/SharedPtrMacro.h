/** 
 * @file SharedPtrMacro.h
 * @brief Preprocessor macros used to setup the shared pointer implementation depending on CMake setup.
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHAREDPTRMACRO_H
#define SHAREDPTRMACRO_H

// Boost Version
#ifdef QUICC_SMARTPTR_BOOST
   // Include the right header
   #include <boost/shared_ptr.hpp>

   /**
    * @def SharedPtrMacro
    * Macro allowing to use different implementations of the shared_ptr smart pointer.
    * Sets up the Boost version.
    */
   #define SharedPtrMacro boost::shared_ptr
#endif //QUICC_SMARTPTR_BOOST

// std (at least C++11) Version
#ifdef QUICC_SMARTPTR_STD
   // Include the right header
   #include <memory>

   /**
    * @def SharedPtrMacro
    * Macro allowing to use different implementations of the shared_ptr smart pointer.
    */
   #define SharedPtrMacro std::shared_ptr
#endif //QUICC_SMARTPTR_STD

#endif // SHAREDPTRMACRO_H
