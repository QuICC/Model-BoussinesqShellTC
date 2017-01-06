/** 
 * @file StaticAssert.hpp
 * @brief Simple static assert implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STATICASSERT_HPP
#define STATICASSERT_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

namespace Debug {

   /**
    * @brief Simple static type assert implementation
    */
   template <typename T1, typename T2> class StaticTypeAssert;

   /**
    * @brief Simple static type assert implementation (OK if same type)
    */
   template <typename T> class StaticTypeAssert<T,T>{};

   /**
    * @brief Simple static assert implementation
    */
   template <bool> class StaticAssert;

   /**
    * @brief Simple static assert implementation (ok if true)
    */
   template <> class StaticAssert<true> {};
}
}

#endif // STATICASSERT_HPP
