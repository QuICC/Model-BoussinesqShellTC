/** 
 * @file IConfigurationPart.hpp 
 * @brief Implementation of the truncation part of the configuration file
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ICONFIGURATIONPART_HPP
#define ICONFIGURATIONPART_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <string>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace IoConfig {

   /**
    * @brief Implementation of the truncation part of the configuration file
    */
   class IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          */
         explicit IConfigurationPart(const std::string& parent);

         /**
          * @brief Destructor
          */
         virtual ~IConfigurationPart();

         /**
          * @brief Get parent node XML tag
          */
         const std::string& parent() const;

         /**
          * @brief Get integer value by name
          *
          * @param name Name of the value
          */
         int iValue(std::string name) const;

         /**
          * @brief Set integer value by name
          *
          * @param name    Name of the value
          * @param value   New value
          */
         void setIValue(std::string name, int value);

         /**
          * @brief Get integer name to value map
          */
         const std::map<std::string, int>& iMap() const;

         /**
          * @brief Get integer name to value map
          */
         std::map<std::string, int>& rIMap();

         /**
          * @brief Get integer name to value map iterator range
          */
         std::pair<std::map<std::string, int>::const_iterator, std::map<std::string, int>::const_iterator> iRange() const;

         /**
          * @brief Get float value by name
          *
          * @param name Name of the value
          */
         MHDFloat fValue(std::string name) const;

         /**
          * @brief Set float value by name
          *
          * @param name    Name of the value
          * @param value   New value
          */
         void setFValue(std::string name, MHDFloat value);

         /**
          * @brief Get float name to value map
          */
         const std::map<std::string, MHDFloat>& fMap() const;

         /**
          * @brief Get float name to value map
          */
         std::map<std::string, MHDFloat>& rFMap();

         /**
          * @brief Get float name to value map iterator range
          */
         std::pair<std::map<std::string, MHDFloat>::const_iterator, std::map<std::string, MHDFloat>::const_iterator> fRange() const;

         /**
          * @brief Check compatibility of data
          */
         virtual void checkData() = 0;

         /**
          * @brief Output run information
          */
         virtual void printInfo() const;
         
      protected:
         /**
          * @brief Add XML tag with integer data
          */
         void addIntegerTag(const std::string& name, int value);

         /**
          * @brief Add XML tag with integer data
          */
         void addFloatTag(const std::string& name, MHDFloat value);

      private:
         /**
          * @brief XML parent tag name
          */
         std::string mParent;

         /**
          * @brief XML tags of the integer data
          */
         std::map<std::string, int> mIData;

         /**
          * @brief XML tags of the float data
          */
         std::map<std::string, MHDFloat> mFData;
   };

   /// Typedef for a shared pointer of a configuration file part
   typedef SharedPtrMacro<IConfigurationPart> SharedIConfigurationPart;

   /// Typedef for a shared pointer of a configuration file part
   typedef SharedPtrMacro<const IConfigurationPart> SharedCIConfigurationPart;

}
}

#endif // ICONFIGURATIONPART_HPP
