/** 
 * @file IConfigurationPart.cpp
 * @brief Source of the implementation of the base of the configuration file parts
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoConfig/ConfigParts/IConfigurationPart.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace QuICC {

namespace IoConfig {

   IConfigurationPart::IConfigurationPart(const std::string& parent)
      : mParent(parent)
   {
   }

   IConfigurationPart::~IConfigurationPart()
   {
   }

   const std::string& IConfigurationPart::parent() const
   {
      return this->mParent;
   }

   int IConfigurationPart::iValue(std::string name) const
   {
      // Make sure initialisation was correct
      assert(this->mIData.find(name) != this->mIData.end());

      return this->mIData.find(name)->second;
   }

   void IConfigurationPart::setIValue(std::string name, int value)
   {
      // Make sure initialisation was correct
      assert(this->mIData.find(name) != this->mIData.end());

      this->mIData.find(name)->second = value;
   }

   const std::map<std::string, int>& IConfigurationPart::iMap() const
   {
      return this->mIData;
   }

   std::map<std::string, int>& IConfigurationPart::rIMap()
   {
      return this->mIData;
   }

   std::pair<std::map<std::string, int>::const_iterator,std::map<std::string, int>::const_iterator> IConfigurationPart::iRange() const
   {
      return std::make_pair(this->mIData.begin(),this->mIData.end());
   }

   MHDFloat IConfigurationPart::fValue(std::string name) const
   {
      // Make sure initialisation was correct
      assert(this->mFData.find(name) != this->mFData.end());

      return this->mFData.find(name)->second;
   }

   void IConfigurationPart::setFValue(std::string name, MHDFloat value)
   {
      // Make sure initialisation was correct
      assert(this->mFData.find(name) != this->mFData.end());

      this->mFData.find(name)->second = value;
   }

   const std::map<std::string, MHDFloat>& IConfigurationPart::fMap() const
   {
      return this->mFData;
   }

   std::map<std::string, MHDFloat>& IConfigurationPart::rFMap()
   {
      return this->mFData;
   }

   std::pair<std::map<std::string, MHDFloat>::const_iterator,std::map<std::string, MHDFloat>::const_iterator> IConfigurationPart::fRange() const
   {
      return std::make_pair(this->mFData.begin(),this->mFData.end());
   }

   void IConfigurationPart::printInfo() const
   {
      // Create header
      IoTools::Formatter::printLine(std::cout, '-');
      IoTools::Formatter::printCentered(std::cout, this->parent(), '*');
      IoTools::Formatter::printLine(std::cout, '-');

      std::stringstream oss;
      std::map<std::string, int>::const_iterator  itI;
      for(itI = this->mIData.begin(); itI != this->mIData.end(); itI++)
      {
         oss << itI->first << ": " << itI->second;
         IoTools::Formatter::printCentered(std::cout, oss.str(), ' ');
         oss.str("");
      }

      std::map<std::string, MHDFloat>::const_iterator  itF;
      for(itF = this->mFData.begin(); itF != this->mFData.end(); itF++)
      {
         oss << itF->first << ": " << itF->second;
         IoTools::Formatter::printCentered(std::cout, oss.str(), ' ');
         oss.str("");
      }
      IoTools::Formatter::printLine(std::cout, '-');
      IoTools::Formatter::printNewline(std::cout);
   }

   void IConfigurationPart::addIntegerTag(const std::string& name, int value)
   {
      this->mIData.insert(std::make_pair(name, value));
   }

   void IConfigurationPart::addFloatTag(const std::string& name, MHDFloat value)
   {
      this->mFData.insert(std::make_pair(name, value));
   }

}
}
