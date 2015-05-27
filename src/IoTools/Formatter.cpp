/** 
 * @file Formatter.cpp
 * @brief Source of the implementation of a few formatting functions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "IoTools/Formatter.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoTools {

   void Formatter::printNewline(std::ostream& stream)
   {
      stream << std::endl;
   }

   void Formatter::printLine(std::ostream& stream, char fill)
   {
      stream << std::setfill(fill) << std::setw(Formatter::TEXT_WIDTH) << "" << std::endl;
   }

   void Formatter::printHalfLine(std::ostream& stream, char fillLeft, char fillRight)
   {
      stream << std::setfill(fillLeft) << std::setw(Formatter::TEXT_WIDTH/2) << "" << std::setfill(fillRight) << std::setw(Formatter::TEXT_WIDTH/2) << "" << std::endl;
   }

   void Formatter::printCentered(std::ostream& stream, const std::string& text, char fill, int base)
   {
      int pad = Formatter::TEXT_WIDTH - std::max(static_cast<int>(text.length()),base) - 2;

      stream << std::setfill(fill) << std::setw(text.length() + pad/2 + 2) << " " + text + " " << std::setw(pad/2 + pad%2) << "" << std::endl;
   }

   void Formatter::printLeftCentered(std::ostream& stream, const std::string& text, char fill, int base, bool hasNewline)
   {
      int pad = Formatter::TEXT_WIDTH - std::max(static_cast<int>(text.length()),base) - 2;

      stream << std::setfill(fill) << std::setw(text.length() + pad/2 + 2) << " " + text;
      if(hasNewline)
      {
         Formatter::printNewline(stream);
      }
   }

   Formatter::Formatter()
   {
   }
}
}
