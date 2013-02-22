#
# Create a configuration selection entry
#
function (geomhdiscc_create_choice choiceList choiceName chosen)
   foreach(choice ${${choiceList}})
      if(DEFINED str)
         set(str "${str}, ${choice}")
      else(DEFINED str)
         set(str "${choice}")
      endif(DEFINED str)
   endforeach(choice  ${${choiceList}})
   set(str "Available ${choiceName}(s): ${str}")
   set(${chosen} ${${chosen}} CACHE STRING ${str} FORCE)
endfunction(geomhdiscc_create_choice)

# 
# Verify that selected setup is valid
#
function (geomhdiscc_check_choice choiceList choiceName chosen out)
   set(${out} OFF PARENT_SCOPE)
   if(NOT ${chosen} STREQUAL "")
      list(FIND ${choiceList} ${${chosen}} test)
      if(test EQUAL -1)
         message(SEND_ERROR "------------>>> Unkown ${${choiceName}} <<<------------")
      else(test EQUAL -1)
         set(${out} ON PARENT_SCOPE)
         message(STATUS " --> ${choiceName}: "${${chosen}})
      endif(test EQUAL -1)
   else(NOT ${chosen} STREQUAL "")
      message(SEND_ERROR "------------>>> ${choiceName} is required <<<------------")
   endif(NOT ${chosen} STREQUAL "")
endfunction(geomhdiscc_check_choice)

#
# Provide a configuration selection entry
#
function (geomhdiscc_provide_choice choiceList choiceName chosen out)
   # Create selection entry
   geomhdiscc_create_choice(${choice} ${choiceList} ${choiceName} ${chosen})
   # Validate selection
   geomhdiscc_check_choice(${choice} ${choiceList} ${choiceName} ${chosen} ${out})
   set(${out} ${${out}} PARENT_SCOPE)
endfunction(geomhdiscc_provide_choice)

#
# Load platform file and setup requested platform
#
function (geomhdiscc_load_platform platform)
   include("cmake.d/platforms/${${platform}}.cmake")
   message(STATUS " --> Platform: "${${platform}})
endfunction (geomhdiscc_load_platform)

#
# Add definition to flags
#
function (geomhdiscc_add_definition base)
   string(TOUPPER "${base}_${${base}}" def)
   add_definitions("-D${def}")
endfunction (geomhdiscc_add_definition base)

#
# Recursive crawler through sources
#
function (geomhdiscc_append_sources MHDAll MHDPath MHDDirs)
   # Loop over all the subdirectories
   foreach(MHDDir ${${MHDDirs}})
      # Set new path including subdirectory
      set(MHDNext ${MHDPath}/${MHDDir})
      # Unset list of subdirectories and sources
      set(MHDSources)
      set(MHDSrcSubDirs)
      # Include SourcesList.cmake file
      include(${MHDNext}/SourcesList.cmake)
      # Check if there are additional subdirectories
      if(DEFINED MHDSrcSubDirs)
         geomhdiscc_append_sources(${MHDAll} ${MHDNext} MHDSrcSubDirs)
      endif()
      # Loop over all sources and add full path source to list
      foreach(MHDSource ${MHDSources})
         list(APPEND ${MHDAll} ${MHDNext}/${MHDSource})
      endforeach(MHDSource)
      # Update the list of all sources in parent scope
      set(${MHDAll} ${${MHDAll}} PARENT_SCOPE)
   endforeach(MHDDir)
endfunction ()
