# This file is NOT licensed under the GPLv3, which is the license for the rest
# of YouCompleteMe.
#
# Here's the license text for this file:
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org/>

from distutils.sysconfig import get_python_inc
import os
import platform
import os.path as p
import subprocess

import logging
log = logging.getLogger('my_logger')

basedir = p.abspath( p.dirname( __file__ ) )
dbdir = basedir+'/.ycm_data'
dbdir_fallback = dbdir+'/fallback'
SOURCE_EXTENSIONS = [ '.cpp', '.cu' ]
include_basedir = 'include/QuICC'

database = None

# These are the compilation flags that will be used in case there's no
# compilation database set (by default, one is not set).
# CHANGE THIS LIST OF FLAGS. YES, THIS IS THE DROID YOU HAVE BEEN LOOKING FOR.
flags = [
  '-x', 
  'c++', 
  '-Wall', 
  '-Wextra', 
  '-Werror', 
  '-std=c++17'
]

def IsHeaderFile( filename ):
  extension = p.splitext( filename )[ 1 ]
  return extension in [ '.hpp', '.h' ]


def FindCorrespondingSourceFile( filename ):
  if IsHeaderFile( filename ):
    basename = p.splitext( filename )[ 0 ]
    src_basename, comp_basename = basename.split(include_basedir,1)
    src_basename += 'src/'
    for extension in SOURCE_EXTENSIONS:
      replacement_file = src_basename + comp_basename.split('/',1)[1] + extension
      if p.exists( replacement_file ):
        return replacement_file
      replacement_file = src_basename + comp_basename + extension
      if p.exists( replacement_file ):
        return replacement_file
  return filename


def Settings( **kwargs ):
  # Do NOT import ycm_core at module scope.
  import ycm_core

  global database
  if database is None and p.exists( dbdir ) and p.exists( dbdir + '/compile_commands.json'):
    log.info('DEBUGGING: trying local DB')
    database = ycm_core.CompilationDatabase( dbdir )
  elif database is None and p.exists( dbdir_fallback ):
    log.info('DEBUGGING: trying fallback DB')
    database = ycm_core.CompilationDatabase( dbdir_fallback )
    
  language = kwargs[ 'language' ]

  if language == 'cfamily':
    # If the file is a header, try to find the corresponding source file and
    # retrieve its flags from the compilation database if using one. This is
    # necessary since compilation databases don't have entries for header files.
    # In addition, use this source file as the translation unit. This makes it
    # possible to jump from a declaration in the header file to its definition
    # in the corresponding source file.
    filename = FindCorrespondingSourceFile( kwargs[ 'filename' ] )

    if not database:
      return {
        'flags': flags,
        'include_paths_relative_to_dir': basedir,
        'override_filename': filename
      }

    compilation_info = database.GetCompilationInfoForFile( filename )
    if not compilation_info.compiler_flags_:
      return {}

    # Bear in mind that compilation_info.compiler_flags_ does NOT return a
    # python list, but a "list-like" StringVec object.
    final_flags = list( compilation_info.compiler_flags_ )

    return {
      'flags': final_flags,
      'include_paths_relative_to_dir': compilation_info.compiler_working_dir_,
      'override_filename': filename
    }

  return {}
