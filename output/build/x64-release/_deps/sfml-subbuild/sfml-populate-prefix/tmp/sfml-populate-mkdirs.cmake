# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-src"
  "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-build"
  "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-subbuild/sfml-populate-prefix"
  "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-subbuild/sfml-populate-prefix/tmp"
  "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-subbuild/sfml-populate-prefix/src/sfml-populate-stamp"
  "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-subbuild/sfml-populate-prefix/src"
  "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-subbuild/sfml-populate-prefix/src/sfml-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "C:/Users/mur-m/OneDrive/Documents/Courseworks/Coursework V/Numerical Methods Coursework/output/build/x64-release/_deps/sfml-subbuild/sfml-populate-prefix/src/sfml-populate-stamp/${subDir}")
endforeach()
