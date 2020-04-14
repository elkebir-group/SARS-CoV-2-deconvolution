cmake_minimum_required( VERSION 2.8 )

project( SARS-COV-2-DECONVOLUTION )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

find_package( GUROBI REQUIRED )

set (simulate_src
  src/simulatemain.cpp
  src/inputinstance.cpp
  src/utils.cpp
  src/solution.cpp
)

set (simulate_hdr
  src/inputinstance.h
  src/utils.h
  src/solution.h
)

set (gradient_src
  src/gradientmain.cpp
  src/solvergradient.cpp
  src/solvergradientb.cpp
  src/solvergradientu.cpp
  src/solver.cpp
  src/inputinstance.cpp
  src/utils.cpp
  src/solution.cpp
)

set (gradient_hdr
  src/solvergradientb.h
  src/solvergradient.h
  src/solvergradientu.h
  src/solver.h
  src/inputinstance.h
  src/utils.h
  src/solution.h
)

set (merge_src
  src/mergemain.cpp
  src/inputinstance.cpp
  src/utils.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set (merge_hdr
  src/inputinstance.h
  src/utils.h
  src/solution.h
  src/solutionset.h
)

set ( relax_src
  src/relaxmain.cpp
  src/inputinstance.cpp
  src/utils.cpp
  src/solver.cpp
  src/solverrelax.cpp
  src/solvercau.cpp
  src/solvercaul1.cpp
  src/solvercaul2.cpp
  src/solution.cpp
)

set ( relax_hdr
  src/inputinstance.h
  src/utils.h
  src/solver.h
  src/solverrelax.h
  src/solvercau.h
  src/solvercaul1.h
  src/solvercaul2.h
  src/solution.h
)

set( deconvolve_src
  src/deconvolvemain.cpp
  src/inputinstance.cpp
  src/utils.cpp
  src/dc.cpp
  src/solver.cpp
  src/solverca.cpp
  src/solvercabl1.cpp
  src/solvercab.cpp
  src/solvercaul1.cpp
  src/solvercau.cpp
  src/solvermilp.cpp
  src/solvermilpl1.cpp
  src/solvermilpbinom.cpp
  src/solution.cpp
)

set( deconvolve_hdr
  src/dkm/dkm.hpp
  src/dkm/dkm_utils.hpp
  src/inputinstance.h
  src/utils.h
  src/dc.h
  src/solver.h
  src/solverca.h
  src/solvercabl1.h
  src/solvercab.h
  src/solvercaul1.h
  src/solvercau.h
  src/solvermilp.h
  src/solvermilpl1.h
  src/solvermilpbinom.h
  src/solution.h
)

# first we can indicate the documentation build as an option and set it to OFF by default
option(BUILD_DOC "Build documentation" OFF)

if (BUILD_DOC)
  if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
    # check if Doxygen is installed
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        # set input and output files
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        # request to configure the file
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        message("Doxygen build started")

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
      message("Doxygen need to be installed to generate the doxygen documentation")
    endif (DOXYGEN_FOUND)
  endif()
endif()

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package( Boost REQUIRED COMPONENTS system filesystem program_options )
MESSAGE( STATUS "found boost:"
  "\nINCLUDE: ${Boost_INCLUDE_DIRS}"
  "\nLIB: ${Boost_LIBRARIES}"
)

include_directories( "src" ${Boost_INCLUDE_DIRS} "src/LBFGSpp" "src/Eigen" )
link_directories( ${Boost_LIBRARY_DIRS} ${GUROBI_LIBRARIES} )

set( CommonLibs
  m
  pthread
  ${Boost_LIBRARIES}
)

set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable( deconvolve EXCLUDE_FROM_ALL ${deconvolve_src} ${deconvolve_hdr} )
target_include_directories( deconvolve PUBLIC "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( deconvolve ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( merge ${merge_src} ${merge_hdr} )
target_include_directories( merge PUBLIC "src" ${Boost_INCLUDE_DIRS} )
target_link_libraries( merge ${CommonLibs} )

add_executable( relax EXCLUDE_FROM_ALL ${relax_src} ${relax_hdr} )
target_include_directories( relax PUBLIC "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( relax ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( gradient ${gradient_src} ${gradient_hdr} )
target_include_directories( gradient PUBLIC "src" ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} )
target_link_libraries( gradient ${CommonLibs} ${GUROBI_LIBRARIES} )

add_executable( simulate ${simulate_src} ${simulate_hdr} )
target_include_directories( simulate PUBLIC "src" ${Boost_INCLUDE_DIRS} )
target_link_libraries( simulate ${CommonLibs} )
