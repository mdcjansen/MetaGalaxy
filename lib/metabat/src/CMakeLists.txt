# find and use Boost
set(Boost_NO_BOOST_CMAKE ON) # for backwards comptibility of cmake and boost >= 1.70
set(Boost_USE_STATIC_LIBS   ON)
find_package(Boost 1.55.0 COMPONENTS program_options filesystem system graph serialization iostreams REQUIRED)
include_directories( ${Boost_INCLUDE_DIRS} )

# find and use OpenMP
find_package(OpenMP)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    message("Using OpenMP")
else()
    message(WARNING "OpenMP must be available to your compiler or it will be REALLLLLLLLLLY slow: ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER}")
endif()

if(APPLE)
	#full Single Unix Standard v3 (SUSv3) conformance (the Unix API)
	add_definitions(-D_DARWIN_C_SOURCE)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -qreport -qhot -qmaxmem=819200 -qoffload -qtgtarch=sm_70 -g")
set(EXTRALIBS ${EXTRALIBS} )
endif()

message("Boost cmake dir: ${Boost_DIR}")
message("Boost include: ${Boost_INCLUDE_DIRS}")
message("Boost libraries:" ${Boost_LIBRARIES})
message("Boost lib ver:" ${Boost_LIB_VERSION})
message("Boost lib dir:" ${Boost_LIBRARY_DIRS})

set(targets metabat1 metabat2)
foreach(target ${targets})

   add_executable(${target} ${target}.cpp)
   add_dependencies(${target} htslib zlib check_git_repository)
   target_link_libraries(${target} ${zlib_LIB} ${Boost_LIBRARIES} ${EXTRALIBS} ) 

endforeach()

install(TARGETS ${targets}
        DESTINATION bin/
       )

set(targets jgi_summarize_bam_contig_depths contigOverlaps)
foreach(target ${targets})
   add_executable(${target} ${target}.cpp)
   add_dependencies(${target} htslib zlib check_git_repository)
   target_link_libraries(${target} ${htslib_LIB} ${zlib_LIB} ${Boost_LIBRARIES} ) 
endforeach()


install(TARGETS ${targets}
        DESTINATION bin/
       )
