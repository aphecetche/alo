cmake_minimum_required(VERSION 3.9)

cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0028 NEW)

find_package(AliRoot REQUIRED)
find_package(FlatBuffers REQUIRED)
find_package(ROOT REQUIRED)
find_package(Boost 1.59 COMPONENTS unit_test_framework filesystem program_options REQUIRED)

add_library(MCHMappingInterface INTERFACE)
target_include_directories(MCHMappingInterface INTERFACE 
$<BUILD_INTERFACE:${O2_MCH_MAPPING_SRC}/Interface/include>
$<INSTALL_INTERFACE:include>)

install(DIRECTORY ${O2_MCH_MAPPING_SRC}/Interface/include/MCHMappingInterface DESTINATION include)

set(SRCS
        ${O2_MCH_MAPPING_SRC}/Impl3/src/GenDetElemId2SegType.cxx
        ${O2_MCH_MAPPING_SRC}/Impl3/src/GenDetElemId2SegType.h
        ${O2_MCH_MAPPING_SRC}/Impl3/src/PadGroup.h
        ${O2_MCH_MAPPING_SRC}/Impl3/src/PadGroupType.cxx
        ${O2_MCH_MAPPING_SRC}/Impl3/src/PadGroupType.h
        ${O2_MCH_MAPPING_SRC}/Impl3/src/SegmentationCImpl3.cxx
        ${O2_MCH_MAPPING_SRC}/Impl3/src/SegmentationCreator.cxx
        ${O2_MCH_MAPPING_SRC}/Impl3/src/SegmentationCreator.h
        ${O2_MCH_MAPPING_SRC}/Impl3/src/SegmentationImpl3.cxx
        ${O2_MCH_MAPPING_SRC}/Impl3/src/SegmentationImpl3.h
        )

# We add all segmentation creators by default,
# but the final goal would be to tailor this
# for each executable reaching a given FLP (so it gets only
# the mapping it needs)
foreach (segtype RANGE 20)
    set(SRCS ${SRCS} ${O2_MCH_MAPPING_SRC}/Impl3/src/GenSegmentationCreatorForSegType${segtype}.cxx)
endforeach ()

add_library(MCHMappingImpl3 SHARED ${SRCS})
target_link_libraries(MCHMappingImpl3 PUBLIC MCHMappingInterface)
target_include_directories(MCHMappingImpl3 SYSTEM PRIVATE ${Boost_INCLUDE_DIR} 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> ${MS_GSL_INCLUDE_DIR})
install(TARGETS MCHMappingImpl3 LIBRARY DESTINATION lib)

include(GenerateExportHeader)
generate_export_header(MCHMappingImpl3)

set_target_properties(MCHMappingImpl3 PROPERTIES CXX_VISIBILITY_PRESET hidden)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(MCHMappingImpl3 PRIVATE -fext-numeric-literals)
endif ()

add_executable(mch-convert-esd-to-digits ${ALO_ESDCONVERTOR_SRC}/mch-convert-esd-to-digits.cxx ${ALO_ESDCONVERTOR_SRC}/FileCreation.cxx)
target_link_libraries(mch-convert-esd-to-digits PRIVATE flatbuffers::flatbuffers ROOT::RIO AliRoot::ESD MCHMappingImpl3 Boost::program_options)
install(TARGETS mch-convert-esd-to-digits RUNTIME DESTINATION bin)

add_executable(mch-dump-digits ${ALO_ESDCONVERTOR_SRC}/mch-dump-digits.cxx)
target_link_libraries(mch-dump-digits PRIVATE flatbuffers::flatbuffers MCHMappingImpl3)
install(TARGETS mch-dump-digits RUNTIME DESTINATION bin)

