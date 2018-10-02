cmake_minimum_required(VERSION 3.6)
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0028 NEW)

find_package(AliRoot REQUIRED)
find_package(O2 REQUIRED)
find_package(FlatBuffers REQUIRED)
find_package(ROOT REQUIRED)
find_package(Boost 1.59 COMPONENTS unit_test_framework filesystem program_options REQUIRED)

add_executable(mch-convert-esd-to-digits mch-convert-esd-to-digits.cxx FileCreation.cxx)
target_link_libraries(mch-convert-esd-to-digits PRIVATE flatbuffers::flatbuffers ROOT::RIO AliRoot::ESD o2::MCHMappingSegContour3 o2::MCHMappingImpl3 Boost::program_options)
install(TARGETS mch-convert-esd-to-digits RUNTIME DESTINATION bin)

add_executable(mch-dump-digits mch-dump-digits.cxx)
target_link_libraries(mch-dump-digits PRIVATE flatbuffers::flatbuffers o2::MCHMappingImpl3)
install(TARGETS mch-dump-digits RUNTIME DESTINATION bin)

add_executable(test_FileCreation test_FileCreation.cxx FileCreation.cxx)
target_link_libraries(test_FileCreation PRIVATE o2::MCHMappingImpl3 Boost::unit_test_framework Boost::filesystem)
add_test(NAME test_FileCreation COMMAND test_FileCreation)
install(TARGETS test_FileCreation RUNTIME DESTINATION bin)
