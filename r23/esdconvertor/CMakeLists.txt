
find_package(FlatBuffers REQUIRED)
find_package(ROOT REQUIRED)
find_package(AliRoot REQUIRED)
find_package(Boost COMPONENTS unit_test_framework filesystem REQUIRED)

add_executable(mch-esd-to-digits mch-esd-to-digits.cxx FileCreation.cxx)
target_link_libraries(mch-esd-to-digits PRIVATE flatbuffers::flatbuffers ROOT::RIO AliRoot::ESD AliRoot::MUONbase o2::MCHMappingSegContour3 o2::MCHMappingImpl3)
install(TARGETS mch-esd-to-digits RUNTIME DESTINATION bin)

add_executable(mch-digits-dump mch-digits-dump.cxx)
target_link_libraries(mch-digits-dump PRIVATE flatbuffers::flatbuffers o2::MCHMappingImpl3)
install(TARGETS mch-digits-dump RUNTIME DESTINATION bin)

add_executable(test_FileCreation test_FileCreation.cxx FileCreation.cxx)
target_link_libraries(test_FileCreation PRIVATE o2::MCHMappingImpl3 Boost::unit_test_framework Boost::filesystem)
add_test(NAME test_FileCreation COMMAND test_FileCreation)
install(TARGETS test_FileCreation RUNTIME DESTINATION bin)
