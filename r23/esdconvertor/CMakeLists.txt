# FIXME: should call flatc to re-generate the header files
# from the .fbs files instead of using the ones commited ?

add_executable(mch-convert-esd 
        mch-convert-esd.cxx 
        FileCreation.cxx
        ConvertESD.cxx
        GetESDClusters.cxx
        SegmentationPair.cxx
        )

target_link_libraries(mch-convert-esd PRIVATE flatbuffers::flatbuffers ROOT::RIO AliRoot::ESD o2::MCHMappingSegContour3 o2::MCHMappingImpl3 Boost::program_options)
install(TARGETS mch-convert-esd RUNTIME DESTINATION bin)

add_executable(mch-dump mch-dump.cxx SegmentationPair.cxx DumpDigits.cxx DumpClusters.cxx)
target_link_libraries(mch-dump PRIVATE flatbuffers::flatbuffers o2::MCHMappingImpl3 Boost::program_options)
install(TARGETS mch-dump RUNTIME DESTINATION bin)

add_executable(test_FileCreation test_FileCreation.cxx FileCreation.cxx)
target_link_libraries(test_FileCreation PRIVATE o2::MCHMappingImpl3 Boost::unit_test_framework Boost::filesystem)
add_test(NAME test_FileCreation COMMAND test_FileCreation)
install(TARGETS test_FileCreation RUNTIME DESTINATION bin)
