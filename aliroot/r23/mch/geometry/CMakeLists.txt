add_executable(mch-export-geometry mch-export-geometry.cxx)
target_link_libraries(mch-export-geometry PRIVATE Boost::program_options ROOT::Geom AliRoot::CDB AliRoot::STEER RapidJSON::RapidJSON)
