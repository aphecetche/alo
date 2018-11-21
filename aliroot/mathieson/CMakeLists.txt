add_executable(test_Mathieson test_Mathieson.cxx)
target_link_libraries(test_Mathieson
                      PUBLIC AliRoot::MUONbase Boost::unit_test_framework)
add_test(NAME test_Mathieson COMMAND test_Mathieson)
