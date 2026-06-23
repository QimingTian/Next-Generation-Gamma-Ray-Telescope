# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(GammaRayTelescope LANGUAGES CXX)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
option(WITH_HADRONIC "Enable FTFP_BERT hadronic physics for background studies" OFF)

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

if(WITH_HADRONIC)
  add_compile_definitions(WITH_HADRONIC)
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(Main Main.cc ${sources} ${headers})
target_link_libraries(Main ${Geant4_LIBRARIES})
target_compile_features(Main PRIVATE cxx_std_17)

if(WITH_HADRONIC)
  add_executable(MainHad Main.cc ${sources} ${headers})
  target_link_libraries(MainHad ${Geant4_LIBRARIES})
  target_compile_features(MainHad PRIVATE cxx_std_17)
endif()

file(GLOB MACRO_FILES "${PROJECT_SOURCE_DIR}/macros/*.mac")
foreach(_macro ${MACRO_FILES})
  get_filename_component(_name ${_macro} NAME)
  configure_file(${_macro} ${PROJECT_BINARY_DIR}/macros/${_name} COPYONLY)
endforeach()

install(TARGETS Main DESTINATION bin)
install(DIRECTORY macros/ DESTINATION share/GammaRayTelescope/macros)
