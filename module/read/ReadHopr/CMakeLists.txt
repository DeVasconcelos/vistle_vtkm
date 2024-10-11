vistle_find_package(HDF5 COMPONENTS C)

add_module(ReadHopr "read HOPR HDF5 data" utils/ReadHDF5.h ReadHopr.h ReadHopr.cpp)

include_directories(${HDF5_INCLUDE_DIRS})

target_link_libraries(ReadHopr ${HDF5_LIBRARIES})
