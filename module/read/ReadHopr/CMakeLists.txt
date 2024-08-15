vistle_find_package(HDF5)

add_module(ReadHopr "read HOPR HDF5 data" ReadHopr.h ReadHopr.cpp)

target_link_libraries(ReadHopr ${HDF5_LIBRARIES})
