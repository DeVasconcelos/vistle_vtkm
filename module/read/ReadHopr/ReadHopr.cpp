#include <hdf5.h>
#include <string>
#include <vector>

#include <vistle/core/unstr.h>

#include "ReadHopr.h"

MODULE_MAIN(ReadHopr)
using vistle::Parameter;
using vistle::Reader;

ReadHopr::ReadHopr(const std::string &name, int moduleID, mpi::communicator comm): Reader(name, moduleID, comm)
{
    m_meshFile = addStringParameter("mesh_file", "HOPR HDF5 (.h5) file containing the mesh information", "",
                                    Parameter::ExistingFilename);
    m_stateFile = addStringParameter("state_file", "HOPR HDF5 (.h5) file containing the state information", "",
                                     Parameter::ExistingFilename);

    auto hoprFormat = "HOPR HDF5 Files (*.h5)";
    setParameterFilters(m_meshFile, hoprFormat);
    setParameterFilters(m_stateFile, hoprFormat);

    observeParameter(m_meshFile);
    observeParameter(m_stateFile);
}

ReadHopr::~ReadHopr()
{}

bool ReadHopr::examine(const vistle::Parameter *param)
{
    return true;
}


vistle::Byte hoprToVistleType(int hoprType)
{
    //TODO: for now we only support linear HOPR grids
    switch (hoprType % 10) {
    case 3:
        return vistle::UnstructuredGrid::TRIANGLE;
    case 4:
        return hoprType < 100 ? vistle::UnstructuredGrid::QUAD : vistle::UnstructuredGrid::TETRAHEDRON;
    case 5:
        return vistle::UnstructuredGrid::PYRAMID;
    case 6:
        return vistle::UnstructuredGrid::PRISM;
    case 8:
        return vistle::UnstructuredGrid::HEXAHEDRON;
    default:
        throw vistle::exception("Encountered unsupported HOPR data type");
    }
}

bool ReadHopr::read(vistle::Reader::Token &token, int timestep, int block)
{
    // ---- READ IN MESH FILE ----
    auto h5Mesh = H5Fopen(m_meshFile->getValue().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // contains element type + list
    auto elemInfo_DId = H5Dopen(h5Mesh, "ElemInfo", H5P_DEFAULT);
    auto elemInfo_SId = H5Dget_space(elemInfo_DId);
    std::vector<hsize_t> elemInfoDim(H5Sget_simple_extent_ndims(elemInfo_SId));
    H5Sget_simple_extent_dims(elemInfo_SId, elemInfoDim.data(), nullptr);

    hsize_t total_size = 1;
    for (hsize_t dim: elemInfoDim) {
        total_size *= dim;
    }

    auto elemInfo_TId = H5Dget_type(elemInfo_DId);
    if (H5Tequal(elemInfo_TId, H5T_NATIVE_INT)) {
        // FIXME: make sure this vector is of the same datatype as used in the .h5 file
        std::vector<int> elemInfo(total_size);
        H5Dread(elemInfo_DId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, elemInfo.data());

        std::vector<int> type_list(elemInfoDim[0]);
        std::vector<int> element_list(elemInfoDim[0]);

        size_t counter = 0;
        for (hsize_t i = 0; i < total_size; i += 6) {
            type_list[counter] = hoprToVistleType(elemInfo[i]);
            element_list[counter] = elemInfo[i + 4];
            counter++;
        }
    } else {
        sendError("ElemInfo data type is not supported!");
    }

    H5Tclose(elemInfo_TId);
    H5Sclose(elemInfo_SId);
    H5Dclose(elemInfo_DId);

    // contains coordinates stored per element, i.e., connectivity list is implied
    auto nodeCoords_DId = H5Dopen(h5Mesh, "NodeCoords", H5P_DEFAULT);
    auto nodeCoords_SId = H5Dget_space(nodeCoords_DId);
    std::vector<hsize_t> nodeCoordsDim(H5Sget_simple_extent_ndims(nodeCoords_SId));
    H5Sget_simple_extent_dims(nodeCoords_SId, nodeCoordsDim.data(), nullptr);

    total_size = 1;
    for (hsize_t dim: nodeCoordsDim) {
        total_size *= dim;
    }

    auto nodeCoords_TId = H5Dget_type(nodeCoords_DId);
    if (H5Tequal(nodeCoords_TId, H5T_NATIVE_DOUBLE)) {
        // FIXME: make sure this vector is of the same datatype as used in the .h5 file
        std::vector<double> nodeCoords(total_size);
        H5Dread(nodeCoords_DId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nodeCoords.data());

        std::vector<double> x(nodeCoordsDim[0]);
        std::vector<double> y(nodeCoordsDim[0]);
        std::vector<double> z(nodeCoordsDim[0]);

        size_t counter = 0;
        for (hsize_t i = 0; i < total_size; i += 3) {
            x[counter] = nodeCoords[i];
            y[counter] = nodeCoords[i + 1];
            z[counter] = nodeCoords[i + 2];
            counter++;
        }

        for (hsize_t i = 0; i < nodeCoordsDim[0]; i++) {
            std::cout << x[i] << " " << y[i] << " " << z[i] << std::endl;
        }
    } else {
        sendError("NodeCoords datatype is not supported!");
    }

    H5Tclose(nodeCoords_TId);
    H5Sclose(nodeCoords_SId);
    H5Dclose(nodeCoords_DId);

    H5Fclose(h5Mesh);

    // ---- READ IN STATE FILE ----
    auto h5State = H5Fopen(m_stateFile->getValue().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    H5Fclose(h5State);

    return true;
}
bool ReadHopr::prepareRead()
{
    return true;
}
bool ReadHopr::finishRead()
{
    return true;
}
