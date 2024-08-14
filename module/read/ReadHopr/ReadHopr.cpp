#include <hdf5.h>
#include <string>

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
bool ReadHopr::read(vistle::Reader::Token &token, int timestep, int block)
{
    // ---- READ IN MESH FILE ----
    auto h5Mesh = H5Fopen(m_meshFile->getValue().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // contains element type + list
    auto elemInfoId = H5Dopen(h5Mesh, "ElemInfo", H5P_DEFAULT);
    H5Dclose(elemInfoId);

    // contains coordinates stored per element, i.e., connectivity list is implied
    auto nodeCoordsId = H5Dopen(h5Mesh, "NodeCoords", H5P_DEFAULT);
    H5Dclose(nodeCoordsId);

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
