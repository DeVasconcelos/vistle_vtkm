#include <hdf5.h>
#include <string>
#include <vector>

#include <vistle/core/unstr.h>
#include <vistle/module/module.h>

#include "ReadHopr.h"

using namespace vistle;
MODULE_MAIN(ReadHopr)

ReadHopr::ReadHopr(const std::string &name, int moduleID, mpi::communicator comm): Reader(name, moduleID, comm)
{
    m_meshFile = addStringParameter("mesh_file", "HOPR HDF5 (.h5) file containing the mesh information", "",
                                    Parameter::ExistingFilename);
    m_stateFile = addStringParameter("state_file", "HOPR HDF5 (.h5) file containing the state information", "",
                                     Parameter::ExistingFilename);

    auto hoprFormat = "HOPR HDF5 Files (*.h5)";
    setParameterFilters(m_meshFile, hoprFormat);
    setParameterFilters(m_stateFile, hoprFormat);

    m_gridOut = createOutputPort("grid_out", "grid");

    observeParameter(m_meshFile);
    observeParameter(m_stateFile);
}

ReadHopr::~ReadHopr()
{}

bool ReadHopr::examine(const Parameter *param)
{
    return true;
}


Byte hoprToVistleType(int hoprType)
{
    //TODO: for now we only support linear HOPR grids
    switch (hoprType % 10) {
    case 3:
        return UnstructuredGrid::TRIANGLE;
    case 4:
        return hoprType < 100 ? UnstructuredGrid::QUAD : UnstructuredGrid::TETRAHEDRON;
    case 5:
        return UnstructuredGrid::PYRAMID;
    case 6:
        return UnstructuredGrid::PRISM;
    case 8:
        return UnstructuredGrid::HEXAHEDRON;
    default:
        throw exception("Encountered unsupported HOPR data type");
    }
}

template<typename T>
void readDataset(hid_t fileId, const char *datasetName, std::vector<T> &result)
{
    auto datasetId = H5Dopen(fileId, datasetName, H5P_DEFAULT);
    auto spaceId = H5Dget_space(datasetId);
    std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(spaceId));
    H5Sget_simple_extent_dims(spaceId, shape.data(), nullptr);

    hsize_t total_size = 1;
    for (hsize_t dim: shape) {
        total_size *= dim;
    }

    auto dtypeId = H5Dget_type(datasetId);
    try {
        result.resize(total_size);
        H5Dread(datasetId, dtypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, result.data());
    } catch (...) {
        std::cerr << "An exception occurred when trying to read in " << datasetName << " of type " << dtypeId << "."
                  << std::endl;
    }

    H5Tclose(dtypeId);
    H5Sclose(spaceId);
    H5Dclose(datasetId);
}

UnstructuredGrid::ptr ReadHopr::readMesh(const char *filename)
{
    // read in information necessary to build an unstructured grid
    auto h5Mesh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    std::vector<int> elemInfo;
    readDataset(h5Mesh, "ElemInfo", elemInfo);

    std::vector<double> nodeCoords;
    readDataset(h5Mesh, "NodeCoords", nodeCoords);

    // create the unstructured grid
    size_t numElements = elemInfo.size() / 6;
    size_t numCorners = nodeCoords.size() / 3;
    size_t numVertices = nodeCoords.size() / 3;

    UnstructuredGrid::ptr result(new UnstructuredGrid(numElements, numCorners, numVertices));

    // Hopr's 'ElemInfo' consists of six columns: the first and fifth column contain the element
    // types (= vistle's type list) and offsets into the connectivity list (= vistle element list),
    // respectively
    if (elemInfo.size() > 0) {
        auto counter = 0;
        for (hsize_t i = 0; i < elemInfo.size(); i += 6) {
            result->tl()[counter] = hoprToVistleType(elemInfo[i]);
            if (i > 0)
                result->el()[counter] = elemInfo[i + 4];
            counter++;
        }
        result->el()[result->el().size() - 1] = numCorners;

    } else {
        sendError("An exception occurred while reading in 'ElemInfo'. Cannot create mesh.");
        return nullptr;
    }

    // Hopr's 'NodeCoords' consists of three columns: the x, y and z coordinates of the points
    // that make up the elements of the grid. The coordinates are ordered by element, i.e.,
    // -> the same point appears multiple times in the array if it belongs to more than one element.
    // -> vistle's connectivity list is simple a range from 0 to #points in nodeCoords

    // BUG: the order used to define an element out of an array of coordinates seems to differ between
    // Hopr and vistle (at least for hexahedra)
    if (nodeCoords.size() > 0) {
        auto counter = 0;
        for (hsize_t i = 0; i < nodeCoords.size(); i += 3) {
            result->x()[counter] = nodeCoords[i];
            result->y()[counter] = nodeCoords[i + 1];
            result->z()[counter] = nodeCoords[i + 2];
            result->cl()[counter] = counter;
            counter++;
        }
    } else {
        sendError("An exception occurred while reading in 'NodeCoords'. Cannot create mesh.");
        return nullptr;
    }

    H5Fclose(h5Mesh);

    return result;
}

bool ReadHopr::read(Reader::Token &token, int timestep, int block)
{
    auto result = readMesh(m_meshFile->getValue().c_str());

    // ---- READ IN STATE FILE ----
    auto h5State = H5Fopen(m_stateFile->getValue().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    std::vector<double> DGSolution;
    readDataset(h5State, "DG_Solution", DGSolution);

    H5Fclose(h5State);

    updateMeta(result);
    addObject(m_gridOut, result);

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
