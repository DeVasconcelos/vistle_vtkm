#include <hdf5.h>
#include <string>
#include <vector>

#include <vistle/core/unstr.h>
#include <vistle/module/module.h>

#include "ReadHopr.h"

using namespace vistle;
MODULE_MAIN(ReadHopr)

// TODO: find out why VTK produces 8x more cells than reading in the .h5 mesh...

// FIXME: Only for single-process mode: When attempting to read in multiple .h5 files in the
//        same Vistle map, vistle will crash UNLESS hdf5 was compiled to be thread-safe (i.e., with
//        ./configure --enable-threadsafe --enable-unsupported).
//        --> Make sure we are using thread-safe hdf5 when compiling vistle in single-process mode
//        (probably through a hdf5 submodule.)

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


template<typename T>
void readH5Dataset(hid_t fileId, const char *datasetName, std::vector<T> &result)
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

Byte hoprToVistleType(int hoprType)
{
    //TODO: add support for bilinear and nonlinear HOPR cell types
    switch (hoprType) {
    case 3:
        return UnstructuredGrid::TRIANGLE;
    case 4:
        return UnstructuredGrid::QUAD;
    case 104:
        return UnstructuredGrid::TETRAHEDRON;
    case 105:
        return UnstructuredGrid::PYRAMID;
    case 106:
        return UnstructuredGrid::PRISM;
    case 108:
        return UnstructuredGrid::HEXAHEDRON;
    default:
        std::stringstream msg;
        msg << "The HOPR data type with the encoding " << hoprType << " is not supported.";

        std::cerr << msg.str() << std::endl;
        throw exception(msg.str());
    }
}

size_t addCellToConnectivityList(UnstructuredGrid::ptr grid, size_t offset, Byte cellType)
{
    // TODO: create tests (read in .h5 meshes including each cell type)
    std::vector<Byte> order;
    switch (cellType) {
    case UnstructuredGrid::TRIANGLE:
        order = {0, 1, 2};
        break;
    case UnstructuredGrid::QUAD:
    case UnstructuredGrid::TETRAHEDRON:
        order = {0, 1, 2, 3};
        break;
    case UnstructuredGrid::PYRAMID:
        order = {0, 1, 3, 2, 4};
        break;
    case UnstructuredGrid::PRISM:
        order = {0, 1, 2, 3, 4, 5};
        break;
    case UnstructuredGrid::HEXAHEDRON:
        order = {0, 1, 3, 2, 4, 5, 7, 6};
        break;
    default:
        throw exception("Found unhandled cell type after converting Hopr to vistle types.");
    }

    for (size_t i = 0; i < order.size(); i++) {
        grid->cl()[offset + i] = offset + order[i];
    }
    return order.size();
}

UnstructuredGrid::ptr ReadHopr::createMeshFromFile(const char *filename)
{
    // read in information necessary to build an unstructured grid
    auto h5Mesh = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (h5Mesh < 0) {
        sendError("An error occurred while reading in mesh file. Cannot create mesh!");
        return nullptr;
    }

    std::vector<int> elemInfo;
    readH5Dataset(h5Mesh, "ElemInfo", elemInfo);

    std::vector<double> nodeCoords;
    readH5Dataset(h5Mesh, "NodeCoords", nodeCoords);

    // create the unstructured grid
    UnstructuredGrid::ptr result(
        new UnstructuredGrid(elemInfo.size() / 6, nodeCoords.size() / 3, nodeCoords.size() / 3));

    // Hopr's 'NodeCoords' consists of three columns: the x, y and z coordinates of the points
    // that make up the elements of the grid. The coordinates are ordered by element/cell, i.e.,
    // the same point appears multiple times in the array if it belongs to more than one element.
    if (nodeCoords.size() > 0) {
        size_t counter = 0;
        for (hsize_t i = 0; i < nodeCoords.size(); i += 3) {
            result->x()[counter] = nodeCoords[i];
            result->y()[counter] = nodeCoords[i + 1];
            result->z()[counter] = nodeCoords[i + 2];
            counter++;
        }
    } else {
        sendError("An exception occurred while reading in 'NodeCoords'. Cannot create mesh.");
        H5Fclose(h5Mesh);
        return nullptr;
    }

    // Hopr's 'ElemInfo' consists of six columns: ...
    if (elemInfo.size() > 0) {
        Byte vistleType;

        size_t counter = 0;
        size_t clSize = 0;
        for (hsize_t i = 0; i < elemInfo.size(); i += 6) {
            // ... the 1st column contains the element types (= vistle's type list 'tl')
            try {
                vistleType = hoprToVistleType(elemInfo[i]);
            } catch (...) {
                sendError("Encountered unsupported HOPR data type. Please note that bilinear and non-linear cell types "
                          "are not supported yet.");
                H5Fclose(h5Mesh);
                return nullptr;
            }
            result->tl()[counter] = vistleType;

            // ... the 5th column contains the offsets into the point coordinates list (which, in this
            // case, corresponds to vistle's element list 'el' because the point coordinates are stored cell-wise)
            if (i > 0)
                result->el()[counter] = elemInfo[i + 4];

            // ... there is no equivalent for vistle's connectivity list. Since the node order is not always the
            // same as in vistle, we have to create it ourselves.
            clSize += addCellToConnectivityList(result, clSize, vistleType);
            counter++;
        }
        result->el()[result->el().size() - 1] = clSize;

    } else {
        sendError("An exception occurred while reading in 'ElemInfo'. Cannot create mesh.");
        H5Fclose(h5Mesh);
        return nullptr;
    }

    H5Fclose(h5Mesh);
    return result;
}

bool ReadHopr::read(Reader::Token &token, int timestep, int block)
{
    auto meshFileName = m_meshFile->getValue();

    UnstructuredGrid::ptr result;
    if (meshFileName.size()) {
        result = createMeshFromFile(meshFileName.c_str());
    } else {
        sendError("No mesh file is given, cannot create mesh.");
        return true;
    }

    // ---- READ IN STATE FILE ----
    /*auto h5State = H5Fopen(m_stateFile->getValue().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    std::vector<double> DGSolution;
    readH5Dataset(h5State, "DG_Solution", DGSolution);

    H5Fclose(h5State);
    */

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
