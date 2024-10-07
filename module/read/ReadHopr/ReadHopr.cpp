#include <hdf5.h>
#include <cctype>
#include <cmath>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <vistle/core/unstr.h>
#include <vistle/module/module.h>

#include "ReadHopr.h"

using namespace vistle;
MODULE_MAIN(ReadHopr)

const std::string Invalid("(NONE)");

// TODO: find out why VTK produces 8x more cells than reading in the .h5 mesh...
//       --> that's because the solution is a polynomial of degree 4 (can be read out of state file)
// TODO: create higher order elements! (right now we are only storing the corner nodes, i.e., pretending
//       the state is linear)

// While the C version of HDF5 can be compiled to be threadsafe (with './configure --enable-threadsafe
// --enable-unsupported'), it is not by default. This is an issue, when compiling vistle in single-process
// mode, as calling the ReadHopr-module multiple times at the same time leads to vistle crashing.
// To allow the user to use any HDF5 package, even if it is no threadsafe, we create a mutex and
// corresponding lock- and unlock-functions that can be used to make sure that the HDF5 library is not
// accessed by two threads at the same time.
#if defined(MODULE_THREAD) // If VISTLE_MULTI_PROCESS is OFF...
static std::mutex hdf5_mutex; // ...avoid simultaneous access to HDF5 library.
#ifdef COLLECTIVE
#define LOCK_HDF5(comm) \
    std::unique_lock<std::mutex> hdf5_guard(hdf5_mutex, std::defer_lock); \
    if ((comm).rank() == 0) \
        hdf5_guard.lock(); \
    (comm).barrier();
#define UNLOCK_HDF5(comm) \
    (comm).barrier(); \
    if (hdf5_guard) \
        hdf5_guard.unlock();
#else
#define LOCK_HDF5(comm) \
    std::unique_lock<std::mutex> hdf5_guard(hdf5_mutex, std::defer_lock); \
    hdf5_guard.lock();
#define UNLOCK_HDF5(comm) \
    if (hdf5_guard) \
        hdf5_guard.unlock();
#endif
#else
#define LOCK_HDF5(comm)
#define UNLOCK_HDF5(comm)
#endif

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

    for (int i = 0; i < NumPorts; i++) {
        std::stringstream choiceFieldName;
        choiceFieldName << "state_field_" << i;

        m_fieldChoice[i] = addStringParameter(
            "state_field_" + std::to_string(i),
            "This data field from the state file will be added to output port field_out_" + std::to_string(i) + ".", "",
            Parameter::Choice);
        m_fieldsOut[i] = createOutputPort("field_out_" + std::to_string(i), "data field");
    }

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
struct H5Dataset {
    std::vector<T> vector;
    std::vector<hsize_t> dimension;
};

//TODO: strings must be handled differently (as they can be variable or fixed size, see 'readH5Attribute<std::string>')
template<typename T>
H5Dataset<T> readH5Dataset(hid_t fileId, const char *datasetName)
{
    H5Dataset<T> result;
    auto datasetId = H5Dopen(fileId, datasetName, H5P_DEFAULT);
    if (datasetId < 0) {
        std::cerr << "Could not open dataset " << datasetName << "!" << std::endl;
        H5Dclose(datasetId);
        return result;
    }
    auto spaceId = H5Dget_space(datasetId);
    std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(spaceId));
    H5Sget_simple_extent_dims(spaceId, shape.data(), nullptr);

    for (hsize_t dim: shape) {
        result.dimension.push_back(dim);
    }

    auto totalSize = std::accumulate(result.dimension.begin(), result.dimension.end(), 1, std::multiplies<hsize_t>());

    auto dtypeId = H5Dget_type(datasetId);
    result.vector.resize(totalSize);
    if (H5Dread(datasetId, dtypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT, result.vector.data()) < 0) {
        std::cerr << "An error occurred when trying to read in " << datasetName << " of type " << dtypeId << "."
                  << std::endl;
        result.vector.clear();
    }

    //TODO: find a way to not always have to close all ids manually
    H5Tclose(dtypeId);
    H5Sclose(spaceId);
    H5Dclose(datasetId);
    return result;
}

//TODO: use struct like H5Dataset here, too!
template<typename T>
std::vector<T> readH5Attribute(hid_t fileId, const char *attrName)
{
    std::vector<T> result;
    hsize_t attrSize = 1;
    if (H5Aexists(fileId, attrName)) {
        auto attrId = H5Aopen(fileId, attrName, H5P_DEFAULT);
        if (attrId > -1) {
            auto spaceId = H5Aget_space(attrId);
            if (spaceId > 1) {
                std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(spaceId));
                H5Sget_simple_extent_dims(spaceId, shape.data(), nullptr);

                for (hsize_t dim: shape) {
                    attrSize *= dim;
                }

                result.resize(attrSize);
                if (H5Aread(attrId, H5Aget_type(attrId), result.data()) < 0) {
                    std::cerr << "An error occurred when trying to read the attribute '" << attrName
                              << "' with H5Aread!" << std::endl;
                    H5Sclose(spaceId);
                    H5Aclose(attrId);
                    return std::vector<T>();
                }
            }
            H5Sclose(spaceId);
        } else {
            std::cerr << "An error occurred when trying to open the attribute '" << attrName << "'with H5Aopen."
                      << std::endl;
            return result;
        }
        H5Aclose(attrId);
    } else {
        std::cerr << "The attribute '" << std::to_string(attrName) << "' does not exist in the given file! (H5Aexists)"
                  << std::endl;
    }
    return result;
}

// HDF5 string attributes can be of fixed or variable length, so they have to be read in differently
// than attributes of other types.
template<>
std::vector<std::string> readH5Attribute<std::string>(hid_t fileId, const char *attrName)
{
    std::vector<std::string> result;
    hsize_t attrSize = 1;
    if (H5Aexists(fileId, attrName)) {
        auto attrId = H5Aopen(fileId, attrName, H5P_DEFAULT);
        if (attrId > -1) {
            auto spaceId = H5Aget_space(attrId);
            if (spaceId > 1) {
                std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(spaceId));
                H5Sget_simple_extent_dims(spaceId, shape.data(), nullptr);

                for (hsize_t dim: shape) {
                    attrSize *= dim;
                }

                auto dtypeId = H5Aget_type(attrId);
                assert(H5Tget_class(dtypeId) == H5T_STRING);

                if (H5Tis_variable_str(dtypeId)) {
                    std::vector<char *> tmpResult(attrSize);

                    if (H5Aread(attrId, dtypeId, tmpResult.data()) < 0) {
                        std::cerr << "An error occurred when trying to read the attribute '" << attrName
                                  << "' with H5Aread!" << std::endl;
                        H5Sclose(spaceId);
                        H5Aclose(attrId);
                        return result;
                    }

                    for (hsize_t i = 0; i < attrSize; ++i) {
                        result.push_back(std::string(tmpResult[i]));
                        free(tmpResult[i]);
                    }
                } else {
                    size_t strSize = H5Tget_size(dtypeId);
                    std::vector<char> temp_data(attrSize * strSize);

                    if (H5Aread(attrId, dtypeId, temp_data.data()) < 0) {
                        std::cerr << "An error occurred when trying to read the attribute '" << attrName
                                  << "' with H5Aread!" << std::endl;
                        H5Sclose(spaceId);
                        H5Aclose(attrId);
                        return result;
                    }

                    for (hsize_t i = 0; i < attrSize; ++i) {
                        auto rawStr = std::string(&temp_data[i * strSize], strSize);
                        // remove whitespaces from string
                        rawStr.erase(std::remove_if(rawStr.begin(), rawStr.end(), ::isspace), rawStr.end());
                        result.push_back(rawStr);
                    }
                }
            }
            H5Sclose(spaceId);
        } else {
            std::cerr << "An error occurred when trying to open the attribute '" << attrName << "'with H5Aopen."
                      << std::endl;
            return result;
        }
        H5Aclose(attrId);
    } else {
        std::cerr << "The attribute '" << std::to_string(attrName) << "' does not exist in the given file! (H5Aexists)"
                  << std::endl;
    }
    return result;
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

    auto elemInfo = readH5Dataset<int>(h5Mesh, "ElemInfo").vector;
    if (elemInfo.size() == 0)
        sendError("Could not read in 'ElemInfo' dataset!");

    auto nodeCoords = readH5Dataset<double>(h5Mesh, "NodeCoords").vector;
    if (nodeCoords.size() == 0)
        sendError("Could not read in 'NodeCoords' dataset!");

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

std::map<std::string, Vec<Scalar, 1>::ptr> ReadHopr::getDGSolutionVariables(const char *filename)
{
    std::map<std::string, Vec<Scalar, 1>::ptr> result;
    /*
        DGSolution contains an array of size n * (N + 1)^3 * m, where n is the number of elements
        defined in the mesh, N is the polynomial degree, and m is the number of variables.
    */
    auto h5State = H5Fopen(m_stateFile->getValue().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    if (h5State < 0) {
        sendError("An error occurred while reading in state file. Cannot add data fields to mesh!");
        return result;
    }
    auto DGDataset = readH5Dataset<double>(h5State, "DG_Solution");
    auto DGSolution = DGDataset.vector;

    if (DGSolution.size() == 0)
        sendError("Could not read in 'DG_Solution' dataset!");

    // From general info get:
    // - number of variables in data + their names (e.g., density, momentumX, ...)
    //   Beforehand: create enough output ports (for now let's just read in one)
    auto varNames = readH5Attribute<std::string>(h5State, "VarNames");

    if (varNames.size() == 0)
        sendError("Could not read in 'VarNames' attribute!");

    // - polynomial degree of the solution (N, Ngeo)
    //TODO: think about if we want to keep it like this (works, but hard to read)
    auto N = (readH5Attribute<int>(h5State, "N"))[0];

    // From DG Solutions get:
    // - use algorithm 9 to get the solution at the corner nodes ONLY (no HO nodes for now)
    // TODO: make this work for all element types!
    auto dim = DGDataset.dimension;
    auto nrCorners = 8;
    for (hsize_t varI = 0; varI < dim[4]; varI++) {
        result[varNames[varI]] = Vec<Scalar, 1>::ptr(new Vec<Scalar, 1>(dim[0] * nrCorners));
        auto counter = 0;
        for (hsize_t elemI = 0; elemI < dim[0]; elemI++) {
            for (hsize_t iX = 0; iX < dim[1]; iX++) {
                for (hsize_t iY = 0; iY < dim[2]; iY++) {
                    for (hsize_t iZ = 0; iZ < dim[3]; iZ++) {
                        auto index = varI + dim[4] * (iZ + dim[3] * (iY + dim[2] * (iX + dim[1] * elemI)));
                        auto nodeNr = iZ + (dim[3] * (iY + dim[2] * iX));
                        if ((nodeNr == 0) || (nodeNr == N) || (nodeNr == pow(N + 1, 2) - 1) ||
                            (nodeNr == N * (N + 1)) || (nodeNr == N * pow(N + 1, 2)) ||
                            (nodeNr == N * pow(N + 1, 2) + N) || (nodeNr == pow(N + 1, 3) - 1) ||
                            (nodeNr == N * (N + 1) * (N + 2))) {
                            result[varNames[varI]]->x()[counter] = DGSolution[index];
                            counter++;
                        }
                    }
                }
            }
        }
    }

    // add variable names to field choice parameter
    varNames.insert(varNames.begin(), Invalid);
    for (int i = 0; i < NumPorts; ++i) {
        setParameterChoices(m_fieldChoice[i], varNames);
    }

    H5Fclose(h5State);
    return result;
}


bool ReadHopr::read(Reader::Token &token, int timestep, int block)
{
    UnstructuredGrid::ptr grid;
    std::map<std::string, Vec<Scalar, 1>::ptr> variables;

    auto meshFileName = m_meshFile->getValue();
    if (meshFileName.size()) {
        LOCK_HDF5(comm());
        grid = createMeshFromFile(meshFileName.c_str());
        UNLOCK_HDF5(comm());
    } else {
        sendError("No mesh file was given, so mesh cannot be created.");
        return true;
    }

    auto stateFileName = m_stateFile->getValue();
    if (stateFileName.size()) {
        LOCK_HDF5(comm());
        variables = getDGSolutionVariables(stateFileName.c_str());
        UNLOCK_HDF5(comm());
    } else {
        sendInfo("No state file was given, so no fields will be added to the mesh.");
    }

    updateMeta(grid);
    addObject(m_gridOut, grid);

    for (int i = 0; i < NumPorts; i++) {
        if (m_fieldChoice[i]->getValue() != Invalid) {
            auto varName = m_fieldChoice[i]->getValue();
            auto field = variables[varName];

            if (field) {
                field->addAttribute("_species", varName);
                field->setMapping(vistle::DataBase::Vertex);
                field->setGrid(grid);

                token.applyMeta(field);
                token.addObject(m_fieldsOut[i], field);
            }
        }
    }

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
