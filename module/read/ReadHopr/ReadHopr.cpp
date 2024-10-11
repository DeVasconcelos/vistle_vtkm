#include <algorithm>
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
#include "utils/ReadHDF5.h"

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

void ReadHopr::setFieldChoices(const std::vector<std::string> &choices)
{
    auto choicesPlusInvalid = choices;
    choicesPlusInvalid.insert(choicesPlusInvalid.begin(), Invalid);
    for (int i = 0; i < NumPorts; i++) {
        setParameterChoices(m_fieldChoice[i], choicesPlusInvalid);
    }
}

class CGNSCorners {
private:
    std::map<vistle::Byte, std::vector<std::vector<hsize_t>>> cornerIndicesMap;

public:
    CGNSCorners(hsize_t N)
    {
        cornerIndicesMap = {
            {UnstructuredGrid::TETRAHEDRON, {{0, 0, 0}, {N, 0, 0}, {0, N, 0}, {0, 0, N}}},
            {UnstructuredGrid::PYRAMID, {{0, 0, 0}, {N, 0, 0}, {N, N, 0}, {0, N, 0}, {0, 0, N}}},
            {UnstructuredGrid::PRISM, {{0, 0, 0}, {N, 0, 0}, {0, N, 0}, {0, 0, N}, {N, 0, N}, {0, N, N}}},
            {UnstructuredGrid::HEXAHEDRON,
             {{0, 0, 0}, {N, 0, 0}, {N, N, 0}, {0, N, 0}, {0, 0, N}, {N, 0, N}, {N, N, N}, {0, N, N}}}};
    }

    bool isCorner(std::vector<hsize_t> index, vistle::Byte type) const
    {
        if (cornerIndicesMap.find(type) != cornerIndicesMap.end()) {
            auto cornerIndices = cornerIndicesMap.at(type);
            return std::find(cornerIndices.begin(), cornerIndices.end(), index) != cornerIndices.end();
        } else {
            throw exception(
                "Encountered unsupported element type when trying to read out corner nodes in the given HOPR "
                "dataset.");
        }
    }

    std::vector<std::vector<hsize_t>> getIndices(vistle::Byte type) const
    {
        if (cornerIndicesMap.find(type) != cornerIndicesMap.end()) {
            return cornerIndicesMap.at(type);
        } else {
            throw exception(
                "Encountered unsupported element type when trying to read out corner nodes in the given HOPR "
                "dataset.");
        }
    }
};


// Using Algorithm 9 to get the solution at the CGNS corner nodes only
// (see Hopr documentation, https://hopr.readthedocs.io/en/latest/userguide/meshformat.html)
std::map<std::string, Vec<Scalar, 1>::ptr>
getSolutionDataAtCornerNodes(std::vector<double> DGSolution, std::vector<hsize_t> DGDim, hsize_t polynomialDegree,
                             const std::vector<std::string> &varNames, const Byte *typeList, Index numCorners)
{
    auto cornerNodes = CGNSCorners(polynomialDegree);

    std::map<std::string, Vec<Scalar, 1>::ptr> result;
    for (hsize_t varI = 0; varI < DGDim[4]; varI++) {
        result[varNames[varI]] = Vec<Scalar, 1>::ptr(new Vec<Scalar, 1>(numCorners));
        auto counter = 0;
        for (hsize_t elemI = 0; elemI < DGDim[0]; elemI++) {
            for (auto &cornerIndex: cornerNodes.getIndices(typeList[elemI])) {
                auto [iX, iY, iZ] = std::tie(cornerIndex[0], cornerIndex[1], cornerIndex[2]);
                result[varNames[varI]]->x()[counter] =
                    DGSolution[varI + DGDim[4] * (iZ + DGDim[3] * (iY + DGDim[2] * (iX + DGDim[1] * elemI)))];
                counter++;
            }
        }
    }

    return result;
}

std::map<std::string, Vec<Scalar, 1>::ptr> ReadHopr::extractFieldsFromStateFile(const char *filename,
                                                                                const Byte *typeList, Index numCorners)
{
    auto h5State = H5Fopen(m_stateFile->getValue().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (h5State < 0) {
        sendError("An error occurred while reading in state file. Cannot add data fields to mesh!");
        return std::map<std::string, Vec<Scalar, 1>::ptr>();
    }

    auto varNames = readH5Attribute<std::string>(h5State, "VarNames");

    if (varNames.size() == 0)
        sendError("Could not read in 'VarNames' attribute!");

    setFieldChoices(varNames);

    //TODO: what's the difference between N and Ngeo?
    auto N = (readH5Attribute<hsize_t>(h5State, "N"))[0];

    auto DGDataset = readH5Dataset<double>(h5State, "DG_Solution");
    auto DGSolution = DGDataset.vector;
    auto DGDim = DGDataset.dimension;

    if (DGSolution.size() == 0)
        sendError("Could not read in 'DG_Solution' dataset!");

    auto result = getSolutionDataAtCornerNodes(DGSolution, DGDim, N, varNames, typeList, numCorners);

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
        variables = extractFieldsFromStateFile(stateFileName.c_str(), grid->tl().data(), grid->getNumCorners());
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
