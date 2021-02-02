/**************************************************************************\
 **                                                                        **
 **                                                                        **
 ** Description: Read module for ChEESE tsunami nc-files                   **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 ** Author:    Marko Djuric                                                **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 ** Date:  25.01.2021                                                      **
\**************************************************************************/

//vistle
#include "vistle/core/database.h"
#include "vistle/core/index.h"
#include "vistle/core/object.h"
#include "vistle/core/parameter.h"
#include "vistle/core/polygons.h"
#include "vistle/core/scalar.h"
#include "vistle/core/vector.h"
#include "vistle/module/module.h"
#include "vistle/module/reader.h"

//header
#include "ReadTsunami.h"

//std
#include <cstddef>
#include <string>
#include <vector>

MODULE_MAIN(ReadTsunami)

using namespace vistle;

ReadTsunami::ReadTsunami(const std::string &name, int moduleID, mpi::communicator comm)
: vistle::Reader("Read ChEESE Tsunami files", name, moduleID, comm)
{
    ncDataFile = nullptr;

    // define parameters

    // file-browser
    m_filedir = addStringParameter("file_dir", "NC File directory", "/data/ChEESE/Tsunami/pelicula_eta.nc",
                                   Parameter::Filename);

    // visualise variables
    m_verticalScale = addFloatParameter("VerticalScale", "Vertical Scale parameter", 1.0);
    /* m_step = addIntParameter("StepWidth", "Timestep step width", 1); */
    /* setParameterRange(m_step, Integer(0), Integer(999999)); */

    // define ports

    // 2D Surface
    m_surface_out = createOutputPort("surfaceOut", "2D Grid output (Polygons)");
    m_seaSurface_out = createOutputPort("seaSurfaceOut", "2D See floor (Polygons)");
    m_waterSurface_out = createOutputPort("waterSurfaceOut", "2D water surface (Polygons)");
    m_maxHeight = createOutputPort("maxHeight", "Max water height (Float)");

    /* m_ghostLayerWidth = addIntParameter("ghost_layers", "number of ghost layers on all sides of a grid", 0); */

    //observer parameters
    observeParameter(m_filedir);
    /* observeParameter(m_blocks[0]); */
    /* observeParameter(m_blocks[1]); */
    /* observeParameter(m_blocks[2]); */
    observeParameter(m_step);

    /* setParallelizationMode(ParallelizeBlocks); */
    setParallelizationMode(Serial);
}

ReadTsunami::~ReadTsunami()
{}

bool ReadTsunami::prepareRead()
{
    if (!openNcFile())
        return false;
    return true;
}

/**
 * Open Nc File and set pointer ncDataFile.
 *
 * @return true if its not empty or cannot be opened.
 */
bool ReadTsunami::openNcFile()
{
    std::string sFileName = m_filedir->getValue();

    if (sFileName.empty()) {
        sendInfo("NetCDF filename is empty!");
        return false;
    } else {
        try {
            ncDataFile = new NcFile(sFileName.c_str(), NcFile::read);
            sendInfo("Reading File: " + sFileName);
        } catch (...) {
            sendInfo("Couldn't open NetCDF file!");
            return false;
        }

        if (ncDataFile->getVarCount() == 0) {
            sendInfo("empty NetCDF file!");
            return false;
        } else {
            return true;
        }
    }
}

bool ReadTsunami::examine(const vistle::Parameter *param)
{
    /* size_t nblocks = m_blocks[0]->getValue() * m_blocks[1]->getValue() * m_blocks[2]->getValue(); */
    /* setPartitions(nblocks); */
    /* setTimesteps(m_step->getValue()); */
    return true;
    /* return false; */
    /* return nblocks > 0; */
}

void ReadTsunami::fillCoords2DimPoly(Polygons::ptr poly, const size_t &dimX, const size_t &dimY,
                                     const std::vector<float *> &coords)
{
    int n = 0;
    auto sx_coord = poly->x().data(), sy_coord = poly->y().data(), sz_coord = poly->z().data();
    for (size_t i = 0; i < dimX; i++)
        for (size_t j = 0; j < dimY; j++, n++) {
            sx_coord[n] = coords.at(0)[i];
            sy_coord[n] = coords.at(1)[j];
            sz_coord[n] = 0;
        }
}

void ReadTsunami::fillCoords3DimPoly(Polygons::ptr poly, const size_t &dimX, const size_t &dimY, const size_t &dimZ,

                                     const std::vector<float *> &coords)
{
    int n = 0;
    auto sx_coord = poly->x().data(), sy_coord = poly->y().data(), sz_coord = poly->z().data();
    for (size_t i = 0; i < dimX; i++)
        for (size_t j = 0; j < dimY; j++)
            for (size_t k = 0; k < dimZ; k++, n++) {
                sx_coord[n] = coords.at(0)[i];
                sy_coord[n] = coords.at(1)[j];
                sz_coord[n] = coords.at(2)[k];
            }
}

void ReadTsunami::fillConnectList2DimPoly(Polygons::ptr poly, const size_t &dimX, const size_t &dimY)
{
    int n = 0;
    auto surfaceConnectivityList = poly->cl().data();
    for (size_t j = 1; j < dimX; j++)
        for (size_t k = 1; k < dimY; k++) {
            surfaceConnectivityList[n++] = (j - 1) * dimY + (k - 1);
            surfaceConnectivityList[n++] = j * dimY + (k - 1);
            surfaceConnectivityList[n++] = j * dimY + k;
            surfaceConnectivityList[n++] = (j - 1) * dimY + k;
        }
}

void ReadTsunami::fillPolyList4Corner(Polygons::ptr poly, const size_t &numPoly)
{
    auto polyList = poly->el().data();
    for (size_t p = 0; p < numPoly; p++)
        polyList[p] = p * 4;
}

/**
 * Generate surface from polygons.
 *
 * @numElem Number of polygons that will be used to great the surface.
 * @numCorner Number of all corners.
 * @numVertices Number of all different corners.
 * @dimension dimension in x,y,z.
 * @coords coordinates for polygons.
 * @return vistle::Polygons::ptr
 */
Polygons::ptr ReadTsunami::generateSurface(const size_t &numElem, const size_t &numCorner, const size_t &numVertices,
                                           const std::vector<size_t> &dimension, const std::vector<float *> &coords)
{
    Polygons::ptr surface(new Polygons(numElem, numCorner, numVertices));

    size_t surfaceDimX = dimension.at(0);
    size_t surfaceDimY = dimension.at(1);
    size_t surfaceDimZ = dimension.at(2);

    if (surfaceDimZ == 0) {
        // fill coords 2D
        fillCoords2DimPoly(surface, surfaceDimX, surfaceDimY, coords);
    } else {
        // fill coords 3D
        fillCoords3DimPoly(surface, surfaceDimX, surfaceDimY, surfaceDimZ, coords);
    }

    // fill vertices
    fillConnectList2DimPoly(surface, surfaceDimX, surfaceDimY);

    // fill the polygon list
    fillPolyList4Corner(surface, numElem);

    return surface;
}

void ReadTsunami::block(Reader::Token &token, Index bx, Index by, Index bz, vistle::Index block,
                        vistle::Index time) const
{
    /* // read variables from NetCDF-File */
    /* NcVar latvar = ncDataFile->getVar("lat"); */
    /* NcVar lonvar = ncDataFile->getVar("lon"); */
    /* NcVar grid_latvar = ncDataFile->getVar("grid_lat"); */
    /* NcVar grid_lonvar = ncDataFile->getVar("grid_lon"); */
    /* NcVar bathymetryvar = ncDataFile->getVar("bathymetry"); */
    /* NcVar max_height = ncDataFile->getVar("max_height"); */
    /* NcVar eta = ncDataFile->getVar("eta"); */

    /* // dimension from lat and lon variables */
    /* size_t surfaceDimX{latvar.getDim(0).getSize()}; */
    /* size_t surfaceDimY{lonvar.getDim(0).getSize()}; */
    /* size_t surfaceDimZ{0}; */
    /* std::vector dimSurface{surfaceDimX, surfaceDimY, surfaceDimZ}; */

    /* // num of polygons */
    /* size_t surfaceNumPoly = (surfaceDimX - 1) * (surfaceDimY - 1); */

    /* // pointer for lat values and coords */
    /* float *latVals = new float[surfaceDimX]; */
    /* float *lonVals = new float[surfaceDimY]; */
    /* std::vector coords{latVals, lonVals}; */

    /* // read in lat var ncdata into float-pointer */
    /* latvar.getVar(latVals); */
    /* lonvar.getVar(lonVals); */

    /* // Now for the 2D variables, we create a surface */
    /* Polygons::ptr surfacePolygons = */
    /*     generateSurface(surfaceNumPoly, surfaceNumPoly * 4, surfaceDimX * surfaceDimY, dimSurface, coords); */

    /* // Delete buffers from grid replication */
    /* delete[] latVals; */
    /* delete[] lonVals; */

    /* // get dim from grid_lon & grid_lat */
    /* size_t gridLatDimX = grid_latvar.getDim(0).getSize(); */
    /* size_t gridLonDimY = grid_lonvar.getDim(0).getSize(); */
    /* latVals = new float[gridLatDimX]; */
    /* lonVals = new float[gridLonDimY]; */

    /* // depth */
    /* float *depthVals = new float[gridLatDimX * gridLonDimY]; */

    /* // set where to stream to data (float pointer) */
    /* grid_latvar.getVar(latVals); */
    /* grid_lonvar.getVar(lonVals); */
    /* bathymetryvar.getVar(depthVals); */

    /* // Now for the 2D variables, we create a surface */
    /* size_t numPolygons = (gridLatDimX - 1) * (gridLonDimY - 1); */
    /* Polygons::ptr polygonSeaSurface(new Polygons(numPolygons, numPolygons * 4, gridLatDimX * gridLonDimY)); */

    /* // Fill the _coord arrays (memcpy faster?) */
    /* int n{0}; */
    /* auto x_coord = polygonSeaSurface->x().data(), y_coord = polygonSeaSurface->y().data(), */
    /*      z_coord = polygonSeaSurface->z().data(); */
    /* for (size_t j = 0; j < gridLatDimX; j++) */
    /*     for (size_t k = 0; k < gridLonDimY; k++, n++) { */
    /*         x_coord[n] = latVals[j]; */
    /*         y_coord[n] = lonVals[k]; */
    /*         z_coord[n] = zCalcSeaSurface(depthVals, j, k, gridLonDimY); */
    /*     } */

    /* // Fill the connectivitylist list = numPolygons * 4 */
    /* fillConnectList2DimPoly(polygonSeaSurface, gridLatDimX, gridLonDimY); */

    /* // Fill the polygon list */
    /* fillPolyList4Corner(polygonSeaSurface, numPolygons); */

    /* // Delete buffers from grid replication */
    /* delete[] latVals; */
    /* delete[] lonVals; */
    /* delete[] depthVals; */

    /* /1* // float for read in max_height *1/ */
    /* /1* Vec<Scalar>::ptr scalarMaxHeight(new Vec<Scalar>(max_height.getDimCount())); *1/ */
    /* /1* vistle::Scalar *maxH{scalarMaxHeight->x().data()}; *1/ */
    /* /1* max_height.getVar(maxH); *1/ */

    /* // read in time */
    /* /1* float *floatData = new float[eta.getDim(0).getSize() * eta.getDim(1).getSize() * eta.getDim(2).getSize()]; *1/ */

    /* /1* int numTimesteps = eta.getDim(0).getSize(); *1/ */
    /* /1* eta.getVar(floatData); *1/ */

    /* // create Object::ptr pointer for each timestep */
    /* /1* std::string baseName = m_waterSurface_out->getName(); *1/ */
    /* /1* sendInfo("numTimesteps %d!", numTimesteps); *1/ */
    /* /1* int i = 0; *1/ */

    /* // get vertical Scale */
    /* /1* float zScale = m_verticalScale->getValue(); *1/ */

    /* /1* // create watersurface with polygons for each timestep *1/ */
    /* /1* for (int t = 0; t < numTimesteps; t += m_step->getValue()) { *1/ */

    /* /1*     // Now for the 2D variables, we create a surface *1/ */
    /* /1*     int snumPolygons = (surfaceDimX - 1) * (surfaceDimY - 1); *1/ */

    /* /1*     Polygons::ptr outSurface(new Polygons(numPolygons, numPolygons * 4, gridLatDimX * gridLonDimY)); *1/ */

    /* /1*     auto x_coord = outSurface->x().data(), y_coord = outSurface->y().data(), z_coord = outSurface->z().data(); *1/ */
    /* /1*     for (auto n = 0; n < surfaceDimX * surfaceDimY; n++) { *1/ */
    /* /1*         x_coord[n] = sx_coord[n]; *1/ */
    /* /1*         y_coord[n] = sy_coord[n]; *1/ */
    /* /1*         z_coord[n] = floatData[t * surfaceDimX * surfaceDimY + n] * zScale; *1/ */
    /* /1*     } *1/ */

    /* /1*     auto vl = outSurface->cl().data(); *1/ */
    /* /1*     for (int j = 0; j < snumPolygons * 4; j++) *1/ */
    /* /1*         vl[j] = surfaceConnectivityList[j]; *1/ */

    /* /1*     auto pl = outSurface->el().data(); *1/ */
    /* /1*     for (int p = 0; p < snumPolygons; p++) *1/ */
    /* /1*         pl[p] = surfacePolygonList[p]; *1/ */

    /* /1*     token.addObject(m_surface_out, outSurface); *1/ */
    /* /1*     i++; *1/ */
    /* /1* } *1/ */

    /* /1* surfacePolygons->updateInternals(); *1/ */
    /* /1* surfacePolygons->setBlock(block); *1/ */
    /* /1* if (time >= 0) *1/ */
    /*     /1* surfacePolygons->setTimestep(time); *1/ */

    /* /1* Vec<Scalar>::ptr dataOutput; *1/ */
    /* /1* dataOutput->setGrid(surfacePolygons); *1/ */
    /* /1* dataOutput->setMapping(DataBase::Vertex); *1/ */
    /* /1* dataOutput->addAttribute("_species", ) *1/ */
    /* /1* token.addObject(m_surface_out, dataOutput); *1/ */
    /* token.addObject(m_surface_out, surfacePolygons); */

    /* token.addObject(m_seaSurface_out, polygonSeaSurface); */
    /* delete floatData; */
    /* delete maxH; */
}

bool ReadTsunami::read(Token &token, int timestep, int block)
{
    /* Index steps = m_step->getValue(); */
    /* if (steps > 0 && timestep < 0) { */
    /*     // don't generate constant data when animation has been requested */
    /*     return true; */
    /* } */

    /* Index blocks[3]; */
    /* for (int i = 0; i < 3; ++i) { */
    /*     blocks[i] = m_blocks[i]->getValue(); */
    /* } */

    /* Index b = blockNum; */
    /* Index bx = b % blocks[0]; */
    /* b /= blocks[0]; */
    /* Index by = b % blocks[1]; */
    /* b /= blocks[1]; */
    /* Index bz = b; */

    /* block(token, bx, by, bz, blockNum, timestep); */

    // read variables from NetCDF-File
    NcVar latvar = ncDataFile->getVar("lat");
    NcVar lonvar = ncDataFile->getVar("lon");
    NcVar grid_latvar = ncDataFile->getVar("grid_lat");
    NcVar grid_lonvar = ncDataFile->getVar("grid_lon");
    NcVar bathymetryvar = ncDataFile->getVar("bathymetry");
    NcVar max_height = ncDataFile->getVar("max_height");
    NcVar eta = ncDataFile->getVar("eta");

    // dimension from lat and lon variables
    size_t surfaceDimX{latvar.getDim(0).getSize()};
    size_t surfaceDimY{lonvar.getDim(0).getSize()};
    size_t surfaceDimZ{0};
    std::vector dimSurface{surfaceDimX, surfaceDimY, surfaceDimZ};

    // num of polygons
    size_t surfaceNumPoly = (surfaceDimX - 1) * (surfaceDimY - 1);

    // pointer for lat values and coords
    float *latVals = new float[surfaceDimX];
    float *lonVals = new float[surfaceDimY];
    std::vector coords{latVals, lonVals};

    // read in lat var ncdata into float-pointer
    latvar.getVar(latVals);
    lonvar.getVar(lonVals);

    // Now for the 2D variables, we create a surface
    Polygons::ptr surfacePolygons =
        generateSurface(surfaceNumPoly, surfaceNumPoly * 4, surfaceDimX * surfaceDimY, dimSurface, coords);

    // Delete buffers from grid replication
    delete[] latVals;
    delete[] lonVals;

    // get dim from grid_lon & grid_lat
    size_t gridLatDimX = grid_latvar.getDim(0).getSize();
    size_t gridLonDimY = grid_lonvar.getDim(0).getSize();
    latVals = new float[gridLatDimX];
    lonVals = new float[gridLonDimY];

    // depth
    float *depthVals = new float[gridLatDimX * gridLonDimY];

    // set where to stream to data (float pointer)
    grid_latvar.getVar(latVals);
    grid_lonvar.getVar(lonVals);
    bathymetryvar.getVar(depthVals);

    // Now for the 2D variables, we create a surface
    int numPolygons = (gridLatDimX - 1) * (gridLonDimY - 1);
    Polygons::ptr polygonSeaSurface(new Polygons(numPolygons, numPolygons * 4, gridLatDimX * gridLonDimY));

    // Fill the _coord arrays (memcpy faster?)
    int n{0};
    auto x_coord = polygonSeaSurface->x().data(), y_coord = polygonSeaSurface->y().data(),
         z_coord = polygonSeaSurface->z().data();
    for (size_t j = 0; j < gridLatDimX; j++)
        for (size_t k = 0; k < gridLonDimY; k++, n++) {
            x_coord[n] = latVals[j];
            y_coord[n] = lonVals[k];
            z_coord[n] = zCalcSeaSurface(depthVals, j, k, gridLonDimY);
        }

    // Fill the connectivitylist list = numPolygons * 4
    fillConnectList2DimPoly(polygonSeaSurface, gridLatDimX, gridLonDimY);

    // Fill the polygon list
    fillPolyList4Corner(polygonSeaSurface, numPolygons);

    // Delete buffers from grid replication
    delete[] latVals;
    delete[] lonVals;
    delete[] depthVals;

    /* // float for read in max_height */
    /* Vec<Scalar>::ptr scalarMaxHeight(new Vec<Scalar>(max_height.getDimCount())); */
    /* vistle::Scalar *maxH{scalarMaxHeight->x().data()}; */
    /* max_height.getVar(maxH); */

    // read in time
    /* float *floatData = new float[eta.getDim(0).getSize() * eta.getDim(1).getSize() * eta.getDim(2).getSize()]; */

    /* int numTimesteps = eta.getDim(0).getSize(); */
    /* eta.getVar(floatData); */

    // create Object::ptr pointer for each timestep
    /* std::string baseName = m_waterSurface_out->getName(); */
    /* sendInfo("numTimesteps %d!", numTimesteps); */
    /* int i = 0; */

    // get vertical Scale
    /* float zScale = m_verticalScale->getValue(); */

    /* // create watersurface with polygons for each timestep */
    /* for (int t = 0; t < numTimesteps; t += m_step->getValue()) { */

    /*     // Now for the 2D variables, we create a surface */
    /*     int snumPolygons = (surfaceDimX - 1) * (surfaceDimY - 1); */

    /*     Polygons::ptr outSurface(new Polygons(numPolygons, numPolygons * 4, gridLatDimX * gridLonDimY)); */

    /*     auto x_coord = outSurface->x().data(), y_coord = outSurface->y().data(), z_coord = outSurface->z().data(); */
    /*     for (auto n = 0; n < surfaceDimX * surfaceDimY; n++) { */
    /*         x_coord[n] = sx_coord[n]; */
    /*         y_coord[n] = sy_coord[n]; */
    /*         z_coord[n] = floatData[t * surfaceDimX * surfaceDimY + n] * zScale; */
    /*     } */

    /*     auto vl = outSurface->cl().data(); */
    /*     for (int j = 0; j < snumPolygons * 4; j++) */
    /*         vl[j] = surfaceConnectivityList[j]; */

    /*     auto pl = outSurface->el().data(); */
    /*     for (int p = 0; p < snumPolygons; p++) */
    /*         pl[p] = surfacePolygonList[p]; */

    /*     token.addObject(m_surface_out, outSurface); */
    /*     i++; */
    /* } */

    /* surfacePolygons->updateInternals(); */
    /* surfacePolygons->setBlock(block); */
    /* if (time >= 0) */
    /* surfacePolygons->setTimestep(time); */

    Vec<Scalar,3>::ptr dataOutput(new Vec<Scalar,3>(surfaceDimY*surfaceDimX));
    dataOutput->setGrid(surfacePolygons);
    dataOutput->setTimestep(timestep);
    dataOutput->setMapping(DataBase::Vertex);
    /* dataOutput->addAttribute("_species", "test"); */
    token.addObject(m_surface_out, dataOutput);
    /* token.addObject(m_surface_out, surfacePolygons); */

    /* token.addObject(m_surface_out, surfacePolygons); */
    /* token.addObject(m_seaSurface_out, polygonSeaSurface); */
    /* delete floatData; */
    /* delete maxH; */

    return true;
}

bool ReadTsunami::finishRead()
{
    if (ncDataFile) {
        delete ncDataFile;
        ncDataFile = nullptr;
    }
    return true;
}