#include <string>
#include <sstream>

#include <core/object.h>
#include <core/vec.h>
#include <core/polygons.h>
#include <core/triangles.h>
#include <core/lines.h>
#include <core/points.h>
#include <core/normals.h>
#include <core/unstr.h>

#include <util/coRestraint.h>

// Includes for the CFX application programming interface (API)
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cfxExport.h>  // linked but still qtcreator doesn't find it
#include <getargs.h>    // linked but still qtcreator doesn't find it
#include <iostream>
//#include <unistd.h>

//#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/cstdint.hpp>

#include "ReadCFX.h"

namespace bf = boost::filesystem;

MODULE_MAIN(ReadCFX)

using namespace vistle;


ReadCFX::ReadCFX(const std::string &shmname, const std::string &name, int moduleID)
   : Module("ReadCFX", shmname, name, moduleID) {

    // file browser parameter
    m_resultfiledir = addStringParameter("resultfiledir", "CFX case directory","/mnt/raid/home/hpcjwint/data/cfx/rohr/hlrs_002.res", Parameter::Directory);
    //m_resultfiledir = addStringParameter("resultfiledir", "CFX case directory","/home/jwinterstein/data/cfx/rohr/hlrs_002.res", Parameter::Directory);

    // time parameters
    m_starttime = addFloatParameter("starttime", "start reading at the first step after this time", 0.);
    setParameterMinimum<Float>(m_starttime, 0.);
    m_stoptime = addFloatParameter("stoptime", "stop reading at the last step before this time",
          std::numeric_limits<double>::max());
    setParameterMinimum<Float>(m_stoptime, 0.);
    m_timeskip = addIntParameter("timeskip", "skip this many timesteps after reading one", 0);
    setParameterMinimum<Integer>(m_timeskip, 0);
    m_readGrid = addIntParameter("read_grid", "load the grid?", 1, Parameter::Boolean);

    //zone selection
    m_zoneSelection = addStringParameter("zones","select zone numbers e.g. 1,4,6-10","all");

    // mesh ports
    m_gridOut = createOutputPort("grid_out1");

    // data ports and data choice parameters
    for (int i=0; i<NumPorts; ++i) {
        {// data ports
            std::stringstream s;
            s << "data_out" << i;
            m_volumeDataOut.push_back(createOutputPort(s.str()));
        }
        {// data choice parameters
            std::stringstream s;
            s << "Data" << i;
            auto p =  addStringParameter(s.str(), "name of field", "(NONE)", Parameter::Choice);
            std::vector<std::string> choices;
            choices.push_back("(NONE)");
            setParameterChoices(p, choices);
            m_fieldOut.push_back(p);
        }
    }
    m_readBoundary = addIntParameter("read_boundary", "load the boundary?", 0, Parameter::Boolean);
    //m_boundaryPatchesAsVariants = addIntParameter("patches_as_variants", "create sub-objects with variant attribute for boundary patches", 1, Parameter::Boolean);

    // 2d data ports and 2d data choice parameters
    for (int i=0; i<NumBoundaryPorts; ++i) {
       {// 2d data ports
          std::stringstream s;
          s << "data_2d_out" << i;
          m_boundaryDataOut.push_back(createOutputPort(s.str()));
       }
       {// 2d data choice parameters
          std::stringstream s;
          s << "Data2d" << i;
          auto p =  addStringParameter(s.str(), "name of field", "(NONE)", Parameter::Choice);
          std::vector<std::string> choices;
          choices.push_back("(NONE)");
          setParameterChoices(p, choices);
          m_boundaryOut.push_back(p);
       }
    }
    //m_buildGhostcellsParam = addIntParameter("build_ghostcells", "whether to build ghost cells", 1, Parameter::Boolean);ll


   /*addIntParameter("indexed_geometry", "create indexed geometry?", 0, Parameter::Boolean);
   addIntParameter("triangulate", "only create triangles", 0, Parameter::Boolean);

   addIntParameter("first_block", "number of first block", 0);
   addIntParameter("last_block", "number of last block", 0);

   addIntParameter("first_step", "number of first timestep", 0);
   addIntParameter("last_step", "number of last timestep", 0);
   addIntParameter("step", "increment", 1);

   addIntParameter("ignore_errors", "ignore files that are not found", 0, Parameter::Boolean);
*/
}

ReadCFX::~ReadCFX() {

}

CaseInfo::CaseInfo()
    : m_valid(false) {

}


bool CaseInfo::checkFile(const char *filename) {
    const int MIN_FILE_SIZE = 1024; // minimal size for .res files [in Byte]

    const int MACIC_LEN = 5; // "magic" at the start
    const char *magic = "*INFO";
    char magicBuf[MACIC_LEN];

    boost::uintmax_t fileSize;
    boost::system::error_code ec;

    FILE *fi = fopen(filename, "r");

    if (!fi) {
        std::cout << filename << strerror(errno) << std::endl;
        return false;
    }
    else {
        fileSize = bf::file_size(filename, ec);
        if (ec)
            std::cout << "error code: " << ec << std::endl;
    }

    if (fileSize < MIN_FILE_SIZE) {
        std::cout << filename << "too small to be a real result file" << std::endl;
        std::cout << fileSize << "filesize" << std::endl;
        return false;
    }

    size_t iret = fread(magicBuf, 1, MACIC_LEN, fi);
    if (iret != MACIC_LEN) {
        std::cout << "checkFile :: error reading MACIC_LEN " << std::endl;
        return false;
    }

    if (strncasecmp(magicBuf, magic, MACIC_LEN) != 0) {
        std::cout << filename << "does not start with '*INFO'" << std::endl;
        return false;
    }

    fclose(fi);
    return true;
}

void CaseInfo::getFieldList() {

    int usr_level = 0;
    int dimension, corrected_boundary_node;
    int length;

    int nvars = cfxExportVariableCount(usr_level);
    //std::cerr << "nvars = " << nvars << std::endl;

    m_boundary_param.clear();
    m_boundary_param.push_back("(NONE)");
    m_field_param.clear();
    m_field_param.push_back("(NONE)");

    for(int varnum=1;varnum<=nvars;varnum++) {   //starts from 1 because cfxExportVariableName(varnum,1) only returnes values from 1 and higher
        if(cfxExportVariableSize(varnum,&dimension,&length,&corrected_boundary_node)) { //cfxExportVariableSize returns 1 if successful or 0 if the variable is out of range
            if(length == 1) {
                m_boundary_param.push_back(cfxExportVariableName(varnum,1)); //0 is short form and 1 is long form of the variable name
                //std::cerr << "cfxExportVariableName("<< varnum << ",1) = " << cfxExportVariableName(varnum,1) << std::endl;
            }
            else {
                m_field_param.push_back(cfxExportVariableName(varnum,1)); //0 is short form and 1 is long form of the variable name
            }
        }
    }
}


/*int ReadFOAM::rankForBlock(int processor) const {

   if (m_case.numblocks == 0)
      return 0;

   if (processor == -1)
      return -1;

   return processor % size();
}*/


/*int ReadCFX::rankForBlock(int block) const {

    const int numBlocks = m_lastBlock-m_firstBlock+1;
    if (numBlocks == 0)
        return 0;

    if (block == -1)
        return -1;

    return block % size();
}*/


bool ReadCFX::parameterChanged(const Parameter *p) {
    auto sp = dynamic_cast<const StringParameter *>(p);
    if (sp == m_resultfiledir) {
        std::string c = sp->getValue();
        const char *resultfiledir;
        resultfiledir = c.c_str();

        m_case.m_valid = m_case.checkFile(resultfiledir);
        if (!m_case.m_valid) {
           std::cerr << resultfiledir << " is not a valid CFX .res file" << std::endl;
           return false;
        }
        else {
            sendInfo("Please wait...");

            if (nzones > 0) {
                cfxExportDone();
            }

            char *resultfileName = strdup(resultfiledir);
            nzones = cfxExportInit(resultfileName, NULL);



            if (nzones < 0) {
                cfxExportDone();
                sendError("cfxExportInit could not open %s", resultfileName);
                return false;
            }

            int timeStepNum = cfxExportTimestepNumGet(1);
            if (timeStepNum < 0) {
                sendInfo("no timesteps");
            }

//            std::cerr << "cfxExportTimestepCount()" << cfxExportTimestepCount() << std::endl;
//            std::cerr << "cfxExportTimestepSet(timeStepNum)" << cfxExportTimestepSet(timeStepNum) << std::endl;
//            std::cerr << "cfxExportTimestepNumGet(1)" << cfxExportTimestepNumGet(1) << std::endl;

            cfxExportTimestepSet(timeStepNum);
            if (cfxExportTimestepSet(timeStepNum) < 0) {
                sendInfo("Invalid timestep %d", timeStepNum);
            }

            sendInfo("Found %d zones", nzones);
            sendInfo("The initialisation was successfully done");

            //fill choice parameter
            m_case.getFieldList();
            for (auto out: m_fieldOut) {
               setParameterChoices(out, m_case.m_field_param);
            }
            for (auto out: m_boundaryOut) {
                setParameterChoices(out, m_case.m_boundary_param);
            }
//            std::vector<std::string>::const_iterator it1;
//            for(it1=m_case.m_field_param.begin();it1!=m_case.m_field_param.end();++it1) {
//                std::cerr << "field_choices = " << *it1 << std::endl;
//            }
//            for(it1=m_case.m_boundary_param.begin();it1!=m_case.m_boundary_param.end();++it1) {
//                std::cerr << "field_choices = " << *it1 << std::endl;
//            }

            //print out zone names
            for(int i=1;i<=nzones;i++) {
                cfxExportZoneSet(i,NULL);
                sendInfo("name of zone no. %d: %s \n",i,cfxExportZoneName(i));
                cfxExportZoneFree();
            }

            free(resultfileName);
        }
    }

    return Module::parameterChanged(p);
}

void ReadCFX::loadGrid() {
    bool readGrid = m_readGrid->getValue();
    if(readGrid) {
        UnstructuredGrid::ptr grid(new UnstructuredGrid(0, 0, 0)); //initialized with number of elements, number of connectivities, number of coordinates

        cfxExportZoneSet(0,NULL); //0 means global zone (all zones)

        //load nodes into unstructured grid
        double *x_coord = new double, *y_coord = new double, *z_coord = new double;

        nnodes = cfxExportNodeCount();
        nodeList = cfxExportNodeList(); //load coordinates into array with structs, each struct containing x,y,z of one node
        for(int nodeid=1;nodeid<=nnodes;++nodeid) {
            if(!cfxExportNodeGet(nodeid,x_coord,y_coord,z_coord)) {  //get access to coordinates
                std::cerr << "error while reading nodes out of .res file: nodeid is out of range" << std::endl;
            }
            grid->x().push_back(*(x_coord));
            grid->y().push_back(*(y_coord));
            grid->z().push_back(*(z_coord));
        }


        //Test, ob Einlesen funktioniert hat
        std::cerr << "nnodes = " << nnodes << std::endl;
        std::cerr << "grid->getNumCoords()" << grid->getNumCoords() << std::endl;
        std::cerr << "x,y,z (1)" << grid->x().at(1) << ", " << grid->y().at(1) << ", " << grid->z().at(1) << std::endl;
        std::cerr << "x,y,z (10)" << grid->x().at(10) << ", " << grid->y().at(10) << ", " << grid->z().at(10) << std::endl;
        std::cerr << "x,y,z (nnodes-1)" << grid->x().at(nnodes-1) << ", " << grid->y().at(nnodes-1) << ", " << grid->z().at(nnodes-1) << std::endl;

        cfxExportNodeFree();
        delete[] x_coord;
        delete[] y_coord;
        delete[] z_coord;
    }




    //grid->tl.

//            if (cfxExportZoneSet (zone, counts) < 0)
//                     cfxExportFatal ("invalid zone number");
//            nnodes = cfxExportNodeCount();
//            nelems = cfxExportElementCount();
//            if (counts[cfxCNT_PYR]) {
//               printf ("%d pyramid elements found - they are being ignored\n",
//                  counts[cfxCNT_PYR]);
//               nelems -= counts[cfxCNT_PYR];
//            }
//            if (!nnodes || !nelems)
//               cfxExportFatal ("no nodes and/or elements");
    return;
}

bool ReadCFX::compute() {

    std::cerr << "Compute Start. \n";


    // read zone selection; m_zonesSelected contains a bool array of which zones are selected

    //m_zone(zone) = 1 Zone ist selektiert
    //m_zone(zone) = 0 Zone ist nicht selektiert
    //group = -1 alle Zonen sind selektiert
    m_zonesSelected.clear();
    m_zonesSelected.add(m_zoneSelection->getValue());
    ssize_t val = m_zonesSelected.getNumGroups();
    ssize_t group;
    m_zonesSelected.get(val,group);

    loadGrid();



    //    //write nodes

//    nnodes = cfxExportNodeCount();
//    std::cerr << "nnodes = " << nnodes << std::endl;

//    nodeList = cfxExportNodeList();

//    for(int n=0;n<10;n++,nodeList++) {

//        std::cerr << "x = " << nodeList->x << " y = " << nodeList->y << " z = " << nodeList->z << std::endl;
//    }

//    cfxExportNodeFree();

//    //write elements
//    nelems = cfxExportElementCount();
//    std::cerr << "nelems = " << nelems << std::endl;

//    if (counts[cfxCNT_TET]) {
//           elems = cfxExportElementList();
//           for (int n = 0; n < 1; n++, elems++) {
//               if (cfxELEM_TET == elems->type) {
//                   for (int i = 0; i < elems->type; i++)
//                       std::cerr << "elems = " << elems->nodeid[i] << std::endl;
//               }
//           }
//       }



//    cfxExportElementFree();



 /*  m_firstBlock = getIntParameter("first_block");
   m_lastBlock = getIntParameter("last_block");
   m_firstStep = getIntParameter("first_step");
   m_lastStep = getIntParameter("last_step");
   m_step = getIntParameter("step");
   if ((m_lastStep-m_firstStep)*m_step < 0)
       m_step *= -1;
   bool reverse = false;
   int step = m_step, first = m_firstStep, last = m_lastStep;
   if (step < 0) {
       reverse = true;
       step *= -1;
       first *= -1;
       last *= -1;
   }

   std::string filename = getStringParameter("filename");

   int timeCounter = 0;
   for (int t=first; t<=last; t += step) {
       int timestep = reverse ? -t : t;
       bool loaded = false;
       for (int b=m_firstBlock; b<=m_lastBlock; ++b) {
           if (rankForBlock(b) == rank()) {
               std::string f;
               try {
                   using namespace boost::io;
                   boost::format fmter(filename);
                   fmter.exceptions(all_error_bits ^ (too_many_args_bit | too_few_args_bit));
                   fmter % b;
                   fmter % timestep;
                   f = fmter.str();
               } catch (boost::io::bad_format_string except) {
                   sendError("bad format string in filename");
                   return true;
               }
               auto obj = load(f);
               if (!obj) {
                   if (!getIntParameter("ignore_errors")) {
                       sendError("failed to load %s", f.c_str());
                       return true;
                   }
               } else {
                   obj->setBlock(b);
                   obj->setTimestep(timeCounter);
                   loaded = true;
                   addObject("grid_out", obj);
               }
           }
       }
       if (loaded)
           ++timeCounter;
   }*/

    cfxExportDone();

   return true;
}
