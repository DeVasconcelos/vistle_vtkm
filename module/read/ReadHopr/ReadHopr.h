#ifndef READHOPR_H
#define READHOPR_H

#include <vistle/module/reader.h>

namespace vistle {

class ReadHopr: public vistle::Reader {
public:
    ReadHopr(const std::string &name, int moduleID, mpi::communicator comm);
    ~ReadHopr() override;

    bool examine(const vistle::Parameter *param) override;
    bool read(vistle::Reader::Token &token, int timestep = -1, int block = -1) override;
    bool prepareRead() override;
    bool finishRead() override;

private:
    vistle::Port *m_gridOut;

    vistle::StringParameter *m_meshFile;
    vistle::StringParameter *m_stateFile;

    // Reads in the 'NodeCoords' and 'ElemInfo' datasets stored in the HOPR file
    // and uses the information inside to create an unstructured vistle grid.
    vistle::UnstructuredGrid::ptr createMeshFromFile(const char *filename);
};
} // namespace vistle

#endif //READHOPR_H
