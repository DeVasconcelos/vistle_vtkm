#include "ReadHopr.h"

MODULE_MAIN(ReadHopr)
using vistle::Parameter;
using vistle::Reader;

ReadHopr::ReadHopr(const std::string &name, int moduleID, mpi::communicator comm): Reader(name, moduleID, comm)
{}

ReadHopr::~ReadHopr()
{}

bool ReadHopr::examine(const vistle::Parameter *param)
{
    return true;
}
bool ReadHopr::read(vistle::Reader::Token &token, int timestep, int block)
{
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
