#ifndef SPHERES_OVERLAP_H
#define SPHERES_OVERLAP_H

#include <string>
#include <vistle/module/module.h>

/*
    Creates lines between all overlapping spheres.
*/
class SpheresOverlap: public vistle::Module{
public:
    SpheresOverlap(const std::string &name, int moduleID, mpi::communicator comm);
    ~SpheresOverlap();

private:
    vistle::Port *m_spheresIn, *m_linesOut;
    vistle::FloatParameter *m_radiusCoefficient;

    bool compute(const std::shared_ptr<vistle::BlockTask> &task) const override;
};

#endif // SPHERES_OVERLAP_H