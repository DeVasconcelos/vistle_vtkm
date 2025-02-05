#ifndef VISTLE_VTKM_CLIP_H
#define VISTLE_VTKM_CLIP_H

#include <vistle/module/module.h>
#include <vistle/vtkm/ImplFuncController.h>

class ClipVtkm: public vistle::Module {
public:
    ClipVtkm(const std::string &name, int moduleID, mpi::communicator comm);
    ~ClipVtkm();

private:
    vistle::Port *m_dataOut = nullptr;

    vistle::IntParameter *m_flip = nullptr;

    bool compute(const std::shared_ptr<vistle::BlockTask> &task) const override;
    bool changeParameter(const vistle::Parameter *param) override;

    ImplFuncController m_implFuncControl;
};

#endif // VISTLE_VTKM_ISOSURFACE_H
