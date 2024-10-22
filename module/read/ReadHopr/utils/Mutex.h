#ifndef READHOPR_MUTEX_H
#define READHOPR_MUTEX_H

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

#endif //READHOPR_MUTEX_H
