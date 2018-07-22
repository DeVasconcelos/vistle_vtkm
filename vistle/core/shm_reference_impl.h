#ifndef SHM_REFERENCE_IMPL_H
#define SHM_REFERENCE_IMPL_H

#include "archives_config.h"

namespace vistle {

template<class T>
template<class Archive>
void shm_ref<T>::save(Archive &ar) const {

    ar & V_NAME(ar, "shm_name", m_name);
    ar.template saveArray<typename T::value_type>(*this);
}

template<class T>
template<class Archive>
void shm_ref<T>::load(Archive &ar) {
   ar & V_NAME(ar, "shm_name", m_name);
   std::string name = m_name;

   unref();
   m_p = nullptr;
   auto obj = ar.currentObject();
   if (obj && !valid()) {
       //std::cerr << "obj " << obj->name << ": unresolved: " << name << std::endl;
       obj->unresolvedReference();
   }
   auto handler = ar.objectCompletionHandler();
   auto ref =  ar.template getArray<typename T::value_type>(name, [this, name, obj, handler]() -> void {
      //std::cerr << "array completion handler: " << name << std::endl;
      auto ref = Shm::the().getArrayFromName<typename T::value_type>(name);
      assert(ref);
      *this = ref;
      if (obj) {
         //std::cerr << "obj " << obj->name << ": RESOLVED: " << name << std::endl;
         obj->referenceResolved(handler);
      } else {
         //std::cerr << "shm_array RESOLVED: " << name << std::endl;
      }
   });
   if (ref) {
      *this = ref;
#if 0
   } else {
      //std::cerr << "waiting for completion of " << name << std::endl;
      auto obj = ar.currentObject();
      if (obj && !valid())
         obj->unresolvedReference();
#endif
   }
}

template<class T>
T &shm_ref<T>::operator*() {
    return *m_p;
}

template<class T>
const T &shm_ref<T>::operator*() const {
    return *m_p;
}

template<class T>
T *shm_ref<T>::operator->() {
#ifdef NO_SHMEM
    return m_p;
#else
    return m_p.get();
#endif
}

template<class T>
const T *shm_ref<T>::operator->() const {
#ifdef NO_SHMEM
    return m_p;
#else
    return m_p.get();
#endif
}

template<class T>
const shm_name_t &shm_ref<T>::name() const {
    return m_name;
}

}
#endif
