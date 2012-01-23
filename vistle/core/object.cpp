#include <iostream>
#include <iomanip>

#include <limits.h>

#include "message.h"
#include "messagequeue.h"
#include "object.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

using namespace boost::interprocess;

namespace vistle {

Shm* Shm::singleton = NULL;

Shm::Shm(const int m, const int r, const size_t & size,
         message::MessageQueue * mq)
   : moduleID(m), rank(r), objectID(0), messageQueue(mq) {

   shm = new managed_shared_memory(open_or_create, "vistle", size);
}

Shm::~Shm() {

   delete shm;
}

Shm & Shm::instance(const int moduleID, const int rank,
                    message::MessageQueue * mq) {

   if (!singleton)
      singleton = new Shm(moduleID, rank, INT_MAX, mq);

   return *singleton;
}

Shm & Shm::instance() {

   if (!singleton)
      singleton = new Shm(-1, -1, INT_MAX, NULL);//34359738368); // 32GB

   return *singleton;
}

managed_shared_memory & Shm::getShm() {

   return *shm;
}

void Shm::publish(const shm_handle_t & handle) {

   if (messageQueue) {
      vistle::message::NewObject n(moduleID, rank, handle);
      messageQueue->getMessageQueue().send(&n, sizeof(n), 0);
   }
}

std::string Shm::createObjectID() {

   std::stringstream name;
   name << "Object_" << std::setw(8) << std::setfill('0') << moduleID
        << "_" << std::setw(8) << std::setfill('0') << rank
        << "_" << std::setw(8) << std::setfill('0') << objectID ++;

   return name.str();
}

Object * Shm::getObjectFromHandle(const shm_handle_t & handle) {

   try {
      Object *object = static_cast<Object *>
         (shm->get_address_from_handle(handle));
      return object;

   } catch (interprocess_exception &ex) { }

   return NULL;
}

Object::Object(const Type type, const std::string & n): id(type) {

   size_t size = MIN(n.size(), 31);
   n.copy(name, size);
   name[size] = 0;
}

Object::~Object() {

}

Object::Type Object::getType() const {

   return id;
}

std::string Object::getName() const {

   return name;
}

Triangles::Triangles(const size_t & numCorners, const size_t & numVertices,
                     const std::string & name)
   : Object(Object::TRIANGLES, name) {

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   const allocator<size_t, managed_shared_memory::segment_manager>
      alloc_inst_size_t(Shm::instance().getShm().get_segment_manager());

   x = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);
   y = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);
   z = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   cl = Shm::instance().getShm().construct<std::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numCorners, size_t(), alloc_inst_size_t);
}


Triangles * Triangles::create(const size_t & numCorners,
                              const size_t & numVertices) {

   const std::string name = Shm::instance().createObjectID();
   Triangles *t = static_cast<Triangles *>
      (Shm::instance().getShm().construct<Triangles>(name.c_str())[1](numCorners, numVertices, name));

   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(t);

   Shm::instance().publish(handle);
   return t;
}

size_t Triangles::getNumCorners() const {

   return cl->size();
}

size_t Triangles::getNumVertices() const {

   return x->size();
}

Lines::Lines(const size_t & numElements, const size_t & numCorners,
             const size_t & numVertices, const std::string & name)
   : Object(Object::LINES, name) {

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   const allocator<size_t, managed_shared_memory::segment_manager>
      alloc_inst_size_t(Shm::instance().getShm().get_segment_manager());

   x = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   y = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   z = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   el = Shm::instance().getShm().construct<std::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, size_t(), alloc_inst_size_t);

   cl = Shm::instance().getShm().construct<std::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, size_t(), alloc_inst_size_t);
}


Lines * Lines::create(const size_t & numElements, const size_t & numCorners,
                      const size_t & numVertices) {

   const std::string name = Shm::instance().createObjectID();
   Lines *l = static_cast<Lines *>
      (Shm::instance().getShm().construct<Lines>(name.c_str())[1](numElements, numCorners, numVertices, name));

   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(l);

   Shm::instance().publish(handle);
   return l;
}

size_t Lines::getNumElements() const {

   return el->size();
}

size_t Lines::getNumCorners() const {

   return cl->size();
}

size_t Lines::getNumVertices() const {

   return x->size();
}


UnstructuredGrid::UnstructuredGrid(const size_t & numElements,
                                   const size_t & numCorners,
                                   const size_t & numVertices,
                                   const std::string & name)
   : Object(Object::UNSTRUCTUREDGRID, name) {

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   const allocator<size_t, managed_shared_memory::segment_manager>
      alloc_inst_size_t(Shm::instance().getShm().get_segment_manager());

   x = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   y = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   z = Shm::instance().getShm().construct<std::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   tl = Shm::instance().getShm().construct<std::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, size_t(), alloc_inst_size_t);

   el = Shm::instance().getShm().construct<std::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, size_t(), alloc_inst_size_t);

   cl = Shm::instance().getShm().construct<std::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numCorners, size_t(), alloc_inst_size_t);
}

UnstructuredGrid * UnstructuredGrid::create(const size_t & numElements,
                                            const size_t & numCorners,
                                            const size_t & numVertices) {

   const std::string name = Shm::instance().createObjectID();
   UnstructuredGrid *u = static_cast<UnstructuredGrid *>
      (Shm::instance().getShm().construct<UnstructuredGrid>(name.c_str())[1](numElements, numCorners, numVertices, name));

   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(u);

   Shm::instance().publish(handle);
   return u;
}

size_t UnstructuredGrid::getNumElements() const {

   return el->size();
}

size_t UnstructuredGrid::getNumCorners() const {

   return cl->size();
}

size_t UnstructuredGrid::getNumVertices() const {

   return x->size();
}

template<> const Object::Type Vec<float>::type  = Object::VECFLOAT;
template<> const Object::Type Vec<int>::type    = Object::VECINT;
template<> const Object::Type Vec<char>::type   = Object::VECCHAR;
template<> const Object::Type Vec3<float>::type = Object::VEC3FLOAT;
template<> const Object::Type Vec3<int>::type   = Object::VEC3INT;
template<> const Object::Type Vec3<char>::type  = Object::VEC3CHAR;

} // namespace vistle
