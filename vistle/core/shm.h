#ifndef SHM_H
#define SHM_H

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>

#include "exception.h"
#include "index.h"
#include "export.h"

//#define SHMDEBUG
//#define SHMPUBLISH

namespace vistle {

typedef boost::interprocess::managed_shared_memory::handle_t shm_handle_t;

struct V_COREEXPORT shm_name_t {
   char name[32];
   shm_name_t(const std::string &s = "INVALID");

   operator const char *() const;
   operator char *();
   operator std::string () const;
};
std::string operator+(const std::string &s, const shm_name_t &n);

namespace message {
   class MessageQueue;
}

class V_COREEXPORT shm_exception: public exception {
   public:
   shm_exception(const std::string &what = "shared memory error") : exception(what) {}
};

class Object;

#ifdef SHMDEBUG
struct ShmDebugInfo {
   shm_name_t name;
   shm_handle_t handle;
   char deleted;
   char type;

   ShmDebugInfo(char type='\0', const std::string &name = "", shm_handle_t handle = 0)
      : handle(handle)
      , deleted(0)
      , type(type)
   {
      memset(this->name, '\0', sizeof(this->name));
      strncpy(this->name, name.c_str(), sizeof(this->name)-1);
   }
};
#endif

template<typename T, class allocator>
class shm_array {

 public: 
   typedef T value_type;
   shm_array(const allocator &alloc = allocator()) : m_allocator(alloc) {}
   shm_array(const size_t size, const allocator &alloc = allocator()) : m_allocator(alloc) { resize(size); }
   shm_array(const size_t size, const T &value, const allocator &alloc = allocator()) : m_allocator(alloc) { resize(size, value); }
   template< class InputIt >
      shm_array( InputIt first, InputIt last, 
            const allocator &alloc = allocator() );
   shm_array(const shm_array &other);
   shm_array &operator=(const shm_array &rhs);
   ~shm_array() { m_allocator.deallocate(m_data, m_capacity); }

   typedef typename allocator::pointer pointer;
   typedef T *iterator;
   typedef const T *const_iterator;

   iterator begin() const { return &*m_data; }
   iterator end() const { return (&*m_data) + m_size; }
   T *data() const { return &*m_data; }

   T &operator[](const size_t idx) { return m_data[idx]; }
   T &operator[](const size_t idx) const { return m_data[idx]; }
   void push_back(const T &v) { if (m_size >= m_capacity) reserve(m_capacity==0 ? 1 : m_capacity*2); assert(m_size < m_capacity); m_data[m_size] = v; ++m_size; }
   T &back() { return m_data[m_size-1]; }
   T &front() { return m_data[0]; }

   bool empty() const { return m_size == 0; }
   void clear() { resize(0); }
   size_t size() const { return m_size; }
   void resize(const size_t size) { reserve(size); m_size = size; }
   void resize(const size_t size, const T &value) { reserve(size); for (size_t i=m_size; i<size; ++i) m_data[i] = value; m_size = size; }
   size_t capacity() const { return m_capacity; }
   void reserve(const size_t size) { pointer new_data = m_allocator.allocate(size); if (m_data) ::memcpy(&*new_data, &*m_data, m_size*sizeof(T)); m_allocator.deallocate(m_data, m_capacity); m_data = new_data; m_capacity = size; }

 private:
   size_t m_size = 0;
   size_t m_capacity = 0;
   pointer m_data = NULL;
   allocator m_allocator;

   friend class boost::serialization::access;
   template<class Archive>
      void serialize(Archive &ar, const unsigned int version);
};

template<typename T>
struct shm {
   typedef boost::interprocess::allocator<T, boost::interprocess::managed_shared_memory::segment_manager> allocator;
   typedef boost::interprocess::basic_string<T, std::char_traits<T>, allocator> string;
   typedef boost::interprocess::vector<T, allocator> vector;
   typedef boost::interprocess::offset_ptr<vector> ptr;
   typedef vistle::shm_array<T, allocator> array;
   typedef boost::interprocess::offset_ptr<array> array_ptr;
   static typename boost::interprocess::managed_shared_memory::segment_manager::template construct_proxy<T>::type construct(const std::string &name);
   static T *find(const std::string &name);
   static void destroy(const std::string &name);
};

class V_COREEXPORT Shm {

 public:
   static Shm & the();
   static Shm & create(const std::string &shmname, const int moduleID, const int rank,
                         message::MessageQueue *messageQueue = NULL);
   static Shm & attach(const std::string &shmname, const int moduleID, const int rank,
                         message::MessageQueue *messageQueue = NULL);
   void detach();

   const std::string &name() const;

   typedef boost::interprocess::allocator<void, boost::interprocess::managed_shared_memory::segment_manager> void_allocator;
   const void_allocator &allocator() const;

   boost::interprocess::managed_shared_memory &shm();
   const boost::interprocess::managed_shared_memory &shm() const;
   std::string createObjectID();

   boost::shared_ptr<const Object> getObjectFromHandle(const shm_handle_t & handle) const;
   shm_handle_t getHandleFromObject(boost::shared_ptr<const Object> object) const;
   shm_handle_t getHandleFromObject(const Object *object) const;
   boost::shared_ptr<const Object> getObjectFromName(const std::string &name) const;

   static std::string shmIdFilename();
   static bool cleanAll();

#ifdef SHMDEBUG
   static vistle::shm<ShmDebugInfo>::vector *s_shmdebug;
   void markAsRemoved(const std::string &name);
#endif

 private:
   Shm(const std::string &name, const int moduleID, const int rank, const size_t size,
       message::MessageQueue *messageQueue, bool create);
   ~Shm();

   void_allocator *m_allocator;
   std::string m_name;
   bool m_created;
   const int m_moduleID;
   const int m_rank;
   int m_objectID;
   static Shm *s_singleton;
   boost::interprocess::managed_shared_memory *m_shm;
};

template<typename T>
typename boost::interprocess::managed_shared_memory::segment_manager::template construct_proxy<T>::type shm<T>::construct(const std::string &name) {
   return Shm::the().shm().construct<T>(name.c_str());
}

template<typename T>
T *shm<T>::find(const std::string &name) {
   return Shm::the().shm().find<T>(name.c_str()).first;
}

template<typename T>
void shm<T>::destroy(const std::string &name) {
      Shm::the().shm().destroy<T>(name.c_str());
#ifdef SHMDEBUG
      Shm::the().markAsRemoved(name);
#endif
}

template<typename T>
class V_COREEXPORT ShmVector {
#ifdef SHMDEBUG
   friend void Shm::markAsRemoved(const std::string &name);
#endif

   public:
      class ptr {
         public:
            ptr(ShmVector *p = NULL);
            ptr(const ptr &ptr);
            ~ptr();
            ptr &operator=(const ptr &other);
            ptr &operator=(ShmVector *p);

            ShmVector &operator*() {
               return *m_p;
            }
            ShmVector *operator->() {
               return &*m_p;
            }

         private:
            boost::interprocess::offset_ptr<ShmVector> m_p;
      };

      ShmVector(Index size = 0);
      int refcount() const;
      void* operator new(size_t size);
      void operator delete(void *p);

      T &operator[](Index i) { return (*m_x)[i]; }
      const T &operator[](Index i) const { return (*m_x)[i]; }

      Index size() const { return m_x->size(); }
      void resize(Index s);

      typename shm<T>::array_ptr &operator()() { return m_x; }
      typename shm<const T>::array_ptr &operator()() const { return m_x; }

      void push_back(const T &d) { m_x->push_back(d); }

   private:
      ~ShmVector();
      void ref();
      void unref();

      friend class boost::serialization::access;
      template<class Archive>
         void serialize(Archive &ar, const unsigned int version);
      template<class Archive>
         void save(Archive &ar, const unsigned int version) const;
      template<class Archive>
         void load(Archive &ar, const unsigned int version);

      boost::interprocess::interprocess_mutex m_mutex;
      int m_refcount;
      shm_name_t m_name;
      typename shm<T>::array_ptr m_x;
};

} // namespace vistle

#endif

#ifdef VISTLE_IMPL
#include "shm_impl.h"
#endif
