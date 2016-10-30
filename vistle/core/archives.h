#ifndef ARCHIVES_H
#define ARCHIVES_H

#include <functional>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/archive/detail/oserializer.hpp>
#include <boost/archive/detail/iserializer.hpp>

#include <boost/mpl/vector.hpp>

#include <util/vecstreambuf.h>
#include "shm.h"
#include "findobjectreferenceoarchive.h"

namespace vistle {

template<class T>
class shm_obj_ref;

class shallow_oarchive;
class shallow_iarchive;

}

namespace boost {
namespace archive {
extern template class V_COREEXPORT basic_binary_oprimitive<
    vistle::shallow_oarchive,
    std::ostream::char_type, 
    std::ostream::traits_type
>;
extern template class V_COREEXPORT basic_binary_iprimitive<
    vistle::shallow_iarchive,
    std::istream::char_type, 
    std::istream::traits_type
>;
} // namespace archive
} // namespace boost

namespace vistle {

class Object;
struct ObjectData;
typedef std::shared_ptr<Object> obj_ptr;
typedef std::shared_ptr<const Object> obj_const_ptr;

class V_COREEXPORT shallow_oarchive: public boost::archive::binary_oarchive_impl<shallow_oarchive, std::ostream::char_type, std::ostream::traits_type> {

    typedef boost::archive::binary_oarchive_impl<shallow_oarchive, std::ostream::char_type, std::ostream::traits_type> Base;
public:
    shallow_oarchive(std::ostream &os, unsigned int flags=0);
    shallow_oarchive(std::streambuf &bsb, unsigned int flags=0);
    ~shallow_oarchive();

};

class V_COREEXPORT deep_oarchive: public shallow_oarchive {

    typedef shallow_oarchive Base;
    void move_subarchives(deep_oarchive &other) {
        for (auto &p: other.objects)
            objects.insert(std::move(p));
        other.objects.clear();
        for (auto &p: other.arrays)
            arrays.insert(std::move(p));
        other.arrays.clear();
    }

public:
    deep_oarchive(std::ostream &os, unsigned int flags=0);
    deep_oarchive(std::streambuf &bsb, unsigned int flags=0);
    ~deep_oarchive();
    std::map<std::string, std::vector<char>> objects;
    std::map<std::string, std::vector<char>> arrays;

    template<class T>
    deep_oarchive &operator<<(vistle::ShmVector<T> &t) {
        vecostreambuf<char> vb;
        deep_oarchive ar(vb);
        ar & *t;
        ar.arrays.emplace(t->name(), vb.get_vector());
        move_subarchives(ar);
        return *this;
    }

    template<class T>
    deep_oarchive &operator<<(vistle::shm_obj_ref<T> &t) {
        vecostreambuf<char> vb;
        deep_oarchive ar(vb);
        ar & *t;
        ar.arrays.emplace(t->name(), vb.get_vector());
        move_subarchives(ar);
        return *this;
    }

    struct directory_entry {
        std::string name;
        bool is_array;
        size_t size;
        char *data;

        directory_entry(): is_array(false), size(0), data(nullptr) {}
        directory_entry(const std::string &name, bool is_array, size_t size, char *data)
            : name(name), is_array(is_array), size(size), data(data) {}

        template<class Archive>
        void serialize(Archive &ar, const unsigned int version) {
            ar & name;
            ar & is_array;
            ar & size;
        }
    };

    typedef std::vector<directory_entry> directory;
    directory get_directory() {
        directory dir;
        for (auto &obj: objects) {
            dir.emplace_back(obj.first, false, obj.second.size(), obj.second.data());
        }
        for (auto &arr: arrays) {
            dir.emplace_back(arr.first, true, arr.second.size(), arr.second.data());
        }
        return dir;
    }
};


class V_COREEXPORT Fetcher {
public:
    virtual ~Fetcher();
    virtual void requestArray(const std::string &name, int type, const std::function<void()> &completeCallback) = 0;
    virtual void requestObject(const std::string &name, const std::function<void()> &completeCallback) = 0;
};

class V_COREEXPORT shallow_iarchive: public boost::archive::binary_iarchive_impl<shallow_iarchive, std::istream::char_type, std::istream::traits_type> {

    typedef boost::archive::binary_iarchive_impl<shallow_iarchive, std::istream::char_type, std::istream::traits_type> Base;
public:
    shallow_iarchive(std::istream &is, unsigned int flags=0);
    shallow_iarchive(std::streambuf &bsb, unsigned int flags=0);
    ~shallow_iarchive();

    void setFetcher(std::shared_ptr<Fetcher> fetcher);
    void setCurrentObject(ObjectData *data);
    ObjectData *currentObject() const;

    template<typename T>
    ShmVector<T> getArray(const std::string &name, const std::function<void()> &completeCallback) const {
        auto arr = Shm::the().getArrayFromName<T>(name);
        if (!arr) {
            assert(m_fetcher);
            m_fetcher->requestArray(name, shm<T>::array::typeId(), completeCallback);
            arr = Shm::the().getArrayFromName<T>(name);
        }
        return arr;
    }

    obj_const_ptr getObject(const std::string &name, const std::function<void()> &completeCallback) const {
        auto obj = Shm::the().getObjectFromName(name);
        if (!obj) {
            assert(m_fetcher);
            m_fetcher->requestObject(name, completeCallback);
            obj = Shm::the().getObjectFromName(name);
        }
        return obj;
    }

    void setObjectCompletionHandler(const std::function<void()> &completer);
    const std::function<void()> &objectCompletionHandler() const;

private:
    std::shared_ptr<Fetcher> m_fetcher;
    ObjectData *m_currentObject;
    std::function<void()> m_completer;
};


class V_COREEXPORT deep_iarchive: public shallow_iarchive {

    typedef shallow_iarchive Base;
public:
    deep_iarchive(std::istream &is, unsigned int flags=0);
    deep_iarchive(std::streambuf &bsb, unsigned int flags=0);
    ~deep_iarchive();

    template<class T>
    deep_iarchive &operator>>(vistle::ShmVector<T> &t) {
        *this >> *t;
        return *this;
    }

    template<class T>
    deep_iarchive &operator>>(vistle::shm_obj_ref<T> &t) {
        *this >> *t;
        return *this;
    }
};



typedef boost::mpl::vector<
   shallow_iarchive,
   deep_iarchive
      > InputArchives;

typedef boost::mpl::vector<
   shallow_oarchive,
   deep_oarchive
      > OutputArchives;

} // namespace vistle


BOOST_SERIALIZATION_REGISTER_ARCHIVE(vistle::shallow_oarchive)
BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(vistle::shallow_oarchive)
BOOST_SERIALIZATION_REGISTER_ARCHIVE(vistle::deep_oarchive)
BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(vistle::deep_oarchive)
BOOST_SERIALIZATION_REGISTER_ARCHIVE(vistle::shallow_iarchive)
BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(vistle::shallow_iarchive)
BOOST_SERIALIZATION_REGISTER_ARCHIVE(vistle::deep_iarchive)
BOOST_SERIALIZATION_USE_ARRAY_OPTIMIZATION(vistle::deep_iarchive)

#endif
