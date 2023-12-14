//
// Created by simon on 08.05.21.
//

#ifndef ROBUSTOPTIMIZATION_HELPERS_H
#define ROBUSTOPTIMIZATION_HELPERS_H

#include "Logger.h"

#include <exception>
#include <utility>
#include <type_traits>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

namespace helpers {

struct MyException : public std::exception {
    explicit MyException(std::string message);

    const char* what() const throw() final;
private:
    const std::string _message;
};

void exception_throw(const std::string& msg, Logger & logger = global_logger);

void exception_check(const bool test, const std::string& msg, Logger & logger = global_logger);

void warning_throw(const std::string& msg, Logger & logger = global_logger);

void warning_check(const bool test, const std::string& msg, Logger & logger = global_logger);

template<class ObjT>
class SmartIndex;

template<class ObjT>
class IndexedObject {
public:
    using Index = SmartIndex<ObjT>;

    explicit IndexedObject(SmartIndex<ObjT> id);

    SmartIndex<ObjT> const& id() const;

private:
    SmartIndex<ObjT> const _id;
};

template<class ObjT>
class IndexedObjectOwner {
    static_assert(std::is_base_of<IndexedObject<ObjT>, ObjT>::value,
                  "IndexedObjectOwner can only be used on objects derived from IndexedObject!");
    using IdT = SmartIndex<ObjT>;
    friend SmartIndex<ObjT>;

public:
    IndexedObjectOwner() = default;

    std::size_t num_objects() const;

protected:
    template<class... Args>
    IdT base_add_object(Args&& ... args);

    std::vector<IdT> const& ids() const;

    std::vector<ObjT> const& objects() const;

    std::vector<ObjT> & non_const_objects();

    ObjT const& object(IdT id) const;

    ObjT & object(IdT id);

    void clear();

private:
    IdT next_id() const;

private:
    std::vector<ObjT> _objects;
    std::vector<IdT> _ids;
};

template<class ObjT>
class SmartIndex {
    friend IndexedObjectOwner<ObjT>;
public:
    using Owner = IndexedObjectOwner<ObjT>;

public:
    static SmartIndex NO_INDEX;

public:
    SmartIndex(size_t id, Owner const* owner);

    Owner const& owner() const;

    ObjT const* operator->() const;

    auto operator==(SmartIndex<ObjT> const& other) const;
    auto operator!=(SmartIndex<ObjT> const& other) const;
    auto operator<=>(SmartIndex<ObjT> const& other) const;

    explicit operator size_t() const;

    size_t raw_id() const;

private:
    size_t _id;
    IndexedObjectOwner<ObjT> const* _owner;
};

template<class ObjT>
SmartIndex<ObjT> SmartIndex<ObjT>::NO_INDEX = SmartIndex<ObjT>(0, nullptr);


template<class T>
struct VectorSlice {
    VectorSlice(std::vector<T> const& vector, size_t const begin_idx, size_t const end_idx) : _vector(&vector),
                                                                                              _begin_idx(begin_idx),
                                                                                              _end_idx(end_idx) {}

    T const& at(size_t i) const {
        exception_check(0 <= i and i <= _end_idx - _begin_idx - 1, "Out of bounds!");
        return _vector->at(i + _begin_idx);
    }

    T const& operator[](size_t i) const {
        return (*_vector)[i + _begin_idx];
    }

    typename std::vector<T>::const_iterator begin() const {
        return _vector->begin() + _begin_idx;
    }

    typename std::vector<T>::const_iterator end() const {
        return _vector->begin() + _end_idx;
    }


private:
    std::vector<T> const* _vector;
    size_t _begin_idx = 0;
    size_t _end_idx = 0;
};

double mean(std::vector<double> const& data);
double variance(std::vector<double> const& data);
double standard_deviation(std::vector<double> const& data);
}

#include "helpers.tplt"

#endif //ROBUSTOPTIMIZATION_HELPERS_H
