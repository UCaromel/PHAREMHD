//#include "python3/pybind_def.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "dict.hpp"

namespace py = pybind11;

//using PHARE::initializer::InitFunction;


using Dict_t = cppdict::Dict<int,  double>;

template<typename T>
void add(std::string const& path, T&& value)
{
    //cppdict::add(path, std::forward<T>(value),
        //         PHARE::initializer::PHAREDictHandler::INSTANCE().dict());
        Dict_t d;
}

PYBIND11_MODULE(pyMHD, m)
{
    // expose dict add function per template type to force casts/error when used from python
    
     m.def("add_size_t", add<std::size_t>, "add_size_t");

    // m.def("addInitFunction1D", add<InitFunction<1>>, "add");
    // m.def("addInitFunction2D", add<InitFunction<2>>, "add");
    // m.def("addInitFunction3D", add<InitFunction<3>>, "add");

    //m.def("stop", []() { PHARE::initializer::PHAREDictHandler::INSTANCE().stop(); });
}