#include "ConservativeVariablesCC.hpp"
#include "PrimitiveVariablesCC.hpp"
#include "PhareMHD.hpp"
#include "Enums.hpp"
#include "ModularityUtils.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

//#include "dict.hpp"

namespace py = pybind11;

py::array_t<double> convert_to_numpy(const std::vector<std::vector<double>>& vec) {
    size_t ny = vec.size();
    size_t nx = (ny > 0) ? vec[0].size() : 0;
    py::array_t<double> numpy_array({ny, nx});
    auto buffer = numpy_array.request();
    double* ptr = static_cast<double*>(buffer.ptr);
    for (size_t i = 0; i < ny; ++i) {
        if (vec[i].size() != nx) {
            throw std::runtime_error("All nested vectors must have the same size");
        }
        std::copy(vec[i].begin(), vec[i].end(), ptr + i * nx);
    }
    return numpy_array;
}

std::vector<std::vector<double>> convert_from_numpy(py::array_t<double> numpy_array) {
    auto buffer = numpy_array.request();
    double* ptr = static_cast<double*>(buffer.ptr);
    size_t ny = buffer.shape[0];
    size_t nx = buffer.shape[1];
    std::vector<std::vector<double>> vec(ny, std::vector<double>(nx));
    for (size_t i = 0; i < ny; ++i) {
        std::copy(ptr + i * nx, ptr + (i + 1) * nx, vec[i].begin());
    }
    return vec;
}

PYBIND11_MODULE(pyMHD, m)
{
    py::class_<PrimitiveVariablesCC>(m, "PrimitiveVariablesCC")

        .def_readwrite("nx", &PrimitiveVariablesCC::nx)
        .def_readwrite("ny", &PrimitiveVariablesCC::ny)
        .def_property("rho",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.rho); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.rho = convert_from_numpy(array); })
        .def_property("vx",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.vx); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.vx = convert_from_numpy(array); })
        .def_property("vy",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.vy); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.vy = convert_from_numpy(array); })
        .def_property("vz",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.vz); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.vz = convert_from_numpy(array); })
        .def_property("Bx",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.Bx); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.Bx = convert_from_numpy(array); })
        .def_property("By",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.By); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.By = convert_from_numpy(array); })
        .def_property("Bz",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.Bz); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.Bz = convert_from_numpy(array); })
        .def_property("P",
            [](PrimitiveVariablesCC& self) { return convert_to_numpy(self.P); },
            [](PrimitiveVariablesCC& self, py::array_t<double> array) { self.P = convert_from_numpy(array); })


        .def(py::init<double,double>())
        .def(py::init<const ConservativeVariablesCC&>())
        .def("set", &PrimitiveVariablesCC::set)
        .def("init", &PrimitiveVariablesCC::init)
        .def("__call__", &PrimitiveVariablesCC::operator())
        .def("assign", [](PrimitiveVariablesCC& self, const PrimitiveVariablesCC& other) {
            self = other;
        });

    py::enum_<BoundaryConditions>(m, "BoundaryConditions")
        .value("Periodic", Periodic)
        .value("ZeroGradient", ZeroGradient);

    py::enum_<Reconstruction>(m, "Reconstruction")
        .value("Constant", Constant)
        .value("Linear", Linear);
    
    py::enum_<Slope>(m, "Slope")
        .value("VanLeer", VanLeer)
        .value("MinMod", MinMod);

    py::enum_<Riemann>(m, "RiemannSolver")
        .value("Rusanov", Rusanov)
        .value("HLL", HLL);

    py::enum_<CTMethod>(m, "CTMethod")
        .value("Average", Average)
        .value("UCT_HLL", UCT_HLL);

    py::enum_<Integrator>(m, "Integrator")
        .value("EulerIntegrator", EulerIntegrator)
        .value("TVDRK2Integrator", TVDRK2Integrator)
        .value("TVDRK3Integrator", TVDRK3Integrator);
    
    py::enum_<dumpVariables>(m, "dumpVariables")
        .value("Primitive", Primitive)
        .value("Conservative", Conservative)
        .value("Both", Both);
    
    m.def("PhareMHD", &PhareMHD, 
          py::arg("primvar0"), py::arg("resultdir"), py::arg("order"), py::arg("nghost"), 
          py::arg("boundaryconditions"), py::arg("reconstruction"), py::arg("slopelimiter"), py::arg("riemannsolver"), py::arg("constainedtransport"), py::arg("timeintegrator"), 
          py::arg("Dx"), py::arg("Dy"), py::arg("FinalTime"), py::arg("Dt") = 0.0, py::arg("dumpvariables") = Conservative, py::arg("dumpfrequency") = 1);

/*
    py::class_<ConservativeVariablesCC>(m, "ConservativeVariablesCC")

        .def_readwrite("nx", &ConservativeVariablesCC::nx)
        .def_readwrite("ny", &ConservativeVariablesCC::ny)
        .def_property("rho",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.rho); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.rho = convert_from_numpy(array); })
        .def_property("rhovx",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.rhovx); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.rhovx = convert_from_numpy(array); })
        .def_property("rhovy",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.rhovy); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.rhovy = convert_from_numpy(array); })
        .def_property("rhovz",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.rhovz); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.rhovz = convert_from_numpy(array); })
        .def_property("Bx",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.Bx); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.Bx = convert_from_numpy(array); })
        .def_property("By",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.By); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.By = convert_from_numpy(array); })
        .def_property("Bz",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.Bz); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.Bz = convert_from_numpy(array); })
        .def_property("Etot",
            [](ConservativeVariablesCC& self) { return convert_to_numpy(self.Etot); },
            [](ConservativeVariablesCC& self, py::array_t<double> array) { self.Etot = convert_from_numpy(array); })


        .def(py::init<double,double>())
        .def(py::init<const PrimitiveVariablesCC&>())
        .def("set", &ConservativeVariablesCC::set)
        .def("__call__", &ConservativeVariablesCC::operator())
        .def(py::self * double())
        .def(py::self - py::self)
        .def(py::self + py::self)
        .def("assign", [](ConservativeVariablesCC& self, const ConservativeVariablesCC& other) {
            self = other;
        });
 
    py::class_<ReconstructedValues>(m, "ReconstructedValues")

        .def_readwrite("rho", &ReconstructedValues::rho)
        .def_readwrite("vx", &ReconstructedValues::vx)
        .def_readwrite("vy", &ReconstructedValues::vy)
        .def_readwrite("vz", &ReconstructedValues::vz)
        .def_readwrite("Bx", &ReconstructedValues::Bx)
        .def_readwrite("By", &ReconstructedValues::By)
        .def_readwrite("Bz", &ReconstructedValues::Bz)
        .def_readwrite("P", &ReconstructedValues::P)

        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def("assign", [](ReconstructedValues& self, const ReconstructedValues& other) {
            self = other;
        });

    
    m.def("InitialiseGhostCells", &InitialiseGhostCells<PrimitiveVariablesCC>);
    m.def("UpdateGhostCells", &UpdateGhostCells<PrimitiveVariablesCC>);
    m.def("UpdateGhostCells", &UpdateGhostCells<ConservativeVariablesCC>);

    m.attr("gam") = &gam;
    m.def("EosEtot", &EosEtot);
    m.def("EosP", &EosP);

    py::enum_<Dir>(m, "Dir")
        .value("X", Dir::X)
        .value("Y", Dir::Y);

    m.def("ComputeFluxVector", &ComputeFluxVector);

    py::class_<Interface>(m, "Interface")

        .def_readwrite("uL", &Interface::uL)
        .def_readwrite("uR", &Interface::uR)
        .def_readwrite("fL", &Interface::fL)
        .def_readwrite("fR", &Interface::fR)
        .def_readwrite("SL", &Interface::SL)
        .def_readwrite("SR", &Interface::SR)
        .def_readwrite("Splus", &Interface::Splus)

        .def(py::init<const PrimitiveVariablesCC&, int, int, int, int, Dir>());
    

    m.def("RusanovRiemannSolver", &RusanovRiemannSolver);

    m.def("GodunovFluxX", &GodunovFluxX);
    m.def("GodunovFluxY", &GodunovFluxY);
    m.def("ComputeFluxDifferenceX", &ComputeFluxDifferenceX);
    m.def("ComputeFluxDifferenceY", &ComputeFluxDifferenceY);

    py::class_<ConstrainedTransport>(m, "ConstrainedTransport")

        .def_readwrite("vx", &ConstrainedTransport::vx)
        .def_readwrite("vy", &ConstrainedTransport::vy)
        .def_readwrite("Bx", &ConstrainedTransport::Bx)
        .def_readwrite("By", &ConstrainedTransport::By)
        .def_readwrite("Ez", &ConstrainedTransport::Ez)
        .def_readwrite("bx", &ConstrainedTransport::bx)
        .def_readwrite("by", &ConstrainedTransport::by)
        .def_readwrite("BX", &ConstrainedTransport::BX)
        .def_readwrite("BY", &ConstrainedTransport::BY)

        .def(py::init<const ConservativeVariablesCC&, double, double, double, int>());
    
    m.def("ApplyConstrainedTransport", &ApplyConstrainedTransport);

    m.def("EulerAdvance", &EulerAdvance);
    m.def("Euler", &Euler);
    m.def("TVDRK2", &TVDRK2);
    m.def("TVDRK3", &TVDRK3);

    m.def("ComPuteNewDt", &ComPuteNewDt);
*/
}

/*
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
    
    //m.def("add_size_t", add<std::size_t>, "add_size_t");

    // m.def("addInitFunction1D", add<InitFunction<1>>, "add");
    // m.def("addInitFunction2D", add<InitFunction<2>>, "add");
    // m.def("addInitFunction3D", add<InitFunction<3>>, "add");

    //m.def("stop", []() { PHARE::initializer::PHAREDictHandler::INSTANCE().stop(); });

}*/