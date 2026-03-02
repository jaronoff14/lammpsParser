#include <pybind11/pybind11.h>
#include <pybind11/stl.h>   // automatic conversion for STL containers
#include "lammps_parser.h"

namespace py = pybind11;

PYBIND11_MODULE(lammpsParser, m) {
    m.doc() = "LAMMPS dump parser (header + flexible atoms format)";

    // BoxBound
    py::class_<lammps_parser::BoxBound>(m, "BoxBound")
        .def_readwrite("lo", &lammps_parser::BoxBound::lo)
        .def_readwrite("hi", &lammps_parser::BoxBound::hi)
        .def_readwrite("tilt", &lammps_parser::BoxBound::tilt);

    // Header
    py::class_<lammps_parser::Header>(m, "Header")
        .def_readwrite("timestep", &lammps_parser::Header::timestep)
        .def_readwrite("n_atoms", &lammps_parser::Header::n_atoms)
        .def_readwrite("box_bounds", &lammps_parser::Header::box_bounds)
        .def_readwrite("box_bounds_flags", &lammps_parser::Header::box_bounds_flags);

    // Atom
    py::class_<lammps_parser::Atom>(m, "Atom")
        .def_readwrite("id", &lammps_parser::Atom::id)
        .def_readwrite("type", &lammps_parser::Atom::type)
        .def_readwrite("xu", &lammps_parser::Atom::xu)
        .def_readwrite("yu", &lammps_parser::Atom::yu)
        .def_readwrite("zu", &lammps_parser::Atom::zu)
        .def_readwrite("mol", &lammps_parser::Atom::mol)
        .def_readwrite("c_KE", &lammps_parser::Atom::c_KE)
        .def_readwrite("c_PE", &lammps_parser::Atom::c_PE)
        .def_readwrite("c_stress", &lammps_parser::Atom::c_stress)
        .def_readwrite("raw_tokens", &lammps_parser::Atom::raw_tokens)
        .def_readwrite("extras", &lammps_parser::Atom::extras);

    // ParsedFrame
    py::class_<lammps_parser::ParsedFrame>(m, "ParsedFrame")
        .def_readwrite("header", &lammps_parser::ParsedFrame::header)
        .def_readwrite("atoms", &lammps_parser::ParsedFrame::atoms);

    // Functions
    m.def("parse_frame_from_file",
          &lammps_parser::parse_frame_from_file,
          py::arg("filepath"),
          py::arg("format_override") = std::string(lammps_parser::DEFAULT_FORMAT),
          "Parse a single file and return a ParsedFrame.");

    m.def("parse_frames_from_directory",
          &lammps_parser::parse_frames_from_directory,
          py::arg("dirpath"),
          py::arg("format_override") = std::string(lammps_parser::DEFAULT_FORMAT),
          "Parse all files matching `dump.*.lammpstrj` in a directory and return a list of ParsedFrame.");
}
