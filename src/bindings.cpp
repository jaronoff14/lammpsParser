#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "lammps_parser.h"

#include <filesystem>
#include <string>
#include <set>
#include <vector>

namespace py = pybind11;
namespace fs = std::filesystem;

PYBIND11_MODULE(lammpsParser, m) {
    m.doc() = "LAMMPS dump parser (header + flexible atoms format) with NumPy column API and installer";

    py::class_<lammps_parser::BoxBound>(m, "BoxBound")
        .def_readwrite("lo", &lammps_parser::BoxBound::lo)
        .def_readwrite("hi", &lammps_parser::BoxBound::hi)
        .def_readwrite("tilt", &lammps_parser::BoxBound::tilt);

    py::class_<lammps_parser::Header>(m, "Header")
        .def_readwrite("timestep", &lammps_parser::Header::timestep)
        .def_readwrite("n_atoms", &lammps_parser::Header::n_atoms)
        .def_readwrite("box_bounds", &lammps_parser::Header::box_bounds)
        .def_readwrite("box_bounds_flags", &lammps_parser::Header::box_bounds_flags);

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
        .def_readwrite("positions", &lammps_parser::Atom::positions)
        .def_readwrite("raw_tokens", &lammps_parser::Atom::raw_tokens)
        .def_readwrite("extras", &lammps_parser::Atom::extras);

    py::class_<lammps_parser::ParsedFrame>(m, "ParsedFrame")
        .def_readwrite("header", &lammps_parser::ParsedFrame::header)
        .def_readwrite("atoms", &lammps_parser::ParsedFrame::atoms);

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

    // parse_frame_to_numpy (fixed array constructors)
    m.def("parse_frame_to_numpy",
          [](const std::string &filepath, const std::string &format_override) {
              lammps_parser::ParsedFrame pf = lammps_parser::parse_frame_from_file(filepath, format_override);
              size_t N = pf.atoms.size();
              ssize_t Ns = static_cast<ssize_t>(N);

              // 1-D int64 arrays: use array_t(ssize_t) constructor
              py::array_t<long long> arr_id(Ns);
              py::array_t<long long> arr_type(Ns);
              py::array_t<long long> arr_mol(Ns);

              auto buf_id = arr_id.mutable_unchecked<1>();
              auto buf_type = arr_type.mutable_unchecked<1>();
              auto buf_mol = arr_mol.mutable_unchecked<1>();

              // 1-D double arrays: use array_t(ssize_t)
              py::array_t<double> arr_xu(Ns);
              py::array_t<double> arr_yu(Ns);
              py::array_t<double> arr_zu(Ns);
              py::array_t<double> arr_cKE(Ns);
              py::array_t<double> arr_cPE(Ns);

              auto buf_xu = arr_xu.mutable_unchecked<1>();
              auto buf_yu = arr_yu.mutable_unchecked<1>();
              auto buf_zu = arr_zu.mutable_unchecked<1>();
              auto buf_cKE = arr_cKE.mutable_unchecked<1>();
              auto buf_cPE = arr_cPE.mutable_unchecked<1>();

              // c_stress as (N,6) array: construct with a shape vector (portable)
              std::vector<ssize_t> shape2 = { Ns, static_cast<ssize_t>(6) };
              py::array_t<double> arr_cstress(shape2);
              auto buf_cstress = arr_cstress.mutable_unchecked<2>();

              // Collect extras keys across atoms
              std::set<std::string> extras_keys_set;
              for (const auto &a : pf.atoms) {
                  for (const auto &kv : a.extras) extras_keys_set.insert(kv.first);
              }
              std::vector<std::string> extras_keys(extras_keys_set.begin(), extras_keys_set.end());
              // Prepare extras arrays (1-D)
              std::vector<py::array_t<double>> extras_arrays;
              extras_arrays.reserve(extras_keys.size());
              for (size_t k = 0; k < extras_keys.size(); ++k) extras_arrays.emplace_back(py::array_t<double>(Ns));

              // Fill arrays
              for (size_t i = 0; i < N; ++i) {
                  const auto &a = pf.atoms[i];
                  ssize_t idx = static_cast<ssize_t>(i);
                  buf_id(idx) = a.id;
                  buf_type(idx) = a.type;
                  buf_mol(idx) = a.mol;
                  buf_xu(idx) = a.xu;
                  buf_yu(idx) = a.yu;
                  buf_zu(idx) = a.zu;
                  buf_cKE(idx) = a.c_KE;
                  buf_cPE(idx) = a.c_PE;
                  for (int s = 0; s < 6; ++s) buf_cstress(idx, s) = a.c_stress[s];
                  // extras
                  for (size_t k = 0; k < extras_keys.size(); ++k) {
                      const std::string &key = extras_keys[k];
                      double val = 0.0;
                      auto it = a.extras.find(key);
                      if (it != a.extras.end()) val = it->second;
                      auto mut = extras_arrays[k].mutable_unchecked<1>();
                      mut(idx) = val;
                  }
              }

              // Build Python dict result
              py::dict result;
              py::dict header;
              header["timestep"] = pf.header.timestep;
              header["n_atoms"] = (long long)pf.header.n_atoms;
              header["box_bounds_flags"] = pf.header.box_bounds_flags;
              // box_bounds as list of (lo, hi, tilt)
              py::list bb_list;
              for (int i = 0; i < 3; ++i) {
                  py::tuple t = py::make_tuple(pf.header.box_bounds[i].lo,
                                               pf.header.box_bounds[i].hi,
                                               pf.header.box_bounds[i].tilt);
                  bb_list.append(t);
              }
              header["box_bounds"] = bb_list;

              result["header"] = header;
              result["n_atoms"] = (long long)N;

              py::dict columns;
              columns["id"] = arr_id;
              columns["type"] = arr_type;
              columns["mol"] = arr_mol;
              columns["xu"] = arr_xu;
              columns["yu"] = arr_yu;
              columns["zu"] = arr_zu;
              columns["c_KE"] = arr_cKE;
              columns["c_PE"] = arr_cPE;
              columns["c_stress"] = arr_cstress;

              for (size_t k = 0; k < extras_keys.size(); ++k) {
                  columns[extras_keys[k].c_str()] = extras_arrays[k];
              }

              result["columns"] = columns;
              return result;
          },
          py::arg("filepath"),
          py::arg("format_override") = std::string(lammps_parser::DEFAULT_FORMAT),
          "Parse a single frame and return a column-major dict suitable for NumPy.");
        }

    // install helper (use hasattr to avoid deprecated operator)
//     m.def("install_into_active_env",
//           [m](bool overwrite) -> std::string {
//               if (!py::hasattr(m, "__file__")) throw std::runtime_error("Cannot determine module __file__");
//               std::string src = py::str(m.attr("__file__"));
//               py::module sysconfig = py::module::import("sysconfig");
//               std::string purelib;
//               try {
//                   py::object pathobj = sysconfig.attr("get_path")(std::string("purelib"));
//                   purelib = py::str(pathobj);
//               } catch (const py::error_already_set &e) {
//                   throw std::runtime_error(std::string("Failed to get purelib path: ") + e.what());
//               }

//               fs::path dest_dir(purelib);
//               if (!fs::exists(dest_dir)) {
//                   py::module site = py::module::import("site");
//                   py::object listobj = site.attr("getsitepackages")();
//                   if (py::len(listobj) > 0) dest_dir = fs::path(py::str(listobj[0]));
//               }

//               if (!fs::exists(dest_dir)) throw std::runtime_error("Destination site-packages does not exist: " + dest_dir.string());
//               fs::path src_path(src);
//               fs::path dest_path = dest_dir / src_path.filename();
//               if (fs::exists(dest_path) && !overwrite) throw std::runtime_error("Destination file exists: " + dest_path.string() + ". Use overwrite=True to replace.");
//               try { fs::copy_file(src_path, dest_path, fs::copy_options::overwrite_existing); }
//               catch (const std::exception &ex) { throw std::runtime_error(std::string("Failed to copy module file: ") + ex.what()); }
//               return dest_path.string();
//           },
//           py::arg("overwrite") = false,
//           "Install this compiled module into the active Python environment's site-packages (purelib).");
// }
