#ifndef LAMMPS_PARSER_H
#define LAMMPS_PARSER_H

#include <array>
#include <string>
#include <vector>
#include <unordered_map>

namespace lammps_parser {

// A single box bound line (lo, hi, optional tilt)
struct BoxBound {
    double lo = 0.0;
    double hi = 0.0;
    double tilt = 0.0;
};

// Header fields parsed from the file prior to the "ITEM: ATOMS" section.
struct Header {
    long long timestep = 0;
    size_t n_atoms = 0;
    std::array<BoxBound,3> box_bounds;
    std::string box_bounds_flags; // flags after BOX BOUNDS line (e.g., "xy xz yz pp pp pp")
};

// Single atom record (flexible)
struct Atom {
    long long id = 0;
    int type = 0;
    double xu = 0.0, yu = 0.0, zu = 0.0;
    long long mol = 0;
    double c_KE = 0.0, c_PE = 0.0;
    std::array<double,6> c_stress = {0,0,0,0,0,0};
    std::array<double,3> positions = {0,0,0};
    std::vector<std::string> raw_tokens; // original tokens from the line
    std::unordered_map<std::string,double> extras; // unknown/extra columns (name -> value)
};

// The combined parse result for a single frame (header + atom list)
struct ParsedFrame {
    Header header;
    std::vector<Atom> atoms; // sorted by mol then id
};

// Default format string (same as original)
extern const char* DEFAULT_FORMAT;

// Top-level parsing functions
// Parse a single file (single frame). If the file contains an ITEM: ATOMS line that specifies columns,
// that will override the format_override. Otherwise, the provided format_override is used.
ParsedFrame parse_frame_from_file(const std::string &filepath,
                                  const std::string &format_override = std::string(DEFAULT_FORMAT));

// Parse multiple files from a directory matching pattern `^dump\..*\.lammpstrj$` (case-insensitive).
// Each matching file is parsed as a single frame (first frame).
std::vector<ParsedFrame> parse_frames_from_directory(const std::string &dirpath,
                                                     const std::string &format_override = std::string(DEFAULT_FORMAT));

} // namespace lammps_parser

#endif // LAMMPS_PARSER_H
