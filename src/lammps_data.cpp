#include "lammps_parser.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace lammps_parser {

const char* DEFAULT_ATOMSTYLE = "molecular";

static std::string trim_copy(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

static std::string lower_copy(const std::string& s) {
    std::string r = s;
    for (char& c : r) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return r;
}

static std::vector<std::string> split_ws(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

static bool is_int_token(const std::string& s) {
    try {
        size_t pos = 0;
        std::stoll(s, &pos);
        return pos == s.size();
    } catch (...) {
        return false;
    }
}

static bool is_known_next_section(const std::string& t) {
    std::string low = lower_copy(trim_copy(t));
    if (low.empty()) return false;

    static const char* sections[] = {
        "masses", "atoms", "velocities", "bonds", "angles",
        "dihedrals", "impropers", "pair coeffs", "bond coeffs",
        "angle coeffs", "dihedral coeffs", "improper coeffs"
    };

    for (const char* s : sections) {
        if (low.rfind(s, 0) == 0) return true;
    }
    return false;
}

static std::string parse_style_from_atoms_header(const std::string& line) {
    const auto hash = line.find('#');
    if (hash == std::string::npos) return "";
    return lower_copy(trim_copy(line.substr(hash + 1)));
}

DataFile parse_data_from_file(const std::string& filepath, const std::string& atom_style) {
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        throw std::runtime_error("Unable to open LAMMPS data file: " + filepath);
    }

    DataFile df;
    df.header.atom_style = lower_copy(atom_style.empty() ? DEFAULT_ATOMSTYLE : atom_style);

    std::string line;
    bool in_atoms = false;
    bool just_saw_atoms_header = false;

    while (std::getline(ifs, line)) {
        std::string t = trim_copy(line);
        if (t.empty() || t[0] == '#') continue;

        std::vector<std::string> tok = split_ws(t);
        if (tok.empty()) continue;

        if (!in_atoms) {
            // counts
            if (tok.size() >= 2 && tok[1] == "atoms" && is_int_token(tok[0])) {
                df.header.natoms = std::stoll(tok[0]);
                continue;
            }
            if (tok.size() >= 3 && tok[1] == "atom" && tok[2] == "types" && is_int_token(tok[0])) {
                df.header.atom_types = std::stoll(tok[0]);
                continue;
            }

            // box bounds
            if (tok.size() >= 4 && lower_copy(tok[2]) == "xlo" && lower_copy(tok[3]) == "xhi") {
                df.header.xlo = std::stod(tok[0]);
                df.header.xhi = std::stod(tok[1]);
                continue;
            }
            if (tok.size() >= 4 && lower_copy(tok[2]) == "ylo" && lower_copy(tok[3]) == "yhi") {
                df.header.ylo = std::stod(tok[0]);
                df.header.yhi = std::stod(tok[1]);
                continue;
            }
            if (tok.size() >= 4 && lower_copy(tok[2]) == "zlo" && lower_copy(tok[3]) == "zhi") {
                df.header.zlo = std::stod(tok[0]);
                df.header.zhi = std::stod(tok[1]);
                continue;
            }

            // triclinic tilt factors
            if (tok.size() >= 6 && lower_copy(tok[3]) == "xy" && lower_copy(tok[4]) == "xz" && lower_copy(tok[5]) == "yz") {
                df.header.xy = std::stod(tok[0]);
                df.header.xz = std::stod(tok[1]);
                df.header.yz = std::stod(tok[2]);
                df.header.triclinic = true;
                continue;
            }

            // atoms section
            if (lower_copy(tok[0]) == "atoms") {
                std::string style_from_file = parse_style_from_atoms_header(t);
                if (!style_from_file.empty()) {
                    df.header.atom_style = style_from_file;
                }
                in_atoms = true;
                just_saw_atoms_header = true;
                continue;
            }

            continue;
        }

        // inside Atoms section
        if (just_saw_atoms_header) {
            if (t.empty()) continue;   // skip blank line after "Atoms"
            just_saw_atoms_header = false;
        }

        if (is_known_next_section(t)) break;

        DataAtom a;
        a.raw_tokens = tok;

        const std::string style = df.header.atom_style;
        if (style == "molecular") {
            if (tok.size() < 6) continue;
            a.id = std::stoll(tok[0]);
            a.mol = std::stoll(tok[1]);
            a.type = std::stoll(tok[2]);
            a.x = std::stod(tok[3]);
            a.y = std::stod(tok[4]);
            a.z = std::stod(tok[5]);
        } else if (style == "atomic") {
            if (tok.size() < 5) continue;
            a.id = std::stoll(tok[0]);
            a.type = std::stoll(tok[1]);
            a.x = std::stod(tok[2]);
            a.y = std::stod(tok[3]);
            a.z = std::stod(tok[4]);
        } else {
            throw std::runtime_error(
                "parse_data_from_file currently supports atom styles 'molecular' and 'atomic' only. "
                "Got: " + style
            );
        }

        df.atoms.push_back(std::move(a));
    }

    return df;
}

} // namespace lammps_parser