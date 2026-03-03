#include "lammps_parser.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <cstring>

namespace lammps_parser {

using std::string;
using std::vector;
using std::array;
using std::unordered_map;
namespace fs = std::filesystem;

const char* DEFAULT_FORMAT = "type id xu yu zu mol c_KE c_PE c_stress[1] c_stress[2] c_stress[3] c_stress[4] c_stress[5] c_stress[6]";

// ----------------- small helpers -----------------
static inline vector<string> split_ws(const string &s) {
    vector<string> out;
    std::istringstream iss(s);
    string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}
static inline string lower_copy(const string &s) {
    string r = s;
    for (char &c : r) c = char(std::tolower((unsigned char)c));
    return r;
}
static inline long long to_ll(const string &s) { try { return std::stoll(s); } catch(...) { return 0; } }
static inline int to_int(const string &s) { try { return std::stoi(s); } catch(...) { return 0; } }
static inline double to_double(const string &s) { try { return std::stod(s); } catch(...) { return 0.0; } }
static inline string normalize_col(const string &col) { return lower_copy(col); }

// parse column with optional index: "c_stress[1]" or "c_stress1"
static inline std::pair<string,int> parse_col_index(const string &col) {
    string s = col;
    auto lb = s.find('[');
    auto rb = s.find(']');
    if (lb != string::npos && rb != string::npos && rb > lb) {
        string base = s.substr(0, lb);
        string idxs = s.substr(lb+1, rb-lb-1);
        int idx = to_int(idxs);
        return {base, idx};
    }
    // trailing digits
    int p = (int)s.size()-1;
    while (p >= 0 && std::isdigit((unsigned char)s[p])) --p;
    if (p < (int)s.size()-1) {
        string base = s.substr(0, p+1);
        string idxs = s.substr(p+1);
        int idx = to_int(idxs);
        if (idx > 0) return {base, idx};
    }
    return {s, 0};
}

// ----------------- header parsing -----------------
// read header until we encounter the `ITEM: ATOMS` line. If we find it, return it in out_item_atoms_line
static Header parseHeaderFromStream(std::istream &in, string &out_item_atoms_line) {
    Header header;
    string line;
    while (std::getline(in, line)) {
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos == string::npos) continue;
        string trimmed = line.substr(pos);
        string low = lower_copy(trimmed);
        if (low.rfind("item:", 0) == 0) {
            vector<string> parts = split_ws(trimmed);
            if (parts.size() >= 2) {
                string itemname = lower_copy(parts[1]);
                if (itemname == "timestep") {
                    if (!std::getline(in, line)) break;
                    auto tparts = split_ws(line);
                    if (!tparts.empty()) header.timestep = to_ll(tparts[0]);
                    continue;
                } else if (itemname == "number") {
                    if (low.find("atoms") != string::npos) {
                        if (!std::getline(in, line)) break;
                        auto nparts = split_ws(line);
                        if (!nparts.empty()) header.n_atoms = (size_t)to_ll(nparts[0]);
                        continue;
                    }
                } else if (low.find("box bounds") != string::npos || (parts.size() >= 3 && lower_copy(parts[1]) == "box" && lower_copy(parts[2]) == "bounds")) {
                    // capture flags after BOX BOUNDS
                    size_t p = lower_copy(trimmed).find("box bounds");
                    if (p != string::npos) {
                        header.box_bounds_flags = trimmed.substr(p + strlen("box bounds"));
                        size_t s2 = header.box_bounds_flags.find_first_not_of(" \t");
                        if (s2 != string::npos) header.box_bounds_flags = header.box_bounds_flags.substr(s2);
                        else header.box_bounds_flags.clear();
                    }
                    for (int i=0;i<3;++i) {
                        if (!std::getline(in, line)) break;
                        auto b = split_ws(line);
                        if (b.size() >= 2) {
                            header.box_bounds[i].lo = to_double(b[0]);
                            header.box_bounds[i].hi = to_double(b[1]);
                            header.box_bounds[i].tilt = (b.size() >= 3 ? to_double(b[2]) : 0.0);
                        } else {
                            header.box_bounds[i].lo = header.box_bounds[i].hi = header.box_bounds[i].tilt = 0.0;
                        }
                    }
                    continue;
                } else if (low.rfind("item: atoms", 0) == 0) {
                    out_item_atoms_line = trimmed;
                    return header;
                }
            }
            if (low.rfind("item: atoms", 0) == 0) {
                out_item_atoms_line = trimmed;
                return header;
            }
        }
    }
    return header;
}

// ----------------- atoms parsing -----------------
static vector<Atom> parseAtomsFromStream_afterHeader(std::istream &in, const string &item_atoms_line, const string &format_override) {
    // determine columns
    vector<string> cols;
    if (!item_atoms_line.empty()) {
        const string prefix = "ITEM: ATOMS";
        string low = lower_copy(item_atoms_line);
        size_t p = low.find(lower_copy(prefix));
        string after;
        if (p != string::npos) after = item_atoms_line.substr(p + prefix.size());
        else {
            auto parts = split_ws(item_atoms_line);
            bool saw_item = false, saw_atoms = false;
            for (auto &t : parts) {
                string tl = lower_copy(t);
                if (!saw_item && tl == "item") { saw_item = true; continue; }
                if (!saw_atoms && tl == "atoms") { saw_atoms = true; continue; }
                if (saw_atoms) { if (!after.empty()) after += " "; after += t; }
            }
        }
        if (!after.empty()) {
            auto toks = split_ws(after);
            for (auto &t : toks) cols.push_back(normalize_col(t));
        }
    }
    if (cols.empty()) {
        auto toks = split_ws(format_override);
        for (auto &t : toks) cols.push_back(normalize_col(t));
    }

    vector<Atom> atoms;
    string line;
    while (std::getline(in, line)) {
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos == string::npos) continue;
        string trimmed = line.substr(pos);
        if (lower_copy(trimmed).rfind("item:", 0) == 0) break; // next block starts
        auto toks = split_ws(trimmed);
        if (toks.empty()) continue;
        Atom a; a.raw_tokens = toks;
        size_t ncols = cols.size();
        for (size_t i=0; i<ncols && i<toks.size(); ++i) {
            string col = cols[i];
            string val = toks[i];
            string colnorm = col;
            if (colnorm == "id") a.id = to_ll(val);
            else if (colnorm == "type") a.type = to_int(val);
            else if (colnorm == "xu" || colnorm == "x") a.xu = to_double(val);
            else if (colnorm == "yu" || colnorm == "y") a.yu = to_double(val);
            else if (colnorm == "zu" || colnorm == "z") a.zu = to_double(val);
            else if (colnorm == "mol" || colnorm == "molecule") a.mol = to_ll(val);
            else if (colnorm == "c_ke" || colnorm == "ke") a.c_KE = to_double(val);
            else if (colnorm == "c_pe" || colnorm == "pe") a.c_PE = to_double(val);
            else {
                auto pi = parse_col_index(colnorm);
                string base = pi.first;
                int idx = pi.second;
                if ((base == "c_stress" || base == "stress" || base == "cstress") && idx >= 1 && idx <= 6)
                    a.c_stress[idx-1] = to_double(val);
                else
                    a.extras[col] = to_double(val);
            }
        }
        if (toks.size() > cols.size()) {
            for (size_t j=cols.size(); j<toks.size(); ++j) {
                string key = "extra_col_" + std::to_string(j - cols.size() + 1);
                a.extras[key] = to_double(toks[j]);
            }
        }
        a.positions = {a.xu, a.yu,a.zu};
        atoms.push_back(std::move(a));
    }

    std::sort(atoms.begin(), atoms.end(), [](const Atom &A, const Atom &B){
        if (A.mol != B.mol) return A.mol < B.mol;
        return A.id < B.id;
    });
    return atoms;
}

// ----------------- public functions -----------------
ParsedFrame parse_frame_from_file(const string &filepath, const string &format_override) {
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) throw std::runtime_error("Unable to open file: " + filepath);
    string item_atoms_line;
    Header header = parseHeaderFromStream(ifs, item_atoms_line);
    if (item_atoms_line.empty()) {
        string line;
        while (std::getline(ifs, line)) {
            size_t pos = line.find_first_not_of(" \t\r\n");
            if (pos == string::npos) continue;
            string trimmed = line.substr(pos);
            if (lower_copy(trimmed).rfind("item: atoms", 0) == 0) { item_atoms_line = trimmed; break; }
        }
    }
    vector<Atom> atoms = parseAtomsFromStream_afterHeader(ifs, item_atoms_line, format_override);
    ParsedFrame pf;
    pf.header = header;
    pf.atoms = std::move(atoms);
    return pf;
}

vector<ParsedFrame> parse_frames_from_directory(const string &dirpath, const string &format_override) {
    vector<std::pair<string, fs::path>> matches;
    std::regex pattern(R"(^dump\..*\.lammpstrj$)", std::regex::ECMAScript | std::regex::icase);
    for (auto &entry : fs::directory_iterator(dirpath)) {
        if (!entry.is_regular_file()) continue;
        string fname = entry.path().filename().string();
        if (std::regex_match(fname, pattern)) matches.emplace_back(fname, entry.path());
    }
    std::sort(matches.begin(), matches.end(), [](auto &a, auto &b){ return a.first < b.first; });
    vector<ParsedFrame> results;
    for (auto &p : matches) {
        try {
            results.push_back(parse_frame_from_file(p.second.string(), format_override));
        } catch (const std::exception &ex) {
            // Issue a warning and continue
            std::cerr << "Warning: failed to parse file " << p.second << ": " << ex.what() << "\n";
        }
    }
    return results;
}

} // namespace lammps_parser
