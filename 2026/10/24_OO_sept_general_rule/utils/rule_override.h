#ifndef OO_RULE_OVERRIDE_H
#define OO_RULE_OVERRIDE_H

// Rule-based policy override for fit2d.C.
// If OO_RULE_DIR is set, fit2d reads <OO_RULE_DIR>/<yTag>_pt<lo>_<hi>.knobs
// (one "key value" pair per line, '#' comments). Every key that is present
// replaces the legacy per-bin setting; absent keys keep the 19_OO behavior.
// The knob files are written by study/rules.py from bin features only
// (prefit ROOT content); nominal B never enters them.

#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TString.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

struct OORuleOverride {
    bool active = false;
    std::string path, rule;
    std::map<std::string, double> values;

    bool has(const char *key) const { return values.count(key); }
    double get(const char *key, double fallback) const {
        auto it = values.find(key);
        return it == values.end() ? fallback : it->second;
    }
    bool flag(const char *key, bool fallback) const { return get(key, fallback ? 1. : 0.) != 0.; }

    void write(TFile &file) const {
        if (!active) return;
        TDirectory::TContext context(&file);
        TNamed("ruleOverrideFile", path.c_str()).Write("", TObject::kOverwrite);
        TNamed("ruleName", rule.c_str()).Write("", TObject::kOverwrite);
        for (auto &kv : values)
            TParameter<double>(Form("rule_%s", kv.first.c_str()), kv.second).Write("", TObject::kOverwrite);
    }
};

inline OORuleOverride oo_load_rule_override(const char *yTag, double ptLow, double ptHigh) {
    OORuleOverride o;
    const char *dir = std::getenv("OO_RULE_DIR");
    if (!dir || !dir[0]) return o;
    o.path = Form("%s/%s_pt%.2f_%.2f.knobs", dir, yTag, ptLow, ptHigh);
    std::ifstream in(o.path);
    if (!in) throw std::runtime_error("OO_RULE_DIR is set but knob file is missing: " + o.path);
    std::string line;
    while (std::getline(in, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        std::istringstream ss(line);
        std::string key, value;
        if (!(ss >> key >> value)) continue;
        if (key == "rule") { o.rule = value; continue; }
        o.values[key] = std::stod(value);
    }
    o.active = true;
    return o;
}

#endif
