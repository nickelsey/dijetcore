#ifndef DIJETCORE_LIB_JSON_H
#define DIJETCORE_LIB_JSON_H

#include "dijetcore/lib/assert.h"
#include "dijetcore/lib/json/json.hpp"
#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/types.h"

#include <fstream>
#include <set>

namespace dijetcore {

// we use a third-party json parser for handling complex
// sets of command line configs to allow for simpler submit
// scripts and easier documentation of analysis settings:
// https://github.com/nlohmann/json

using json = nlohmann::json;

// helper function for validating a configuration file for a
// specific analysis. Checks that cfg has at least one value for
// each key in required
std::set<string> ValidateConfig(dijetcore::json cfg,
                                std::set<string> required) {
  std::set<string> missing_parameters;

  for (auto& key : required) {
    if (cfg.count(key) == 0) missing_parameters.insert(key);
  }

  return missing_parameters;
}

// attempt to load the config file into a json parser
// if loading succeeds, then check that all required parameters
// are present
json LoadConfig(string config_filename, std::set<string> required_params) {
  json config;
  std::ifstream config_file(config_filename);
  if (!config_file.is_open()) {
    LOG(ERROR) << "configuration file: " << config_filename
               << " could not be opened, exiting";
    DIJETCORE_THROW("could not load config file: ", config_filename);
  }
  config_file >> config;

  // validate config file
  auto missing_config_params =
      dijetcore::ValidateConfig(config, required_params);
  if (missing_config_params.size() > 0) {
    DIJETCORE_THROW(
        "error: configuration file incompatible. Missing required parameters: ",
        missing_config_params);
  }
  return config;
}

}  // namespace dijetcore

#endif  // DIJETCORE_LIB_JSON_H
