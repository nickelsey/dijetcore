#ifndef DIJETCORE_LIB_JSON_H
#define DIJETCORE_LIB_JSON_H

#include <dijetcore/lib/json/json.hpp>

namespace dijetcore {
    
    // we use a third-party json parser for handling complex
    // sets of command line configs to allow for simpler submit
    // scripts and easier documentation of analysis settings:
    // https://github.com/nlohmann/json
    
    using json = nlohmann::json;

} // namespace dijetcore

#endif // DIJETCORE_LIB_JSON_H
