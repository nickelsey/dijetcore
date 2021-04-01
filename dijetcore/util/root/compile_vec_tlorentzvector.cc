#include "dijetcore/util/root/compile_vec_tlorentzvector.h"

#include "dijetcore/lib/logging.h"
#include "dijetcore/lib/string/string_utils.h"
#include "dijetcore/lib/random.h"

#include <fstream>
#include <string>

#include "boost/filesystem.hpp"

#include "TROOT.h"

namespace dijetcore {

void GenerateVectorTLorentzVectorDictionary() {

  std::string dijetcore_tlorentzvecdict_base =
      dijetcore::MakeString("vectlorentzvecdict", Counter::instance().counter());
  std::string dijetcore_tlorentzvecdict_name =
      dijetcore::MakeString(dijetcore_tlorentzvecdict_base, ".C");

  std::ofstream writer;
  writer.open(dijetcore_tlorentzvecdict_name.c_str());
  if (writer.is_open()) {
    writer << "#include \"TLorentzVector.h\"\n";
    writer << "#include <vector>\n";
    writer << "#ifdef __MAKECINT__\n";
    writer << "#pragma link C++ class vector<TLorentzVector>+;\n";
    writer << "#endif";
    writer.close();
  } else {
    LOG(ERROR) << "could not open " << dijetcore_tlorentzvecdict_name;
  }
  std::string line_to_process =
      dijetcore::MakeString(".L ", dijetcore_tlorentzvecdict_name, "+");
  gROOT->ProcessLine(line_to_process.c_str());

  // now, do cleanup
  // build the file names that need to be deleted
  boost::filesystem::path base_file = dijetcore_tlorentzvecdict_name;
  boost::filesystem::path c_file =
      dijetcore::MakeString(dijetcore_tlorentzvecdict_base, "_C.d");
  boost::filesystem::path so_file =
      dijetcore::MakeString(dijetcore_tlorentzvecdict_base, "_C.so");
  boost::filesystem::path pcm_file = dijetcore::MakeString(
      dijetcore_tlorentzvecdict_base, "_C_ACLiC_dict_rdict.pcm");

  // attempt to remove them all
  boost::filesystem::remove(base_file);
  boost::filesystem::remove(c_file);
  boost::filesystem::remove(so_file);
  boost::filesystem::remove(pcm_file);
}

} // namespace dijetcore