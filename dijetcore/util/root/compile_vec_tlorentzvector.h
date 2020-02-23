#ifndef DIJETCORE_UTIL_ROOT_COMPILE_VEC_TLORENTZVECTOR_H
#define DIJETCORE_UTIL_ROOT_COMPILE_VEC_TLORENTZVECTOR_H

namespace dijetcore {
  // use to generate a dictionary during binary execution 
  // to allow TTree's to store vector<TLorentzVector>
  // automatically cleans up the files that are generated

  void GenerateVectorTLorentzVectorDictionary();

}

#endif // DIJETCORE_UTIL_ROOT_COMPILE_VEC_TLORENTZVECTOR_H