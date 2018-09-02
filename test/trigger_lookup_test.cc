#include "jet_analysis/util/trigger_lookup.hh"

#include <set>
#include <algorithm>
#include <iterator>
using std::set;

// testing that all trigger sets defined
// in trigger_lookup.hh are defined correctly

int main() {
  
  // Au+Au
  
  set<unsigned> y14ht = GetTriggerIDs("y14ht");
  if (y14ht != set<unsigned>{450203, 450213, 450202, 450212}) return 1;
  
  set<unsigned> y14mb = GetTriggerIDs("y14mb");
  if (y14mb != set<unsigned>{450008, 450018, 450010, 450020, 450011, 450021}) return 1;
  
  set<unsigned> y14ht3 = GetTriggerIDs("y14ht3");
  if (y14ht3 != set<unsigned>{450203, 450213}) return 1;
  
  set<unsigned> y14ht2 = GetTriggerIDs("y14ht2");
  if (y14ht2 != set<unsigned>{450202, 450212}) return 1;
  
  set<unsigned> y14all = GetTriggerIDs("y14all");
  if (y14all != set<unsigned>{450203, 450213, 450202, 450212, 450008, 450018, 450010, 450020, 450011, 450021}) return 1;
  
  set<unsigned> y7ht = GetTriggerIDs("y7ht");
  if (y7ht != set<unsigned>{200620, 200621, 200211, 200212, 200220, 200221, 200222}) return 1;

  set<unsigned> y7mb = GetTriggerIDs("y7mb");
  if (y7mb != set<unsigned>{200001, 200003, 200013}) return 1;

  set<unsigned> y7all = GetTriggerIDs("y7all");
  if (y7all != set<unsigned>{200620, 200621, 200211, 200212, 200220, 200221, 200222, 200001, 200003, 200013}) return 1;

  set<unsigned> y10ht = GetTriggerIDs("y10ht");
  if (y10ht != set<unsigned>{260504, 260514, 260524}) return 1;

  set<unsigned> y10all = GetTriggerIDs("y10all");
  if (y10all != set<unsigned>{260504, 260514, 260524}) return 1;

  set<unsigned> y11ht = GetTriggerIDs("y11ht");
  if (y11ht != set<unsigned>{350512, 350502, 350513, 350503, 350514, 350504}) return 1;

  set<unsigned> y11mb = GetTriggerIDs("y11mb");
  if (y11mb != set<unsigned>{}) return 1;

  set<unsigned> y11npe15 = GetTriggerIDs("y11npe15");
  if (y11npe15 != set<unsigned>{350512, 350502}) return 1;

  set<unsigned> y11npe18 = GetTriggerIDs("y11npe18");
  if (y11npe18 != set<unsigned>{350513, 350503}) return 1;

  set<unsigned> y11npe25 = GetTriggerIDs("y11npe25");
  if (y11npe25 != set<unsigned>{350514, 350504}) return 1;

  set<unsigned> y11all = GetTriggerIDs("y11all");
  if (y11all != set<unsigned>{350512, 350502, 350513, 350503, 350514, 350504}) return 1;

  // p+p
  set<unsigned> y6ht = GetTriggerIDs("y6ppht");
  if (y6ht != set<unsigned>{117211, 117212, 127212, 127213, 137213}) return 1;

  set<unsigned> y6jp = GetTriggerIDs("y6ppjp");
  if (y6jp != set<unsigned>{117221, 127221, 137221, 137222}) return 1;

  set<unsigned> y6all = GetTriggerIDs("y6ppall");
  if (y6all != set<unsigned>{117211, 117212, 127212, 127213, 137213, 117221, 127221, 137221, 137222}) return 1;

  set<unsigned> y8ppht = GetTriggerIDs("y8ppht");
  if (y8ppht != set<unsigned>{220500, 220510, 220520}) return 1;

  set<unsigned> y8ppmb = GetTriggerIDs("y8ppmb");
  if (y8ppmb != set<unsigned>{220000}) return 1;

  set<unsigned> y8ppht0 = GetTriggerIDs("y8ppht0");
  if (y8ppht0 != set<unsigned>{220500}) return 1;

  set<unsigned> y8ppht1 = GetTriggerIDs("y8ppht1");
  if (y8ppht1 != set<unsigned>{220510}) return 1;

  set<unsigned> y8ppht2 = GetTriggerIDs("y8ppht2");
  if (y8ppht2 != set<unsigned>{220520}) return 1;

  set<unsigned> y8ppall = GetTriggerIDs("y8ppall");
  if (y8ppall != set<unsigned>{220500, 220510, 220520, 220000}) return 1;

  set<unsigned> y9ppht = GetTriggerIDs("y9ppht");
  if (y9ppht != set<unsigned>{240530, 240540, 240550, 240560, 240570}) return 1;

  set<unsigned> y9ppjp = GetTriggerIDs("y9ppjp");
  if (y9ppjp != set<unsigned>{240410, 240411, 240650, 240651, 250652}) return 1;
//
  set<unsigned> y9ppall = GetTriggerIDs("y9ppall");
  if (y9ppall != set<unsigned>{240530, 240540, 240550, 240560, 240570, 240410, 240411, 240650, 240651, 250652}) return 1;

  set<unsigned> y12ppht = GetTriggerIDs("y12ppht");
  if (y12ppht != set<unsigned>{370541, 370542, 370351}) return 1;

  set<unsigned> y12ppjp = GetTriggerIDs("y12ppjp");
  if (y12ppjp != set<unsigned>{370621, 370601, 370611}) return 1;

  set<unsigned> y12ppjp2 = GetTriggerIDs("y12ppjp2");
  if (y12ppjp2 != set<unsigned>{370621}) return 1;

  set<unsigned> y12pphm = GetTriggerIDs("y12pphm");
  if (y12pphm != set<unsigned>{370341}) return 1;

  set<unsigned> y12ppmb = GetTriggerIDs("y12ppmb");
  if (y12ppmb != set<unsigned>{370011}) return 1;

  set<unsigned> y12ppall = GetTriggerIDs("y12ppall");
  if (y12ppall != set<unsigned>{370541, 370542, 370351, 370621, 370601, 370611, 370341, 370011}) return 1;

  // d+Au

  set<unsigned> y8dauht = GetTriggerIDs("y8dauht");
  if (y8dauht != set<unsigned>{210500, 210501, 210510, 210511, 210520, 210521, 210541}) return 1;

  set<unsigned> y8dauht0 = GetTriggerIDs("y8dauht0");
  if (y8dauht0 != set<unsigned>{210500, 210501}) return 1;

  set<unsigned> y8dauht1 = GetTriggerIDs("y8dauht1");
  if (y8dauht1 != set<unsigned>{210510, 210511}) return 1;

  set<unsigned> y8dauht2 = GetTriggerIDs("y8dauht2");
  if (y8dauht2 != set<unsigned>{210520, 210521}) return 1;

  set<unsigned> y8dauht4 = GetTriggerIDs("y8dauht4");
  if (y8dauht4 != set<unsigned>{210541}) return 1;

  set<unsigned> y8daumb = GetTriggerIDs("y8daumb");
  if (y8daumb != set<unsigned>{210020}) return 1;

  set<unsigned> y8dauall = GetTriggerIDs("y8dauall");
  if (y8dauall != set<unsigned>{210500, 210501, 210510, 210511, 210520, 210521, 210541, 210020}) return 1;

  return 0;
}
