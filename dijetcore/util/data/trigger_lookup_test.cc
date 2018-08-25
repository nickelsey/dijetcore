#include "gtest/gtest.h"

#include <set>
#include <algorithm>
using std::set;

#include "dijetcore/util/data/trigger_lookup.h"

// testing that all trigger sets defined
// in trigger_lookup.h are defined correctly

TEST(TriggerLookup, AuAu) { 
  
  set<unsigned> y14ht = dijetcore::GetTriggerIDs("y14ht");
  EXPECT_EQ(y14ht, (set<unsigned>{450203, 450213, 450202, 450212}));
  
  set<unsigned> y14mb = dijetcore::GetTriggerIDs("y14mb");
  EXPECT_EQ(y14mb, (set<unsigned>{450008, 450018, 450010, 450020, 450011, 450021}));
  
  set<unsigned> y14ht3 = dijetcore::GetTriggerIDs("y14ht3");
  EXPECT_EQ(y14ht3, (set<unsigned>{450203, 450213}));
  
  set<unsigned> y14ht2 = dijetcore::GetTriggerIDs("y14ht2");
  EXPECT_EQ(y14ht2, (set<unsigned>{450202, 450212}));
  
  set<unsigned> y14all = dijetcore::GetTriggerIDs("y14all");
  EXPECT_EQ(y14all, (set<unsigned>{450203, 450213, 450202, 450212, 450008, 450018, 450010, 450020, 450011, 450021}));
  
  set<unsigned> y7ht = dijetcore::GetTriggerIDs("y7ht");
  EXPECT_EQ(y7ht, (set<unsigned>{200620, 200621, 200211, 200212, 200220, 200221, 200222}));

  set<unsigned> y7mb = dijetcore::GetTriggerIDs("y7mb");
  EXPECT_EQ(y7mb, (set<unsigned>{200001, 200003, 200013}));

  set<unsigned> y7all = dijetcore::GetTriggerIDs("y7all");
  EXPECT_EQ(y7all, (set<unsigned>{200620, 200621, 200211, 200212, 200220, 200221, 200222, 200001, 200003, 200013}));

  set<unsigned> y10ht = dijetcore::GetTriggerIDs("y10ht");
  EXPECT_EQ(y10ht, (set<unsigned>{260504, 260514, 260524}));

  set<unsigned> y10all = dijetcore::GetTriggerIDs("y10all");
  EXPECT_EQ(y10all, (set<unsigned>{260504, 260514, 260524}));

  set<unsigned> y11ht = dijetcore::GetTriggerIDs("y11ht");
  EXPECT_EQ(y11ht, (set<unsigned>{350512, 350502, 350513, 350503, 350514, 350504}));

  set<unsigned> y11mb = dijetcore::GetTriggerIDs("y11mb");
  EXPECT_EQ(y11mb, (set<unsigned>{}));

  set<unsigned> y11npe15 = dijetcore::GetTriggerIDs("y11npe15");
  EXPECT_EQ(y11npe15, (set<unsigned>{350512, 350502}));

  set<unsigned> y11npe18 = dijetcore::GetTriggerIDs("y11npe18");
  EXPECT_EQ(y11npe18, (set<unsigned>{350513, 350503}));

  set<unsigned> y11npe25 = dijetcore::GetTriggerIDs("y11npe25");
  EXPECT_EQ(y11npe25, (set<unsigned>{350514, 350504}));

  set<unsigned> y11all = dijetcore::GetTriggerIDs("y11all");
  EXPECT_EQ(y11all, (set<unsigned>{350512, 350502, 350513, 350503, 350514, 350504}));
}

TEST(TriggerLookup, pp) {
  // p+p
  set<unsigned> y6ht = dijetcore::GetTriggerIDs("y6ppht");
  EXPECT_EQ(y6ht, (set<unsigned>{117211, 117212, 127212, 127213, 137213}));

  set<unsigned> y6jp = dijetcore::GetTriggerIDs("y6ppjp");
  EXPECT_EQ(y6jp, (set<unsigned>{117221, 127221, 137221, 137222}));

  set<unsigned> y6all = dijetcore::GetTriggerIDs("y6ppall");
  EXPECT_EQ(y6all, (set<unsigned>{117211, 117212, 127212, 127213, 137213, 117221, 127221, 137221, 137222}));

  set<unsigned> y8ppht = dijetcore::GetTriggerIDs("y8ppht");
  EXPECT_EQ(y8ppht, (set<unsigned>{220500, 220510, 220520}));

  set<unsigned> y8ppmb = dijetcore::GetTriggerIDs("y8ppmb");
  EXPECT_EQ(y8ppmb, (set<unsigned>{220000}));

  set<unsigned> y8ppht0 = dijetcore::GetTriggerIDs("y8ppht0");
  EXPECT_EQ(y8ppht0, (set<unsigned>{220500}));

  set<unsigned> y8ppht1 = dijetcore::GetTriggerIDs("y8ppht1");
  EXPECT_EQ(y8ppht1, (set<unsigned>{220510}));

  set<unsigned> y8ppht2 = dijetcore::GetTriggerIDs("y8ppht2");
  EXPECT_EQ(y8ppht2, (set<unsigned>{220520}));

  set<unsigned> y8ppall = dijetcore::GetTriggerIDs("y8ppall");
  EXPECT_EQ(y8ppall, (set<unsigned>{220500, 220510, 220520, 220000}));

  set<unsigned> y9ppht = dijetcore::GetTriggerIDs("y9ppht");
  EXPECT_EQ(y9ppht, (set<unsigned>{240530, 240540, 240550, 240560, 240570}));

  set<unsigned> y9ppjp = dijetcore::GetTriggerIDs("y9ppjp");
  EXPECT_EQ(y9ppjp, (set<unsigned>{240410, 240411, 240650, 240651, 250652}));

  set<unsigned> y9ppall = dijetcore::GetTriggerIDs("y9ppall");
  EXPECT_EQ(y9ppall, (set<unsigned>{240530, 240540, 240550, 240560, 240570, 240410, 240411, 240650, 240651, 250652}));

  set<unsigned> y12ppht = dijetcore::GetTriggerIDs("y12ppht");
  EXPECT_EQ(y12ppht, (set<unsigned>{370541, 370542, 370351}));

  set<unsigned> y12ppjp = dijetcore::GetTriggerIDs("y12ppjp");
  EXPECT_EQ(y12ppjp, (set<unsigned>{370621, 370601, 370611}));

  set<unsigned> y12ppjp2 = dijetcore::GetTriggerIDs("y12ppjp2");
  EXPECT_EQ(y12ppjp2, (set<unsigned>{370621}));

  set<unsigned> y12pphm = dijetcore::GetTriggerIDs("y12pphm");
  EXPECT_EQ(y12pphm, (set<unsigned>{370341}));

  set<unsigned> y12ppmb = dijetcore::GetTriggerIDs("y12ppmb");
  EXPECT_EQ(y12ppmb, (set<unsigned>{370011}));

  set<unsigned> y12ppall = dijetcore::GetTriggerIDs("y12ppall");
  EXPECT_EQ(y12ppall, (set<unsigned>{370541, 370542, 370351, 370621, 370601, 370611, 370341, 370011}));
  
}

TEST(TriggerLookup, dAu) {
  set<unsigned> y8dauht = dijetcore::GetTriggerIDs("y8dauht");
  EXPECT_EQ(y8dauht, (set<unsigned>{210500, 210501, 210510, 210511, 210520, 210521, 210541}));

  set<unsigned> y8dauht0 = dijetcore::GetTriggerIDs("y8dauht0");
  EXPECT_EQ(y8dauht0, (set<unsigned>{210500, 210501}));

  set<unsigned> y8dauht1 = dijetcore::GetTriggerIDs("y8dauht1");
  EXPECT_EQ(y8dauht1, (set<unsigned>{210510, 210511}));

  set<unsigned> y8dauht2 = dijetcore::GetTriggerIDs("y8dauht2");
  EXPECT_EQ(y8dauht2, (set<unsigned>{210520, 210521}));

  set<unsigned> y8dauht4 = dijetcore::GetTriggerIDs("y8dauht4");
  EXPECT_EQ(y8dauht4, (set<unsigned>{210541}));

  set<unsigned> y8daumb = dijetcore::GetTriggerIDs("y8daumb");
  EXPECT_EQ(y8daumb, (set<unsigned>{210020}));

  set<unsigned> y8dauall = dijetcore::GetTriggerIDs("y8dauall");
  EXPECT_EQ(y8dauall, (set<unsigned>{210500, 210501, 210510, 210511, 210520, 210521, 210541, 210020}));
}
