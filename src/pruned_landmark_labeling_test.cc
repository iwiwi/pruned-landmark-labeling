#include <gtest/gtest.h>
#include "pruned_landmark_labeling.h"

using namespace std;
using testing::Types;

typedef Types<
  PrunedLandmarkLabeling<0>,
  PrunedLandmarkLabeling<1>,
  PrunedLandmarkLabeling<10>,
  PrunedLandmarkLabeling<50>
  > PllTypes;

template<typename T> class PllTest : public testing::Test {};
TYPED_TEST_CASE(PllTest, PllTypes);

void BFS(const vector<pair<int, int> > &es, int s, vector<int> *distance) {
  int V = s + 1;
  for (size_t i = 0; i < es.size(); ++i) {
    V = max(V, max(es[i].first, es[i].second) + 1);
  }
  vector<vector<int> > adj(V);
  for (size_t i = 0; i < es.size(); ++i) {
    int v = es[i].first, w = es[i].second;
    adj[v].push_back(w);
    adj[w].push_back(v);
  }
  distance->clear();
  distance->resize(V, INT_MAX);

  queue<int> que;
  que.push(s);
  distance->at(s) = 0;
  while (!que.empty()) {
    int v = que.front();
    que.pop();
    for (size_t i = 0; i < adj[v].size(); ++i) {
      int w = adj[v][i];
      if (distance->at(w) == INT_MAX) {
        que.push(w);
        distance->at(w) = distance->at(v) + 1;
      }
    }
  }
}

template<typename TypeParam>
void Test(int V, const vector<pair<int, int> > &es) {
  TypeParam pll1;
  ASSERT_TRUE(pll1.ConstructIndex(es));

  ostringstream oss1;
  ASSERT_TRUE(pll1.StoreIndex(oss1));

  TypeParam pll2;
  istringstream iss1(oss1.str());
  ASSERT_TRUE(pll2.LoadIndex(iss1));

  ostringstream oss2;
  ASSERT_TRUE(pll2.StoreIndex(oss2));
  ASSERT_EQ(oss1.str().length(), oss2.str().length());

  for (int v = 0; v < V; ++v) {
    vector<int> ds;
    BFS(es, v, &ds);
    ds.resize(V, INT_MAX);
    for (int w = 0; w < V; ++w) {
      ASSERT_EQ(ds[w], pll1.QueryDistance(v, w)) << v << " -> " << w;
      ASSERT_EQ(ds[w], pll2.QueryDistance(v, w)) << v << " -> " << w;
    }
  }
}

// Two vertices (0 -- 1)
TYPED_TEST(PllTest, Trivial) {
  vector<pair<int, int> > es;
  es.push_back(make_pair(0, 1));
  Test<TypeParam>(2, es);
}

// Path (0 -- 1 -- 2 --...-- |L-2| -- |L-1|)
TYPED_TEST(PllTest, Path) {
  const int L = 30;
  vector<pair<int, int> > es;
  for (int i = 0; i + 1 < L; ++i) {
    es.push_back(make_pair(i, i + 1));
  }
  Test<TypeParam>(L, es);
}

// Circle (0 -- 1 -- 2 --...-- |L-2| -- |L-1| -- 0)
TYPED_TEST(PllTest, Ring) {
  const int L = 30;
  vector<pair<int, int> > es;
  for (int i = 0; i < L; ++i) {
    es.push_back(make_pair(i, (i + 1) % L));
  }
  Test<TypeParam>(L, es);
}

// Almost empty (for testing disconnected pairs)
TYPED_TEST(PllTest, AlmostEmpty) {
  vector<pair<int, int> > es;
  es.push_back(make_pair(0, 3));
  es.push_back(make_pair(3, 6));
  es.push_back(make_pair(6, 0));
  Test<TypeParam>(10, es);
}

// Erdos-Renyi random graph
TYPED_TEST(PllTest, ErdosRenyiGraph) {
  const int N = 50;
  vector<pair<int, int> > es;
  for (int v = 0; v < N; ++v) {
    for (int w = v + 1; w < N; ++w) {
      double x = rand() / double(RAND_MAX);
      if (x < 0.5) {
        es.push_back(make_pair(v, w));
      }
    }
  }
  Test<TypeParam>(N, es);
}

