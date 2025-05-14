#pragma once
#include "graph.h"
#include "hashbag.h"
#include "parlay/internal/get_time.h"
using namespace std;
// using namespace parlay;

constexpr int NUM_SRC = 1;
constexpr int NUM_ROUND = 1;
int true_work = 0;

constexpr size_t LOCAL_QUEUE_SIZE = 4096;
constexpr size_t DEG_THLD = 20;
constexpr size_t SSSP_SAMPLES = 1000;

enum Algorithm { rho_stepping = 0, delta_stepping, bellman_ford };

class SSSP {
 protected:
  const Graph &G;
  bool sparse;
  int sd_scale;
  size_t frontier_size;
  hashbag<NodeId> bag;
  sequence<EdgeTy> dist;
  sequence<NodeId> frontier;
  sequence<atomic<bool>> in_frontier;
  sequence<atomic<bool>> in_next_frontier;
  // std::vector<size_t> dist_update_count(G.num_nodes, 0); // [NEW]


  void add_to_frontier(NodeId v) {
    if (sparse) {
      if (!in_frontier[v] &&
          compare_and_swap(&in_next_frontier[v], false, true)) {
        bag.insert(v);
      }
    } else {  // dense
      if (!in_frontier[v] && !in_next_frontier[v]) {
        in_next_frontier[v] = true;
      }
    }
  }

  size_t estimate_size() {
    static uint32_t seed = 10086;
    size_t hits = 0;
    for (size_t i = 0; i < SSSP_SAMPLES; i++) {
      NodeId u = hash32(seed) % G.n;
      if (in_frontier[u]) {
        hits++;
      }
      seed++;
    }
    return hits * G.n / SSSP_SAMPLES;
  }


  size_t sparse_relax() {
    // static uint32_t seed = 353442899;
    // size_t sum_deg = 0;
    // for (size_t i = 0; i < SSSP_SAMPLES; i++) {
    // NodeId u = frontier[hash32(seed) % size];
    // sum_deg += G.offset[u + 1] - G.offset[u];
    // seed++;
    //}
    // size_t avg_deg = sum_deg / SSSP_SAMPLES;
    // bool super_sparse = (avg_deg <= DEG_THLD);
    bool super_sparse = true;

    EdgeTy th = get_threshold();

    // for(int i = 0; i < frontier_size; i++) {
    for(size_t i = 0; i < frontier_size; i++) {
      // std::cout << "frontier[i] : " << frontier[i] << std::endl;
      NodeId f = frontier[i];
      in_frontier[f] = false;
      if (dist[f] > th) {
        std::cout << "f : " << f;
        add_to_frontier(f);
        true_work++;
      } else {
        size_t _n = G.offset[f + 1] - G.offset[f];
        if (super_sparse && _n < LOCAL_QUEUE_SIZE) {
          NodeId local_queue[LOCAL_QUEUE_SIZE];
          size_t front = 0, rear = 0;
          local_queue[rear++] = f;
          while (front < rear && rear < LOCAL_QUEUE_SIZE) {
            NodeId u = local_queue[front++];
            size_t deg = G.offset[u + 1] - G.offset[u];
            if (deg >= LOCAL_QUEUE_SIZE || dist[u] > th) {
              std::cout << "u : " << u;
              std::cout << "dist[u] : " << dist[u] << std::endl;
              add_to_frontier(u);
              true_work++;
              continue;
            }
            if (G.symmetrized) {
              EdgeTy temp_dis = dist[u];
              for (EdgeId es = G.offset[u]; es < G.offset[u + 1]; es++) {
                NodeId v = G.edge[es].v;
                EdgeTy w = G.edge[es].w;
                temp_dis = min(temp_dis, dist[v] + w);
              }
              write_min(&dist[u], temp_dis,
                        [](EdgeTy w1, EdgeTy w2) { return w1 < w2; });
            }
            for (EdgeId es = G.offset[u]; es < G.offset[u + 1]; es++) {
              true_work++;
              NodeId v = G.edge[es].v;
              EdgeTy w = G.edge[es].w;
              // std::cout << "u: " << u << " v: " << v << std::endl;
              // std::cout << "dist[u]: " << dist[u] << " dist[v]: " << dist[v] << std::endl;
              // std::cout << "w: " << w << std::endl;
              if (write_min(&dist[v], dist[u] + w,
                            [](EdgeTy w1, EdgeTy w2) { return w1 < w2; })) {
                if (rear < LOCAL_QUEUE_SIZE) {
                  local_queue[rear++] = v;
                } else {
                  add_to_frontier(v);
                  std::cout << "v : " << v;
                }
              }
            }
          }
          for (size_t j = front; j < rear; j++) {
            true_work++;
            add_to_frontier(local_queue[j]);
            std::cout << "local_queue[j] : " << local_queue[j] << std::endl;
          }
        } else {
          // if (G.symmetrized) {
          // size_t deg = G.offset[f + 1] - G.offset[f];
          // auto _dist = delayed_seq<EdgeTy>(deg, [&](size_t es) {
          // NodeId v = G.edge[G.offset[f] + es].v;
          // EdgeTy w = G.edge[G.offset[f] + es].w;
          // return dist[v] + w;
          //});
          // EdgeTy temp_dist = *min_element(_dist);
          // write_min(&dist[f], temp_dist,
          //[](EdgeTy w1, EdgeTy w2) { return w1 < w2; });
          //}
          size_t start_index = G.offset[f];
          size_t end_index = G.offset[f + 1];
          std::cout << "start_index : " << start_index;

          for (size_t es = start_index; es < end_index; ++es) {
              true_work++;
              if (G.symmetrized) {
                  EdgeTy temp_dist = dist[f];
                  for (size_t inner_es = start_index; inner_es < end_index; ++inner_es) { // Corrected loop variable
                      NodeId v = G.edge[inner_es].v; // Use inner_es here
                      EdgeTy w = G.edge[inner_es].w;  // Use inner_es here
                      temp_dist = std::min(temp_dist, dist[v] + w);
                  }
                  if (temp_dist < dist[f]) { // Simulate write_min
                                      std::cout << "v : " << f;
                      std::cout << "dist[v] : " << dist[f] << std::endl;
                      dist[f] = temp_dist;
                      add_to_frontier(f);
                  }
              }

              NodeId v = G.edge[es].v;
              EdgeTy w = G.edge[es].w;
              if (dist[f] + w < dist[v]) { // Simulate write_min
                  dist[v] = dist[f] + w;
                  std::cout << "v : " << v;
                  std::cout << "dist[v] : " << dist[v] << std::endl;
                  add_to_frontier(v);
              }
          }
        }
      }
    };
    swap(in_frontier, in_next_frontier);
    return bag.pack_into(make_slice(frontier));
  }

  size_t dense_relax() {
    while (estimate_size() >= G.n / sd_scale) {
      EdgeTy th = get_threshold();
      // parallel_for(0, G.n, [&](NodeId u) {
      for (NodeId u = 0; u < G.n; u++) {
        if (in_frontier[u]) {
          in_frontier[u] = false;
          if (dist[u] > th) {
            in_next_frontier[u] = true;
            true_work++;
          } else {
            // if (G.symmetrized) {
            // size_t deg = G.offset[u + 1] - G.offset[u];
            // auto _dist = delayed_seq<EdgeTy>(deg, [&](size_t es) {
            // NodeId v = G.edge[G.offset[u] + es].v;
            // EdgeTy w = G.edge[G.offset[u] + es].w;
            // return dist[v] + w;
            //});
            // EdgeTy temp_dist = *min_element(_dist);
            // write_min(&dist[u], temp_dist,
            //[](EdgeTy w1, EdgeTy w2) { return w1 < w2; });
            //}
            size_t start_index = G.offset[u];
            size_t end_index = G.offset[u + 1];

            for (size_t es = start_index; es < end_index; ++es) {
                true_work++;
                if (G.symmetrized) {
                    EdgeTy temp_dist = dist[u];
                    for (size_t inner_es = start_index; inner_es < end_index; ++inner_es) {
                        NodeId v = G.edge[inner_es].v;
                        EdgeTy w = G.edge[inner_es].w;
                        temp_dist = std::min(temp_dist, dist[v] + w);
                    }
                    if (temp_dist < dist[u]) {
                        dist[u] = temp_dist;
                        std::cout << "v : " << u;
                        std::cout << "dist[v] : " << dist[u] << std::endl;
                        add_to_frontier(u);
                    }
                }
                NodeId v = G.edge[es].v;
                EdgeTy w = G.edge[es].w;
                if (dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    std::cout << "v : " << v;
                    std::cout << "dist[v] : " << dist[v] << std::endl;
                    add_to_frontier(v);
                }
            }

          }
        }
      };
      swap(in_frontier, in_next_frontier);
    }
    return count(in_frontier, true);
  }

  void sparse2dense() {
    // parallel_for(0, frontier_size, [&](size_t i) {
    // NodeId u = frontier[i];
    // assert(in_frontier[u] == true);
    //  in_frontier[u] = true;
    //});
  }

  void dense2sparse() {
    auto identity = delayed_seq<NodeId>(G.n, [&](NodeId i) { return i; });
    pack_into_uninitialized(identity, in_frontier, frontier);
  }

  virtual void init() = 0;
  virtual EdgeTy get_threshold() = 0;

 public:
  SSSP() = delete;
  SSSP(const Graph &_G) : G(_G), bag(G.n) {
    dist = sequence<EdgeTy>::uninitialized(G.n);
    frontier = sequence<NodeId>::uninitialized(G.n);
    in_frontier = sequence<atomic<bool>>::uninitialized(G.n);
    in_next_frontier = sequence<atomic<bool>>::uninitialized(G.n);
  }
  sequence<EdgeTy> sssp(NodeId s) {
    if (!G.weighted) {
      fprintf(stderr, "Error: Input graph is unweighted\n");
      exit(EXIT_FAILURE);
    }

    init();
    // parallel_for(0, G.n, [&](NodeId i) {
    for (NodeId i = 0; i < G.n; i++) {
      dist[i] = numeric_limits<EdgeTy>::max() / 2;
      in_frontier[i] = in_next_frontier[i] = false;
    };
    frontier_size = 1;
    dist[s] = 0;
    frontier[0] = s;
    in_frontier[s] = true;
    sparse = true;

    int round = 0;
    int work = 0;
    while (frontier_size) {
      // printf("Round %d: %s, size: %zu, ", round++, sparse ? "sparse" :
      // "dense", frontier_size); internal::timer t;
      if (round % 10 == 0) {
        // printf("Round %d: %s, size: %zu, ", round++, sparse ? "sparse" :
        // "dense", frontier_size);
        std::cout << "Round " << round << ": " << (sparse ? "sparse" : "dense") << ", size: " << frontier_size << std::endl;
      }
      round++;
      if (sparse) {
        printf("sparse\n");
        frontier_size = sparse_relax();
      } else {
        printf("dense\n");
        frontier_size = dense_relax();
      }
      // frontier_size = dense_relax();
      work += frontier_size;
      // printf("relax: %f, ", t.next_time());
      bool next_sparse = (frontier_size < G.n / sd_scale) ? true : false;
      if (sparse && !next_sparse) {
        sparse2dense();
      } else if (!sparse && next_sparse) {
        dense2sparse();
      }
      // printf("pack: %f\n", t.next_time());
      sparse = next_sparse;
    }
    std::cout << "Total work: " << work << ", Total round: " << round << std::endl;
    std::cout << "True work: " << true_work << std::endl;
    true_work = 0;
    return dist;
  }

  void set_sd_scale(int x) { sd_scale = x; }
};

class Rho_Stepping : public SSSP {
  size_t rho;
  uint32_t seed;

 public:
  Rho_Stepping(const Graph &_G, size_t _rho = 1 << 20) : SSSP(_G), rho(_rho) {
    seed = 0;
  }
  void init() override {}
  EdgeTy get_threshold() override {
    if (frontier_size <= rho) {
      if (sparse) {
        auto _dist = delayed_seq<EdgeTy>(
            frontier_size, [&](size_t i) { return dist[frontier[i]]; });
        return *max_element(_dist);
      } else {
        return DIST_MAX;
      }
    }
    EdgeTy sample_dist[SSSP_SAMPLES + 1];
    for (size_t i = 0; i <= SSSP_SAMPLES; i++) {
      if (sparse) {
        NodeId v = frontier[hash32(seed + i) % frontier_size];
        sample_dist[i] = dist[v];
      } else {
        NodeId v = hash32(seed + i) % G.n;
        if (in_frontier[v]) {
          sample_dist[i] = dist[v];
        } else {
          sample_dist[i] = DIST_MAX;
        }
      }
    }
    seed += SSSP_SAMPLES + 1;
    size_t id = 1.0 * rho / frontier_size * SSSP_SAMPLES;
    sort(sample_dist, sample_dist + SSSP_SAMPLES + 1);
    return sample_dist[id];
  }
};

class Delta_Stepping : public SSSP {
  EdgeTy delta;
  EdgeTy thres;

 public:
  Delta_Stepping(const Graph &_G, EdgeTy _delta = 1 << 15)
      : SSSP(_G), delta(_delta) {}
  void init() override { thres = 0; }
  EdgeTy get_threshold() override {
    thres += delta;
    return thres;
  }
};

class Bellman_Ford : public SSSP {
 public:
  Bellman_Ford(const Graph &_G) : SSSP(_G) {}
  void init() override {}
  EdgeTy get_threshold() override { return DIST_MAX; }
};