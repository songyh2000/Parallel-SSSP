#pragma once
#include <queue>

#include "graph.h"
#include "parlay/internal/get_time.h"

sequence<EdgeTy> dijkstra(size_t s, const Graph &G) {
  sequence<EdgeTy> dist(G.n, DIST_MAX);
  dist[s] = 0;
  priority_queue<pair<EdgeTy, NodeId>, vector<pair<EdgeTy, NodeId>>,
                 greater<pair<EdgeTy, NodeId>>>
      pq;
  pq.push(make_pair(dist[s], s));
  while (!pq.empty()) {
    pair<EdgeTy, NodeId> dist_and_node = pq.top();
    pq.pop();
    EdgeTy d = dist_and_node.first;
    NodeId u = dist_and_node.second;
    if (dist[u] < d) continue;
    for (size_t j = G.offset[u]; j < G.offset[u + 1]; j++) {
      NodeId v = G.edge[j].v;
      EdgeTy w = G.edge[j].w;
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w;
        pq.push(make_pair(dist[v], v));
      }
    }
  }
  return dist;
}

#define BELLMAN_FORD_LOCAL_QUEUE_SIZE 4096

sequence<EdgeTy> bellman_ford_with_local_opt(size_t s_idx, const Graph &G) {
  sequence<EdgeTy> dist(G.n, DIST_MAX);
  
  if (G.n == 0) {
    return dist;
  }

  if (s_idx >= G.n) {
    fprintf(stderr, "Error: Source node %zu is out of bounds (Graph has %zu nodes).\n", s_idx, G.n);
    return dist;
  }

  dist[s_idx] = 0;

  std::queue<NodeId> main_q; // Main SPFA queue
  main_q.push(static_cast<NodeId>(s_idx));

  std::vector<bool> is_in_main_q(G.n, false); 
  is_in_main_q[s_idx] = true;

  // This flag enables the local queue optimization strategy.
  bool super_sparse_mode_enabled = true; 

  // Fixed-size array for the local queue, as in the snippet.
  NodeId local_q_arr[BELLMAN_FORD_LOCAL_QUEUE_SIZE];
  // Outer loop: Iterates for a maximum of G.n "passes" or "rounds".
  for (size_t pass_count = 0; pass_count < G.n; ++pass_count) {
    if (main_q.empty()) {
      break;
    }

    // Inner loop as per your request:
    // WARNING: This loop structure, when combined with main_q.pop() inside,
    // will only process approximately half the elements initially in main_q for this pass.
    for (size_t i = 0; i < main_q.size(); ++i) {
      NodeId f_node_id = main_q.front();
      main_q.pop(); // This reduces main_q.size()
      
      size_t f_idx = static_cast<size_t>(f_node_id);
      is_in_main_q[f_idx] = false;

      size_t degree_f = G.offset[f_idx + 1] - G.offset[f_idx];
      
      if (super_sparse_mode_enabled && degree_f < BELLMAN_FORD_LOCAL_QUEUE_SIZE) {
        size_t lq_front = 0, lq_rear = 0;
        local_q_arr[lq_rear++] = f_node_id; 

        while (lq_front < lq_rear) { 
          NodeId u_local_node_id = local_q_arr[lq_front++];
          size_t u_local_idx = static_cast<size_t>(u_local_node_id);
          size_t deg_u_local = G.offset[u_local_idx + 1] - G.offset[u_local_idx];

          if (deg_u_local >= BELLMAN_FORD_LOCAL_QUEUE_SIZE) {
            if (!is_in_main_q[u_local_idx]) {
              main_q.push(u_local_node_id);
              is_in_main_q[u_local_idx] = true;
            }
            continue; 
          }

          for (EdgeId es = G.offset[u_local_idx]; es < G.offset[u_local_idx + 1]; ++es) {

            NodeId v_neighbor_node_id = G.edge[es].v;
            EdgeTy edge_w = G.edge[es].w;
            size_t v_neighbor_idx = static_cast<size_t>(v_neighbor_node_id);

            if (dist[u_local_idx] + edge_w < dist[v_neighbor_idx]) {
              dist[v_neighbor_idx] = dist[u_local_idx] + edge_w;
              
              if (lq_rear < BELLMAN_FORD_LOCAL_QUEUE_SIZE) {
                local_q_arr[lq_rear++] = v_neighbor_node_id;
              } else { 
                if (!is_in_main_q[v_neighbor_idx]) {
                  main_q.push(v_neighbor_node_id);
                  is_in_main_q[v_neighbor_idx] = true;
                }
              }
            }
          } 
        }

        for (size_t j = lq_front; j < lq_rear; ++j) {
          NodeId node_to_main_q = local_q_arr[j];
          size_t node_to_main_q_idx = static_cast<size_t>(node_to_main_q);
          if (!is_in_main_q[node_to_main_q_idx]) {
            main_q.push(node_to_main_q);
            is_in_main_q[node_to_main_q_idx] = true;
          }
        }

      } else { 
        for (EdgeId es = G.offset[f_idx]; es < G.offset[f_idx + 1]; ++es) {

          NodeId v_neighbor_node_id = G.edge[es].v;
          EdgeTy edge_w = G.edge[es].w;
          size_t v_neighbor_idx = static_cast<size_t>(v_neighbor_node_id);

          if (dist[f_idx] + edge_w < dist[v_neighbor_idx]) {
            dist[v_neighbor_idx] = dist[f_idx] + edge_w;
            if (!is_in_main_q[v_neighbor_idx]) {
              main_q.push(v_neighbor_node_id);
              is_in_main_q[v_neighbor_idx] = true;
            }
          }
        }
      }
    } // End of inner loop (processing nodes from main_q for this pass)
  } // End of outer for loop (passes)
  return dist;
}

void verifier(size_t s, const Graph &G, const sequence<EdgeTy> &act_dist) {
  internal::timer tm;
  auto exp_dist = dijkstra(s, G);
  tm.stop();
  printf("dijkstra running time: %-10f\n", tm.total_time());
  parallel_for(0, G.n, [&](size_t i) {
    if (exp_dist[i] != act_dist[i]) {
      string out =
          "exp_dist[" + to_string(i) + "]: " + to_string(exp_dist[i]) + " vs. ";
      out += "act_dist[" + to_string(i) + "]: " + to_string(act_dist[i]) + "\n";
      cout << out;
    }
    assert(exp_dist[i] == act_dist[i]);
  });
}
