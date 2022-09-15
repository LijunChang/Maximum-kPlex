#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "Timer.h"
#include "LinearHeap.h"

class Graph {
private:
	std::string dir; //input graph directory
	ui n; //number of nodes of the graph
	ept m; //number of edges of the graph
	ui K; //the value of k in k-plex

	ept *pstart; //offset of neighbors of nodes
	ept *pend; //used in search
	ept *pend_buf;
	ui *edges; //adjacent ids of edges

	std::vector<ui> kplex;

	ui *s_degree;
	ept *s_pstart;
	ept *s_pend;
	ui *s_edges;
	ui *s_peel_sequence;
	ui *s_core;
	char *s_vis;
	ListLinearHeap *s_heap;

	ui *s_edgelist_pointer;
	ui *s_tri_cnt;
	ui *s_edge_list;
	ui *s_active_edgelist;
	char *s_deleted;

public:
	Graph(const char *_dir, const int _K) ;
	~Graph() ;

	void read_graph_binary() ;
	void read_graph() ;

	void output_one_kplex() ;
	void verify_kplex() ;

	void kPlex_exact() ;

private:
	void write_subgraph(ui n, const std::vector<std::pair<int,int> > &edge_list) ;
	void heuristic_kplex_max_degree(ui processed_threshold) ;
	void extract_subgraph(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<int,int> > &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer) ;
	void extract_subgraph_and_prune(ui u, ui *ids, ui &ids_n, ui *rid, std::vector<std::pair<int,int> > &vp, ui *Q, ui* degree, char *exists, ept *pend, char *deleted, ui *edgelist_pointer) ;
	void extract_subgraph_full(const ui *ids, ui ids_n, ui *rid, std::vector<std::pair<int,int> > &vp, char *exists, ept *pstart, ept *pend, ui *edges, char *deleted, ui *edgelist_pointer) ;

	ui degen(ui n, ui *peel_sequence, ui *core, ept *pstart, ui *edges, ui *degree, char *vis, ListLinearHeap *heap, bool output) ;
	void shrink_graph(ui &n, ept &m, ui *peel_sequence, ui *core, ui *out_mapping, ui *in_mapping, ui *rid, ept *pstart, ui *edges, bool output) ;
	void oriented_triangle_counting(ui n, ui m, ui *peel_sequence, ept *pstart, ept *pend, ui *edges, ui *tri_cnt, ui *adj) ;
	void reorganize_oriented_graph(ui n, ui *tri_cnt, ui *edge_list, ept *pstart, ept *pend, ept *pend2, ui *edges, ui *edgelist_pointer, ui *buf) ;
	ept peeling(ui critical_vertex, ListLinearHeap *linear_heap, ui *Qv, ui &Qv_n, ui d_threshold, ui *Qe, bool initialize_Qe, ui t_threshold, ui *tri_cnt, ui *active_edgelist, ui &active_edgelist_n, ui *edge_list, ui *edgelist_pointer, char *deleted, ui *degree, ept *pstart, ept *pend, ui *edges, char *exists) ;
	char find(ui u, ui w, ept &b, ept e, char *deleted, ept &idx, ui *edgelist_pointer, ui *edges) ;

	// functions for subgraph processing
	void load_graph_from_edgelist(ui _n, const std::vector<std::pair<int,int> > &edge_list, ui &n, ept &m, ui *degree, ept *pstart, ui *edges) ;
	void subgraph_prune(ui *ids, ui &_n, std::vector<std::pair<int,int> > &edge_list, ui *rid, ui *Qv, ui *Qe, char *exists) ;
};
#endif
