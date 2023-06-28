#ifndef _KPLEX_BB_MATRIX_
#define _KPLEX_BB_MATRIX_

#include "Utility.h"
#include "Timer.h"

#define _SECOND_ORDER_PRUNING_

//#define set_bit(array, pos) (((array)[pos]) = 1)
//#define reverse_bit(array, pos) (((array)[pos]) = 1- ((array)[pos]))
//#define test_bit(array, pos) ((array)[pos])

class KPLEX_BB_MATRIX {
private:
    long long n;

    char *matrix;
    long long matrix_size;

#ifdef _SECOND_ORDER_PRUNING_
    ui *cn;
    std::queue<std::pair<ui,ui> > Qe;
    std::vector<std::pair<ui, ui> > removed_edges;
    long long removed_edges_n;
#endif

    ui *degree;
    ui *degree_in_S;

    ui K;
    ui *best_solution;
    ui best_solution_size;

    ui *neighbors;
    ui *nonneighbors;

    ui *SR; // union of S and R, where S is at the front
    ui *SR_rid; // reverse ID for SR
    std::queue<ui> Qv;
    ui *level_id;

    std::vector<ui> must_include_vertices;

public:
    KPLEX_BB_MATRIX() {
    	n = 0;
        matrix = NULL;
        matrix_size = 0;

#ifdef _SECOND_ORDER_PRUNING_
        cn = NULL;
        removed_edges_n = 0;
#endif

        degree = degree_in_S = NULL;
        
        best_solution = NULL;
        K = best_solution_size = 0;

        neighbors = nonneighbors = NULL;

        SR = SR_rid = NULL;
        level_id = NULL;
    }

    ~KPLEX_BB_MATRIX() {
        if(matrix != NULL) {
            delete[] matrix;
            matrix = NULL;
        }
#ifdef _SECOND_ORDER_PRUNING_
        if(cn != NULL) {
        	delete[] cn;
        	cn = NULL;
        }
#endif
        if(degree != NULL) {
            delete[] degree;
            degree = NULL;
        }
        if(degree_in_S != NULL) {
            delete[] degree_in_S;
            degree_in_S = NULL;
        }
        if(best_solution != NULL) {
            delete[] best_solution;
            best_solution = NULL;
        }
        if(SR != NULL) {
            delete[] SR;
            SR = NULL;
        }
        if(SR_rid != NULL) {
            delete[] SR_rid;
            SR_rid = NULL;
        }
        if(neighbors != NULL) {
            delete[] neighbors;
            neighbors = NULL;
        }
        if(nonneighbors != NULL) {
            delete[] nonneighbors;
            nonneighbors = NULL;
        }
        if(level_id != NULL) {
        	delete[] level_id;
        	level_id = NULL;
        }
    }

    void allocateMemory(ui n, ui m) {
	    if(n <= 0||m <= 0) return ;

        matrix = new char[m*2];
#ifdef _SECOND_ORDER_PRUNING_
        cn = new ui[m*2];
#endif
        matrix_size = m*2;
        //printf("matrix size: %lldMB\n", (((long long)UB)*v_n*4)/1024/1024);
        
        degree = new ui[n];
        degree_in_S = new ui[n];
        best_solution = new ui[n];
        SR = new ui[n];
        SR_rid = new ui[n];
        neighbors = new ui[n];
        nonneighbors = new ui[n];
        level_id = new ui[n];
    }

    void load_graph(ui _n, const std::vector<std::pair<int,int> > &vp) {
        n = _n;
        if(((long long)n)*n > matrix_size) {
        	do {
        		matrix_size *= 2;
        	} while(((long long)n)*n > matrix_size);
        	delete[] matrix; matrix = new char[matrix_size];
#ifdef _SECOND_ORDER_PRUNING_
        	delete[] cn; cn = new ui[matrix_size];
#endif
        }

#ifdef _SECOND_ORDER_PRUNING_
        memset(cn, 0, sizeof(ui)*((long long)n)*n);
#endif
        memset(matrix, 0, sizeof(char)*((long long)n)*n);
        for(ui i = 0; i < n; i++) degree[i] = 0;
        for(ui i = 0;i < vp.size();i ++) {
            assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
        	ui a = vp[i].first, b = vp[i].second;
            degree[a] ++;
            degree[b] ++;
            matrix[a*n + b] = matrix[b*n + a] = 1;
        }

#ifndef NDEBUG
        printf("load graph of size n=%lld, m=%lu\n", n, vp.size());
        //for(ui i = 0;i < vp.size();i ++) printf("%d %d\n", vp[i].first, vp[i].second);
#endif
    }

    void kPlex(ui K_, std::vector<ui> &kplex, bool must_include_0) {
        K = K_;
        if(K == 1) {
        	printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
        	return ;
        }
        best_solution_size = kplex.size();
        ui R_end;
        initialization(R_end, must_include_0);
        if(R_end) BB_search(0, R_end, 1, must_include_0, 0);
        if(best_solution_size > kplex.size()) {
            kplex.clear();
            for(int i = 0;i < best_solution_size;i ++) kplex.push_back(best_solution[i]);
        }
    }

    int main(int argc, char *argv[]) {
	    if(argc < 3) {
		    printf("Usage: [1]exe [2]dir [3]k\n");
		    return 0;
	    }
        readGraph_binary(argv[1]);
        printf("Finish reading graph\n");
        K = atoi(argv[2]);
        if(K == 1) {
        	printf("For the special case of computing maximum clique, please invoke SOTA maximum clique solver!\n");
        	return 0;
        }
        best_solution_size = 1;
        ui R_end;
        Timer t;
        initialization(R_end, false);
        if(R_end) BB_search(0, R_end, 1, 0, 0);
        printf("Maximum %u-plex size: %u, time excluding reading: %s (micro seconds)\n", K, best_solution_size, Utility::integer_to_string(t.elapsed()).c_str());
        return 0;
    }

private:
    void readGraph_binary(char* dir) {
	    FILE *f = Utility::open_file( (std::string(dir) + std::string("/b_degree.bin")).c_str(), "rb");

	    int tt;
	    fread(&tt, sizeof(int), 1, f);
	    if(tt != sizeof(int)) {
		    printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		    return ;
	    }
	    ui m;
	    fread(&n, sizeof(int), 1, f);
	    fread(&m, sizeof(int), 1, f);
	    printf("n = %lld, m = %u\n", n, m/2);

        allocateMemory(n, (n*n+1)/2);

	    fread(degree, sizeof(ui), n, f);
	    fclose(f);

	    f = Utility::open_file( (std::string(dir) + std::string("/b_adj.bin")).c_str(), "rb");

#ifdef _SECOND_ORDER_PRUNING_
        memset(cn, 0, sizeof(ui)*n*n);
#endif
        memset(matrix, 0, sizeof(char)*n*n);
	    for(ui i = 0;i < n;i ++) if(degree[i] > 0) {
		    fread(neighbors, sizeof(ui), degree[i], f);
		    for(ui j = 0;j < degree[i];j ++) matrix[i*n + neighbors[j]] = 1;
		    ui d = 0;
		    for(ui j = 0;j < n;j ++) if(matrix[i*n + j]) ++ d;
		    if(d != degree[i]) {
		    	printf("%u may have duplicate edges\n", i);
		    	degree[i] = d;
		    }
		}
	    fclose(f);
    }

    void initialization(ui &R_end, bool must_include_0) {
    	// the following computes a degeneracy ordering and a heuristic solution
    	ui *peel_sequence = neighbors;
    	ui *core = nonneighbors;
    	ui *vis = SR;
    	memset(vis, 0, sizeof(ui)*n);
    	ui max_core = 0, UB = 0, idx = n;
    	for(ui i = 0;i < n;i ++) {
			ui u, min_degree = n;
			for(ui j = 0;j < n;j ++) if(!vis[j]&&degree[j] < min_degree) {
				u = j;
				min_degree = degree[j];
			}
			if(min_degree > max_core) max_core = min_degree;
			core[u] = max_core;
			peel_sequence[i] = u;
			vis[u] = 1;

			ui t_UB = core[u] + K;
			if(n - i < t_UB) t_UB = n - i;
			if(t_UB > UB) UB = t_UB;

			if(idx == n&&min_degree + K >= n - i) idx = i;

			for(ui j = 0;j < n;j ++) if(!vis[j]&&matrix[u*n + j]) -- degree[j];
		}
    	if(n - idx > best_solution_size) {
    		best_solution_size = n - idx;
    		for(ui i = idx;i < n;i ++) best_solution[i-idx] = peel_sequence[i];
    		printf("Degen find a solution of size %u\n", best_solution_size);
    	}

        memset(degree_in_S, 0, sizeof(ui)*n);
        R_end = 0;
        for(ui i = 0;i < n;i ++) SR_rid[i] = n;
        for(ui i = 0;i < n;i ++) if(core[i] + K > best_solution_size) {
            SR[R_end] = i; SR_rid[i] = R_end;
            ++ R_end;
        }

        if(must_include_0&&SR_rid[0] == n) {
        	R_end = 0;
        	return ;
        }

        for(ui i = 0;i < R_end;i ++) {
        	ui u = SR[i];
        	degree[u] = 0;
        	for(ui j = 0;j < R_end;j ++) if(matrix[u*n + SR[j]]) ++ degree[u];
        }

        memset(level_id, 0, sizeof(ui)*n);
        for(ui i = 0;i < R_end;i ++) level_id[SR[i]] = n;

        assert(Qv.empty());

#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = 0;i < R_end;i ++) {
        	ui neighbors_n = 0;
        	char *t_matrix = matrix + SR[i]*n;
        	for(ui j = 0;j < R_end;j ++) if(t_matrix[SR[j]]) neighbors[neighbors_n ++] = SR[j];
        	for(ui j = 0;j < neighbors_n;j ++) for(ui k = j+1;k < neighbors_n;k ++) {
        		++ cn[neighbors[j]*n + neighbors[k]];
        		++ cn[neighbors[k]*n + neighbors[j]];
        	}
        }

        while(!Qe.empty()) Qe.pop();
        for(ui i = 0;i < R_end;i ++) for(ui j = i+1;j < R_end;j ++) {
        	if(matrix[SR[i]*n + SR[j]]&&upper_bound_based_prune(0, SR[i], SR[j])) {
        		Qe.push(std::make_pair(SR[i], SR[j]));
        	}
        }
        removed_edges_n = 0;
#endif

        must_include_vertices.resize(n);

        if(remove_vertices_and_edges_with_prune(0, R_end, 0)) R_end = 0;
    }

    void BB_search(ui S_end, ui R_end, ui level, bool choose_zero, ui must_include_vertices_n) {
    	assert(!choose_zero||!must_include_vertices_n);
        if(S_end > best_solution_size) { // find a larger solution
            best_solution_size = S_end;
            for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
#ifndef NDEBUG
            for(ui i = 0;i < best_solution_size;i ++) assert(degree_in_S[best_solution[i]] + K >= best_solution_size);
            printf("Find a k-plex of size: %u\n", best_solution_size);
#endif
        }
        if(R_end <= best_solution_size) return ;

#ifndef NDEBUG
        for(int i = 0;i < S_end;i ++) assert(degree[SR[i]] + K > best_solution_size&&degree_in_S[SR[i]] + K >= S_end);
#endif

        bool is_kplex = true;
        for(ui i = 0;i < R_end;i ++) if(degree[SR[i]] + K < R_end) {
            is_kplex = false;
            break;
        }
        if(is_kplex) {
            best_solution_size = R_end;
            for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
#ifndef NDEBUG
        	printf("Greedy find a k-plex of size: %u\n", R_end);
#endif
        	return ;
        }
        else if(R_end == best_solution_size+1) return ;

#ifndef NDEBUG
#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = 0;i < R_end;i ++) for(ui j = i+1;j < R_end;j ++) {
        	ui v = SR[i], w = SR[j];
        	ui common_neighbors = 0;
        	for(ui k = S_end;k < R_end;k ++) if(matrix[SR[k]*n + v]&&matrix[SR[k]*n + w]) ++ common_neighbors;
        	assert(cn[v*n + w] == common_neighbors);
        	assert(cn[w*n + v] == common_neighbors);
        }
#endif
        for(ui i = 0;i < R_end;i ++) {
        	ui d1 = 0, d2 = 0;
        	for(ui j = 0;j < S_end;j ++) if(matrix[SR[i]*n + SR[j]]) ++ d1;
        	d2 = d1;
        	for(ui j = S_end;j < R_end;j ++) if(matrix[SR[i]*n + SR[j]]) ++ d2;
        	assert(d1 == degree_in_S[SR[i]]);
        	assert(d2 == degree[SR[i]]);
        }
#endif

        bool must_include = false;
        ui u = n; // u is the branching vertex
        if(choose_zero) {
        	assert(S_end == 0&&must_include_vertices_n == 0);
        	if(SR_rid[0] >= R_end) return ;
        	u = 0;
        	must_include = true;
        }
        else if(must_include_vertices_n) {
        	must_include = true;
        	u = must_include_vertices[-- must_include_vertices_n];
        	assert(SR_rid[u] >= S_end);
        	if(SR_rid[u] >= R_end) return ;
        }
        else {
        	for(ui i = S_end;i < R_end;i ++) if(degree[SR[i]] + 2 >= R_end) {
        		u = SR[i];
        		must_include = true;
        		break;
        	}
        	if(u == n) {
        		// u = choose_branch_vertex(S_end, R_end);
        		u = choose_branch_vertex_based_on_non_neighbors(S_end, R_end);
        	}
        }
        assert(degree[u] + K > best_solution_size&&degree[u] + K > S_end);

        // the first branch includes u into S
        assert(SR[SR_rid[u]] == u&&SR_rid[u] >= S_end&&SR_rid[u] < R_end);
        swap_pos(S_end, SR_rid[u]);
        ++ S_end;

        ui pre_best_solution_size = best_solution_size, old_R_end = R_end;
        ui old_removed_edges_n = 0;
#ifdef  _SECOND_ORDER_PRUNING_
        old_removed_edges_n = removed_edges_n;
#endif
        if(!move_u_to_S_with_prune(S_end, R_end, level)) BB_search(S_end, R_end, level+1, false, must_include_vertices_n);
        restore_SR_and_edges(S_end, R_end, old_R_end, level, old_removed_edges_n);
#ifdef  _SECOND_ORDER_PRUNING_
        assert(removed_edges_n == old_removed_edges_n);
#endif

        if(must_include) {
        	move_u_to_R_wo_prune(S_end, R_end, level);
        	return ;
        }

        // the second branch exclude u from S
        assert(Qv.empty());
#ifdef _SECOND_ORDER_PRUNING_
        while(!Qe.empty()) Qe.pop();
#endif
        bool pruned = remove_u_from_S_with_prune(S_end, R_end, level);
        if(!pruned&&best_solution_size > pre_best_solution_size) pruned = collect_removable_vertices_and_edges(S_end, R_end, level);
        if(!pruned) {
            if(!remove_vertices_and_edges_with_prune(S_end, R_end, level)) {
            	assert(must_include_vertices_n == 0);
            	if(degree[u] + K+1 == R_end) {
            		char *t_matrix = matrix + u*n;
            		bool reduction = true;
            		for(ui i = 0;i < R_end;i ++) if(SR[i] != u&&!t_matrix[SR[i]]) {
            			if(i >= S_end) must_include_vertices[must_include_vertices_n++] = SR[i];
            			if(degree[SR[i]] + K < R_end) {
            				reduction = false;
            				break;
            			}
            		}
            		if(!reduction) must_include_vertices_n = 0;
            	}
            	BB_search(S_end, R_end, level+1, false, must_include_vertices_n);
            }
        }
        restore_SR_and_edges(S_end, R_end, old_R_end, level, old_removed_edges_n);
    }

    bool move_u_to_S_with_prune(ui S_end, ui &R_end, ui level) {
    	assert(S_end > 0);
        ui u = SR[S_end-1];
        char *t_matrix = matrix + u*n;

        ui neighbors_n = 0, nonneighbors_n = 0;
        for(ui i = 0;i < R_end;i ++) if(i != S_end-1) {
            if(t_matrix[SR[i]]) neighbors[neighbors_n++] = SR[i];
            else nonneighbors[nonneighbors_n++] = SR[i];
        }

        for(ui i = 0;i < neighbors_n;i ++) ++ degree_in_S[neighbors[i]];

        assert(Qv.empty());
        if(degree_in_S[u] + K == S_end) { // only neighbors of u in R can be candidates --- RR2
        	ui i = 0;
        	while(i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end) ++ i;
            for(;i < nonneighbors_n;i ++) { // remove non-neighbors from R
            	assert(level_id[nonneighbors[i]] > level);
            	level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }
        else { // only non-neighbors of u may change their allowance --- RR1
        	ui i = 0;
        	while(i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end) ++ i;
            for(;i < nonneighbors_n;i ++) if(S_end - degree_in_S[nonneighbors[i]] >= K) {
            	assert(level_id[nonneighbors[i]] > level);
            	level_id[nonneighbors[i]] = level;
                Qv.push(nonneighbors[i]);
            }
        }

        // RR2
        for(ui i = 0;i < nonneighbors_n&&SR_rid[nonneighbors[i]] < S_end;i ++) if(degree_in_S[nonneighbors[i]] + K == S_end) {
            char *tt_matrix = matrix + nonneighbors[i]*n;
            for(ui j = S_end;j < R_end;j ++) if(level_id[SR[j]] > level&&!tt_matrix[SR[j]]) {
            	level_id[SR[j]] = level;
                Qv.push(SR[j]);
            }
        }

#ifdef _SECOND_ORDER_PRUNING_
        // RR4
        for(ui i = 0;i < nonneighbors_n;i ++) {
            int v = nonneighbors[i];
            assert(!t_matrix[v]);
            if(SR_rid[v] < S_end||level_id[v] == level||t_matrix[v]) continue;
            if(upper_bound_based_prune(S_end, u, v)) {
            	level_id[v] = level;
                Qv.push(v);
            }
        }

        // update cn(.,.)
        for(ui i = 0;i < neighbors_n;i ++) { // process common neighbors of u
            for(ui j = i+1;j < neighbors_n;j ++) {
#ifndef NDEBUG
            	if(!cn[neighbors[i]*n + neighbors[j]]) {
            		printf("cn[neighbors[i]*n + neighbors[j]]: %u %u\n", cn[neighbors[i]*n + neighbors[j]], cn[neighbors[j]*n + neighbors[i]]);
            	}
#endif
                assert(cn[neighbors[i]*n + neighbors[j]]);
                -- cn[neighbors[i]*n + neighbors[j]];
                -- cn[neighbors[j]*n + neighbors[i]];
            }
        }

        while(!Qe.empty()) Qe.pop();
        int new_n = 0;
        for(ui i = 0;i < nonneighbors_n;i ++) if(level_id[nonneighbors[i]] > level) nonneighbors[new_n ++] = nonneighbors[i];
        nonneighbors_n = new_n;
        for(ui i = 1;i < nonneighbors_n;i ++) { // process common non-neighbors of u
            ui w = nonneighbors[i];
            for(ui j = 0;j < i;j ++) {
            	ui v = nonneighbors[j];
            	if(!upper_bound_based_prune(S_end, v, w)) continue;
            	if(SR_rid[w] < S_end) return true; // v, w \in S --- UB2
            	else if(SR_rid[v] >= S_end) { // v, w, \in R --- RR5
            		if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
            	}
            	else { // RR4
            		assert(level_id[w] > level);
            		level_id[w] = level;
            		Qv.push(w);
            		break;
            	}
            }
        }
#endif

        return remove_vertices_and_edges_with_prune(S_end, R_end, level);
    }

    bool remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui level) {
#ifdef _SECOND_ORDER_PRUNING_
        while(!Qv.empty()||!Qe.empty()) {
#else
        while(!Qv.empty()) {
#endif
        	while(!Qv.empty()) {
				ui u = Qv.front(); Qv.pop(); // remove u
				assert(SR[SR_rid[u]] == u);
				assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
				-- R_end;
				swap_pos(SR_rid[u], R_end);

				bool terminate = false;
				ui neighbors_n = 0;
				char *t_matrix = matrix + u*n;
				for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
					ui w = SR[i];
					neighbors[neighbors_n++] = w;
					-- degree[w];
					if(degree[w] + K <= best_solution_size) {
						if(i < S_end) terminate = true; // UB1
						else if(level_id[w] > level) { // RR3
							level_id[w] = level;
							Qv.push(w);
						}
					}
				}
				// UB1
				if(terminate) {
					for(ui i = 0;i < neighbors_n;i ++) ++ degree[neighbors[i]];
					level_id[u] = n;
					++ R_end;
					return true;
				}

#ifdef _SECOND_ORDER_PRUNING_
				for(ui i = 1;i < neighbors_n;i ++) {
					ui w = neighbors[i];
					for(ui j = 0;j < i;j ++) {
						ui v = neighbors[j];
						assert(cn[v*n+w]);
#ifndef NDEBUG
						ui common_neighbors = 0;
						for(ui k = S_end;k <= R_end;k ++) if(matrix[SR[k]*n + v]&&matrix[SR[k]*n + w]) ++ common_neighbors;
						assert(cn[v*n + w] == common_neighbors);
						assert(cn[w*n + v] == common_neighbors);
#endif
						-- cn[v*n + w];
						-- cn[w*n + v];
#ifndef NDEBUG
						common_neighbors = 0;
						for(ui k = S_end;k < R_end;k ++) if(matrix[SR[k]*n + v]&&matrix[SR[k]*n + w]) ++ common_neighbors;
						assert(cn[v*n + w] == common_neighbors);
						assert(cn[w*n + v] == common_neighbors);
#endif

						if(!upper_bound_based_prune(S_end, v, w)) continue;

						if(SR_rid[w] < S_end) terminate = true; // v, w \in S --- UB2
						else if(SR_rid[v] >= S_end) { // v, w, \in R --- RR5
							if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
						}
						else if(level_id[w] > level) { // RR4
							level_id[w] = level;
							Qv.push(w);
						}
					}
				}
				if(terminate) return true;
#endif
        	}

#ifdef _SECOND_ORDER_PRUNING_
        	if(Qe.empty()) break;

        	ui v = Qe.front().first, w =  Qe.front().second; Qe.pop();
        	if(level_id[v] <= level||level_id[w] <= level||!matrix[v*n + w]) continue;
        	assert(SR_rid[v] >= S_end&&SR_rid[v] < R_end&&SR_rid[w] >= S_end&&SR_rid[w] < R_end);

        	if(degree[v] + K <= best_solution_size + 1) {
        		level_id[v] = level;
        		Qv.push(v);
        	}
        	if(degree[w] + K <= best_solution_size + 1) {
        		level_id[w] = level;
        		Qv.push(w);
        	}
        	if(!Qv.empty()) continue;

#ifndef NDEBUG
        	//printf("remove edge between %u and %u\n", v, w);
#endif

        	assert(matrix[v*n + w]);
        	matrix[v*n + w] = matrix[w*n + v] = 0;
        	-- degree[v]; -- degree[w];

        	if(removed_edges.size() == removed_edges_n) {
        		removed_edges.push_back(std::make_pair(v,w));
        		++ removed_edges_n;
        	}
        	else removed_edges[removed_edges_n ++] = std::make_pair(v,w);

        	char *t_matrix = matrix + v*n;
        	for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
        		-- cn[w*n + SR[i]];
        		-- cn[SR[i]*n + w];
        		if(!upper_bound_based_prune(S_end, w, SR[i])) continue;
        		if(i < S_end) {
        			if(level_id[w] > level) {
        				level_id[w] = level;
        				Qv.push(w);
        			}
        		}
        		else if(matrix[w*n + SR[i]]) Qe.push(std::make_pair(w, SR[i]));
        	}
        	t_matrix = matrix + w*n;
        	for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
        		-- cn[v*n + SR[i]];
        		-- cn[SR[i]*n + v];
        		if(!upper_bound_based_prune(S_end, v, SR[i])) continue;
        		if(i < S_end) {
        			if(level_id[v] > level) {
        				level_id[v] = level;
        				Qv.push(v);
        			}
        		}
        		else if(matrix[v*n + SR[i]]) Qe.push(std::make_pair(v, SR[i]));
        	}
#endif
        }

        return false;
    }

    void restore_SR_and_edges(ui S_end, ui &R_end, ui old_R_end, ui level, ui old_removed_edges_n) {
        while(!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(level_id[u] == level&&SR_rid[u] < R_end);
            level_id[u] = n;
        }
        while(R_end < old_R_end) { // insert u back into R
            ui u = SR[R_end];
            assert(level_id[u] == level&&SR_rid[u] == R_end);
            level_id[u] = n;

            ui neighbors_n = 0;
            char *t_matrix = matrix + u*n;
            for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
                ui w = SR[i];
                neighbors[neighbors_n ++] = w;
                ++ degree[w];
            }
#ifdef _SECOND_ORDER_PRUNING_
            for(ui i = 0;i < neighbors_n;i ++) {
                ui v = neighbors[i];
                for(ui j = i + 1;j < neighbors_n;j ++) {
                    ui w = neighbors[j];
                    ++ cn[v*n + w];
                    ++ cn[w*n + v];
                }
            }
            ui *t_cn = cn + u*n;
            for(ui i = 0;i < R_end;i ++) t_cn[SR[i]] = 0;
            for(ui i = 0;i < neighbors_n;i ++) if(SR_rid[neighbors[i]] >= S_end) {
            	ui v = neighbors[i];
            	char *t_matrix = matrix + v*n;
            	for(ui j = 0;j < R_end;j ++) if(t_matrix[SR[j]]) ++ t_cn[SR[j]];
            }
            for(ui i = 0;i < R_end;i ++) {
            	cn[SR[i]*n + u] = t_cn[SR[i]];
#ifndef NDEBUG
            	ui common_neighbors = 0, v = SR[i], w = u;
            	for(ui k = S_end;k < R_end;k ++) if(matrix[SR[k]*n + v]&&matrix[SR[k]*n + w]) ++ common_neighbors;
            	if(t_cn[SR[i]] != common_neighbors) printf("t_cn[SR[i]] = %u, comon_neighbors = %u\n", t_cn[SR[i]], common_neighbors);
            	assert(t_cn[SR[i]] == common_neighbors);
#endif
            }
#endif

            ++ R_end;
        }

#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = old_removed_edges_n;i < removed_edges_n; i ++) { // insert edge back into matrix
        	ui v = removed_edges[i].first, w = removed_edges[i].second;
        	assert(SR_rid[v] >= S_end&&SR_rid[v] < R_end&&SR_rid[w] >= S_end&&SR_rid[w] < R_end);
        	if(matrix[v*n + w]) continue;

#ifndef NDEBUG
        	//printf("restore edge between %u and %u\n", v, w);
#endif
			matrix[v*n + w] = matrix[w*n + v] = 1;
			++ degree[v]; ++ degree[w];

			char *t_matrix = matrix + v*n;
			for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
				++ cn[w*n + SR[i]];
				++ cn[SR[i]*n + w];
			}
			t_matrix = matrix + w*n;
			for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
				++ cn[v*n + SR[i]];
				++ cn[SR[i]*n + v];
			}
        }
        removed_edges_n = old_removed_edges_n;
#endif
    }
    
    void move_u_to_R_wo_prune(ui &S_end, ui &R_end, ui level) {
    	assert(S_end);
        ui u = SR[-- S_end];
        ui neighbors_n = 0;
        char *t_matrix = matrix + u*n;
        for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) neighbors[neighbors_n ++] = SR[i];
        for(ui i = 0;i < neighbors_n;i ++) -- degree_in_S[neighbors[i]];

#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = 0;i < neighbors_n;i ++) {
        	ui v = neighbors[i];
        	for(ui j = i+1;j < neighbors_n;j ++) {
        		ui w = neighbors[j];
        		++ cn[v*n + w];
        		++ cn[w*n + v];
        	}
        }
#endif
    }

    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui level) {
    	assert(S_end);
		ui u = SR[S_end-1];
		-- S_end; -- R_end;
		swap_pos(S_end, R_end);
		level_id[u] = level;

		bool ret = false;
        ui neighbors_n = 0;
        char *t_matrix = matrix + u*n;
        for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) neighbors[neighbors_n ++] = SR[i];
        for(ui i = 0;i < neighbors_n;i ++) {
        	-- degree_in_S[neighbors[i]];
        	-- degree[neighbors[i]];
        	if(degree[neighbors[i]] + K <= best_solution_size) {
        		if(SR_rid[neighbors[i]] < S_end) ret =  true;
        		else {
        			assert(level_id[neighbors[i]] > level);
        			level_id[neighbors[i]] = level;
        			Qv.push(neighbors[i]);
        		}
        	}
        }
        if(ret) return true;

#ifdef _SECOND_ORDER_PRUNING_
        for(ui i = 1;i < neighbors_n;i ++) if(level_id[neighbors[i]] > level) {
        	ui w = neighbors[i];
        	for(ui j = 0;j < i;j ++) {
        		ui v = neighbors[j];
        		if(!upper_bound_based_prune(S_end, v, w)) continue;

        		if(SR_rid[w] < S_end) return true; // v, w \in S
				else if(SR_rid[v] >= S_end) { // v, w, \in R
					if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
				}
				else {
					assert(level_id[w] > level);
					level_id[w] = level;
					Qv.push(w);
					break;
				}
        	}
        }
#endif
		return false;
	}

    bool collect_removable_vertices_and_edges(ui S_end, ui R_end, ui level) {
    	for(ui i = 0;i < S_end;i ++) if(degree[SR[i]] + K <= best_solution_size) return true;

#ifdef _SECOND_ORDER_PRUNING_
    	for(ui i = 0;i < S_end;i ++) for(ui j = i+1;j < S_end;j ++) {
    		if(upper_bound_based_prune(S_end, SR[i], SR[j])) return true;
    	}
#endif

    	for(ui i = S_end;i < R_end;i ++) if(level_id[SR[i]] > level){
    		if(S_end - degree_in_S[SR[i]] >= K||degree[SR[i]] + K <= best_solution_size) {
    			assert(level_id[SR[i]] > level);
    			level_id[SR[i]] = level;
    			Qv.push(SR[i]);
    			continue;
    		}
    		char *t_matrix = matrix + SR[i]*n;
    		for(ui j = 0;j < S_end;j ++) {
#ifdef _SECOND_ORDER_PRUNING_
    			if((S_end - degree_in_S[SR[j]] == K&&!t_matrix[SR[j]])||upper_bound_based_prune(S_end, SR[i], SR[j]))
#else
    			if(S_end - degree_in_S[SR[j]] == K&&!t_matrix[SR[j]])
#endif
    			{
    				assert(level_id[SR[i]] > level);
    				level_id[SR[i]] = level;
    				Qv.push(SR[i]);
    				break;
    			}
    		}
    	}

#ifdef _SECOND_ORDER_PRUNING_
    	for(ui i = S_end;i < R_end;i ++) if(level_id[SR[i]] > level) {
    		for(ui j = i+1;j < R_end;j ++) if(level_id[SR[j]] > level&&matrix[SR[i]*n + SR[j]]) {
    			if(upper_bound_based_prune(S_end, SR[i], SR[j])) Qe.push(std::make_pair(SR[i], SR[j]));
    		}
    	}
#endif

        return false;
    }

#ifdef _SECOND_ORDER_PRUNING_
    bool upper_bound_based_prune(ui S_end, ui u, ui v) {
    	// ui ub = S_end + 2*K - (S_end - degree_in_S[u]) - (S_end - degree_in_S[v]) + cn[u*n + v];
    	ui ub = 2*K + degree_in_S[u] - S_end + degree_in_S[v] + cn[u*n + v];
    	if(SR_rid[u] >= S_end) {
    		-- ub; // S_end ++
    		if(matrix[u*n+v]) ++ ub; // degree_in_S[v] ++
    	}
    	if(SR_rid[v] >= S_end) {
    		-- ub;
    		if(matrix[v*n+u]) ++ ub;
    	}
    	return ub <= best_solution_size;
    }
#endif

    bool inRange(ui u, ui a, ui b) {
        return(SR_rid[u] >= a&&SR_rid[u] < b);
    }

    void swap_pos(ui i, ui j) {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex_based_on_non_neighbors(ui S_end, ui R_end) {
    	ui u = n, min_degree_in_S = n;
    	for(ui i = 0;i < R_end;i ++) {
    		ui v = SR[i];
    		if(degree[v] + K >= R_end) continue;
    		if(degree_in_S[v] < min_degree_in_S) {
    			u = v;
    			min_degree_in_S = degree_in_S[v];
    		}
    		else if(degree_in_S[v] == min_degree_in_S) {
    			if((i < S_end&&degree[v] < degree[u])||(i >= S_end&&degree[v] > degree[u])) {
    				u = v;
    				min_degree_in_S = degree_in_S[v];
    			}
    		}
    	}
    	assert(u != n);
    	if(SR_rid[u] < S_end) {
    		char *t_matrix = matrix+u*n;
    		u = n;
    		ui max_degree = 0;
    		for(ui i = S_end;i < R_end;i ++) if(!t_matrix[SR[i]]&&degree[SR[i]] > max_degree) {
    			max_degree = degree[SR[i]];
    			u = SR[i];
    		}
    	}
    	assert(u != n);
    	return u;
    }
};

#endif
