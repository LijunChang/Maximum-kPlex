/************************************************

** This is an exact solver for maximum k-plexes in massive graphs.
** Date:	2017.12.31
** Revision 2018.06.16
** Author:	Jiejiang chen, Jian Gao

** To compile the solver use command: g++ -std=c++11 -O3 kplex.cpp -o kplex

vertex[i].neighbor[] stores the neighbors, it dynamically changes after removing vertices

U stores the vertices whose neighborhood does not change after a reduction

************************************************/


#include <iostream>
#include <cstring>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <queue>
#include <vector>
#include <algorithm>
#include<unordered_set>
#include<sys/time.h>
#include<sys/types.h>
#include<sys/resource.h>
#include<limits.h>
#include "Array.h"

#define insert_v(end, value) *(end++) = value;
#define delete_i(index, end) *index = *(--end);

class KPlex_IJCAI18 {
public:

struct Vertex {
    int state, degree;
    int *neighbor;
};

int v_n = 0;//number of vertex
int e_n = 0;//number of edge

int *neighbor_len = NULL;//
int **neighbor = NULL;//

int *temp_array = NULL;//
int *temp_mark = NULL;//
int *temp_index = NULL;//
int *neighbor_in_solution = NULL;//
int *temp_remove = NULL;
int temp_remove_size;
Vertex *vertex = NULL;//

Array *remaining_vertex = NULL;//
int para_k = 1;//

Array *U = NULL;//
Array *S = NULL;//

int *best_solution = NULL;//
int best_solution_size = 0;//

int max_degree = -1;//

double read_time, init_time, search_time;//
int *t_U = NULL;//

std::unordered_set<unsigned long long> edge_hash;//

bool output;
int *edge_buf1 = NULL;
int *edge_buf2 = NULL;

KPlex_IJCAI18(bool output_) {
	output = output_;
}

~KPlex_IJCAI18() {
    free(neighbor_len);
    free(temp_array);
    free(temp_mark);
    free(temp_index);
    free(neighbor_in_solution);
    free(best_solution);
    free(temp_remove);
    free(t_U);
    for(int i = 0; i < v_n; i++) {
        vertex[i].neighbor = NULL;
        neighbor[i] = NULL;
    }
	free(edge_buf1);
	free(edge_buf2);
    free(vertex);
    free(neighbor);

    free(remaining_vertex);
    free(S);
    free(U);
}

double get_utime() {
    struct rusage utime;
    getrusage(RUSAGE_SELF, &utime);
    return (double) (utime.ru_utime.tv_sec + (double)utime.ru_utime.tv_usec / 1000000);
}

unsigned long long encode_pairID(unsigned long long v1, unsigned long long v2) {
    unsigned long long n1, n2;
    unsigned long long pairID;
    if(v1 < v2) {
        n1 = v1; n2 = v2;
    }
    else {
        n1 = v2; n2 = v1;
    }
    pairID = ((n1 + n2 + 1) * (n1 + n2) >> 1) + n2;
    return pairID;
}

int edge_is(int m, int n) {
    unsigned long long id = encode_pairID(m, n);
    if(edge_hash.count(id)) return 1;
    else return 0;
}

void allocateMemory(int v_n, int e_n) {
	if(v_n <= 0||e_n <= 0) {
		//printf("v_n, e_n are not correct!\n");
		//exit(1);
        return ;
	}

	neighbor = (int **)malloc(v_n * sizeof(int**));
    neighbor_len = (int *)malloc(v_n * sizeof(int));
    edge_buf1 = (int *)malloc(2 * e_n * sizeof(int));
    edge_buf2 = (int *)malloc(2 * e_n * sizeof(int));
    temp_array = (int *)malloc(v_n * sizeof(int));
    temp_mark = (int *)malloc(v_n * sizeof(int));
    temp_index = (int *)malloc(v_n * sizeof(int));
    temp_remove = (int *)malloc(v_n * sizeof(int));
    vertex = (Vertex *)malloc(v_n * sizeof(Vertex));
    neighbor_in_solution = (int *)malloc(v_n * sizeof(int));
    t_U = (int *)malloc(v_n * sizeof(int));
    best_solution = (int *)malloc(v_n *sizeof(int));

    remaining_vertex = new Array(v_n);
    S = new Array(v_n);
    U = new Array(v_n);
}

void readGraph(char* File_Name) {
	std::ifstream FIC;
    FIC.open(File_Name);
    if(FIC.fail()) printf("### Error open, File_Name %s\n", File_Name);
    FIC >> v_n >> e_n;
    
    int a, b;
    unsigned long long id;

	std::vector<std::pair<int,int> > vp;
	vp.reserve(e_n);

    allocateMemory(v_n, e_n);
	edge_hash.clear();

    for(int i = 0; i < v_n; i++) {
         vertex[i].degree = 0;
         vertex[i].state = 0;
    }

    while(FIC >> a >> b) {
            vertex[a].degree++;
            vertex[b].degree++;
			vp.push_back(std::make_pair(a,b));
            id = encode_pairID(a, b);
            edge_hash.insert(id);
    }
    FIC.close();

    organize_edges(vp);
}

void load_graph(int v_n_, const std::vector<std::pair<int,int> > &vp) {
    v_n = v_n_;
    e_n = vp.size();
    for(int i = 0; i < v_n; i++) {
         vertex[i].degree = 0;
         vertex[i].state = 0;
    }
    edge_hash.clear();
    for(int i = 0;i < e_n;i ++) {
        int a = vp[i].first;
        int b = vp[i].second;
        assert(a >= 0&&a < v_n&&b >= 0&&b < v_n);
        vertex[a].degree ++;
        vertex[b].degree ++;
        unsigned long long id = encode_pairID(a,b);
        edge_hash.insert(id);
    }
    organize_edges(vp);
}

void organize_edges(const std::vector<std::pair<int,int> > &vp) {
    remaining_vertex->clear();
    S->clear();
    U->clear();

	max_degree = 0;
	int count = 0;
    for(int i = 0; i < v_n; i++) {
        if(vertex[i].degree > max_degree) max_degree = vertex[i].degree;
        vertex[i].neighbor = edge_buf1+count;
        neighbor[i] = edge_buf2+count;
		count += vertex[i].degree;
        remaining_vertex->insert_element(i);
    }

    for(int i = 0; i < v_n; i++) vertex[i].degree = 0;
    for(int e  = 0; e < e_n; e++) {
        int a = vp[e].first;
        int b = vp[e].second;
        vertex[a].neighbor[vertex[a].degree] = b;
        vertex[b].neighbor[vertex[b].degree] = a;
        vertex[a].degree++;
        vertex[b].degree++;
    }
    if(output) printf("### the number of vertex:\t%d\n### the number of edge:\t\t%d\n### maximum degree:\t\t%d\n", v_n, e_n, max_degree);
}

void readGraph_binary(char* dir) {
	FILE *f = fopen( (std::string(dir) + std::string("/b_degree.bin")).c_str(), "rb");
	if(f == NULL) {
		printf("!!! Cannot open file %s\n", (std::string(dir) + std::string("/b_degree.bin")).c_str());
		return ;
	}

	int tt;
	fread(&tt, sizeof(int), 1, f);
	if(tt != sizeof(int)) {
		printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
		return ;
	}
	fread(&v_n, sizeof(int), 1, f);
	fread(&e_n, sizeof(int), 1, f);
	e_n /= 2;
	// printf("%d %d\n", v_n, e_n);

	int *degree = new int[v_n];
	fread(degree, sizeof(int), v_n, f);
	fclose(f);

    std::vector<std::pair<int,int> > vp;
    vp.reserve(e_n);

    int a,b;
    unsigned long long id;

	f = fopen( (std::string(dir) + std::string("/b_adj.bin")).c_str(), "rb");
	if(f == NULL) {
		printf("!!! Cannot open file %s\n", (std::string(dir) + std::string("/b_adj.bin")).c_str());
		return ;
	}

    allocateMemory(v_n,e_n);
    edge_hash.clear();

    for(int i = 0; i < v_n; i++) {
         vertex[i].degree = 0;
         vertex[i].state = 0;
    }

	int *buf = new int[v_n];
	for(int i = 0;i < v_n;i ++) if(degree[i] > 0) {
		fread(buf, sizeof(int), degree[i], f);
		for(int j = 0;j < degree[i];j ++) if(buf[j] > i) {
			a = i; b = buf[j];
			vertex[a].degree++;
            vertex[b].degree++;
            vp.push_back(std::make_pair(a,b));
            id = encode_pairID(a, b);
            edge_hash.insert(id);
		}
	}
	fclose(f);

    assert(vp.size() == e_n);
    if(vp.size() != e_n) printf("vp.size != e_n\n");

	delete[] buf;
	delete[] degree;

    organize_edges(vp);
}

void update_best_solution() {
    if(output) {
        printf("%d %d ", best_solution_size, S->size());
        if(S->is_in_array(0)) printf("contain 0\n");
        else printf("does not contain 0\n");
    }
    assert(S->size() > best_solution_size);
    best_solution_size = S->size();
    for(int i = S->begin(); i < S->size(); i++)
        best_solution[i] = S->element_at(i);
}

bool check_solution()//
{
    int i, j, adj_num, v, u;
    for(i = 0; i < best_solution_size; i++)
    {
        adj_num = 0;
        v = best_solution[i];
        for(j = 0; j < best_solution_size; j++)
        {
            u = best_solution[j];
            if(v ==u)
                continue;
            if(edge_is(v,u))
                adj_num++;
        }
        if(adj_num < best_solution_size - para_k)
            return false;
    }
    return true;
}

void printf_solution()//
{
    if(check_solution())
    {
        if(output) printf("after checking , the solution is correct,  solution size: %d, time: %f\n", best_solution_size, get_utime());
    }
    else
    {
        printf("the solution found is wrong\n");
    }

}

void reduce_graph_1()// vertex-reduction
{
    int *degree_counter, *where, *candidate;
    int node, p1, i, j, h, k, t, neighbors,_tsize;
    int u,v;
    int cur_degree = 0;
    degree_counter = temp_mark;
    where = temp_index;
    candidate = temp_array;

    for(i = 0; i <= max_degree; i++)
        degree_counter[i] = 0;
    for(node = remaining_vertex->begin(); node < remaining_vertex->size(); node++)
    {
        v = remaining_vertex->element_at(node);
        degree_counter[vertex[v].degree]++;
    }
    j = 0;
    for(i = 0; i <= max_degree; i++)
    {
        k = degree_counter[i];
        degree_counter[i] = j;
        j += k;
    }

    for(node = remaining_vertex->begin(); node < remaining_vertex->size(); node++)
    {
        v = remaining_vertex->element_at(node);
        t = degree_counter[vertex[v].degree];
        degree_counter[vertex[v].degree]++;
        candidate[t] = v;
        where[v]  = t;
    }

    for(i = max_degree; i > 0; i--)
        degree_counter[i] = degree_counter[i - 1];

    degree_counter[0] = 0;
    p1 = 0;
    cur_degree = vertex[candidate[0]].degree;

    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        neighbor_len[v] = vertex[v].degree;
    }

    while(p1 < remaining_vertex->size())
    {
        node = candidate[p1];
        if(p1 < remaining_vertex->size() - 1  && neighbor_len[node] == neighbor_len[candidate[p1+1]])
            degree_counter[neighbor_len[node]] = p1 + 1;
        if(neighbor_len[node] + para_k - 1 >= remaining_vertex->size() - p1 - 1)
        {
            S->clear();
            for(i = p1; i < remaining_vertex->size(); i++)
                S->insert_element(candidate[i]);
            break;
        }

        for(i = 0; i < vertex[node].degree; i++)
        {
            neighbors = vertex[node].neighbor[i];
            if(where[neighbors] > p1)
            {
                t = where[neighbors];
                h = degree_counter[neighbor_len[neighbors]];
                k = candidate[h];

                candidate[h] = neighbors;
                where[neighbors] = h;
                candidate[t] = k;
                where[k]= t;

                degree_counter[neighbor_len[neighbors]]++;
                neighbor_len[neighbors]--;
                if(neighbor_len[neighbors] != neighbor_len[candidate[h-1]])
                    degree_counter[neighbor_len[neighbors]] = h;
            }
        }
        p1++;
    }
	// the above heuristically computes a k-plex based on the degeneracy ordering

    if(S->size() > best_solution_size)
        update_best_solution();

	std::queue<int> remove_que;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        neighbor_len[v] = vertex[v].degree;
        if(vertex[v].degree + para_k <= best_solution_size ||  remaining_vertex->size() <= best_solution_size)
        {
            vertex[v].state = -1;
            remove_que.push(v);
        }
    }

    while(!remove_que.empty())//
    {
            v = remove_que.front();
            remove_que.pop();
            remaining_vertex->delete_element(v);
            for(i = 0; i < vertex[v].degree; i++)
            {
                u = vertex[v].neighbor[i];
                neighbor_len[u]--;
                if(!vertex[u].state && (neighbor_len[u] + para_k <= best_solution_size ||  remaining_vertex->size() <= best_solution_size))
                {
                    vertex[u].state = -1;
                    remove_que.push(u);
                }
            }
    }

    int mark = 1;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)// reorganize neighbors
    {
        v = remaining_vertex->element_at(i);
        mark = 1;
        for(j = 0; j < vertex[v].degree; j++)
        {
            u = vertex[v].neighbor[j];
            if(vertex[u].state)
            {
                mark = 0;
                vertex[v].degree--;
                vertex[v].neighbor[j] = vertex[v].neighbor[vertex[v].degree];
                j--;
            }
        }
        if(!mark && U->is_in_array(v))// U is used to keep track the verticese whose neighbors does not change
            U->delete_element(v);
    }
}

void reduce_graph_2(Array *Uset)// subgraph-reduction
{
    int i,v,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int ssize;
    ssize = 0;//
	std::queue<int> remove_que;
    for(i = Uset->begin(); i < Uset->size(); i++)
    {
        u1 = Uset->element_at(i);
        neighbor_len[u1] = 0;
    }

	// printf("U.size(): %d\n", U->size());

    for(i = Uset->begin(); i < Uset->size(); i++)
    {
        v = Uset->element_at(i);//
        if(U->is_in_array(v))//
            continue;
        kv = para_k - 1;//
        rv = best_solution_size - kv - ssize;
        rd = vertex[v].degree;
        neighbor_len[v] = 0;
        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))
                continue;
            neighbor_len[u1] = 1;
            neighbor[u1][0] = v;
            //neighbor_len[v]++;
            neighbor[v][neighbor_len[v]++] = u1;
        }
        if(neighbor_len[v]  < rv)//
        {
            Uset->delete_element(v);
            i--;
            continue;
        }
        while(!remove_que.empty())//
            remove_que.pop();
        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))//
                continue;
            for(k = j + 1; k < rd; k++)
            {
                u2 = vertex[v].neighbor[k];
                if(!Uset->is_in_array(u2))
                    continue;
                if(edge_is(u1, u2))
                {
                    neighbor[u1][neighbor_len[u1]++] = u2;
                    neighbor[u2][neighbor_len[u2]++] = u1;
                }
            }
            ku = para_k - 1;
            cu = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
            temp_index[u1] = neighbor_len[u1];//
            if(cu < rv)//
            {
                temp_index[u1] = 0;//
                remove_que.push(u1);
            }
        }

        temp_index[v] = neighbor_len[v];//
        while(!remove_que.empty())//
        {
            u1 = remove_que.front();
            remove_que.pop();
            for(j = 0; j < neighbor_len[u1]; j++)
            {
                u2 = neighbor[u1][j];
                if(!temp_index[u2])//
                    continue;
                temp_index[u2]--;//
                if(temp_index[u2] + ku < rv && u2 != v)//
                {
                    temp_index[u2] = 0;//
                    remove_que.push(u2);
                }
            }
            if(temp_index[v] < rv)//
            {
                Uset->delete_element(v);
                i--;
                break;
            }
        }
    }

    U->clear();//
    int mark = 1;
    for(i = Uset->begin(); i < Uset->size(); i++)//
    {
        v = Uset->element_at(i);
        mark = 1;
        for(j = 0; j < vertex[v].degree; j++)
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))
            {
                mark = 0;
                vertex[v].degree--;
                vertex[v].neighbor[j] = vertex[v].neighbor[vertex[v].degree];
                j--;
            }
        }
        if(mark)//
            U->insert_element(v);
    }
}

void  reduce_graph_in_BB(int* &begin, int* &end, int rev)// subgraph-reduction in BB
{
    int i,v,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int mark = 1;
    int ssize;
    int *ii;
	std::queue<int> remove_que;
    ssize = S->size();//

    for(ii = begin; ii < end; ii++)
        temp_index[*ii] = 1;//

    for(i = 0; i < vertex[rev].degree; i++)
    {
        u1 = vertex[rev].neighbor[i];
        if(temp_index[u1])//
            temp_index[u1] = 2;
    }

    for(ii = begin; ii < end; ii++)
    {
        mark = 1;
        v  = *ii;
        if(temp_index[v] < 2)
            continue;
        kv = para_k - 1 - (ssize - neighbor_in_solution[v]);//
        rv = best_solution_size - kv - ssize;//
        rd = vertex[v].degree;
        neighbor_len[v] = 0;
        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!temp_index[u1])
                continue;
            neighbor[u1][0] = v;
            neighbor_len[u1] = 1;
            neighbor[v][neighbor_len[v]++] = u1;
        }

        if(neighbor_len[v]  < rv)//
        {
            delete_i(ii, end);
            ii--;
            temp_remove[temp_remove_size++] = v;
            temp_index[v] = 0;
          /*  for(i  = 0; i < neighbor_len[v]; i++)
            {
                u1 = neighbor[v][i];
                vertex[u1].state--;
            }*/
            continue;
        }

        while(!remove_que.empty())
            remove_que.pop();

        for(j = 0; j < rd; j++)//
        {
            u1 = vertex[v].neighbor[j];
            if(!temp_index[u1])
                continue;
            for(k = j + 1; k < rd; k++)
            {
                u2 = vertex[v].neighbor[k];
                if(!temp_index[u2])
                    continue;
                if(edge_is(u1, u2))
                {
                    neighbor[u1][neighbor_len[u1]++] = u2;
                    neighbor[u2][neighbor_len[u2]++] = u1;
                }
            }

            ku = para_k - 1 - (ssize - neighbor_in_solution[u1]);
            temp_array[u1] = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
            //temp_array[u1] = neighbor_len[u1] + ku;
            if(temp_array[u1] < rv )
            {
                temp_array[u1] = 0;//
                remove_que.push(u1);
            }
        }

        temp_array[v] = neighbor_len[v];//
        while(!remove_que.empty())//
        {
            u1 = remove_que.front();
            remove_que.pop();
            temp_array[v]--;//
            for(j = 1; j < neighbor_len[u1]; j++)//
            {
                u2 = neighbor[u1][j];
                if(!temp_array[u2])//
                    continue;
                temp_array[u2]--;//
                if(temp_array[u2] > temp_array[v])//
                    temp_array[u2] = temp_array[v];
                if(temp_array[u2] < rv)//
                {
                    temp_array[u2] =0;
                    remove_que.push(u2);
                }
            }
            if(temp_array[v] < rv)//
            {
                delete_i(ii, end);
                temp_remove[temp_remove_size++] = v;
                ii--;//
                temp_index[v] = 0;//
            /*  for(i = 0; i < neighbor_len[v]; i++)
                {
                    u1 = neighbor[v][i];
                    vertex[u1].state--;
                }*/
                //mark = 0;
                break;
            }
        }
    }

    for(ii = begin; ii < end; ii++)//
    {
        v = *ii;
        vertex[v].state = 0;
        for(i = 0; i < vertex[v].degree; i++)
        {
            if(temp_index[vertex[v].neighbor[i]])
                vertex[v].state++;
        }
    }
}

void reduce_graph() {
    int rn = remaining_vertex->size();
    if(output) printf("best solution\tcurrent solution\tremaining vertex\n");
    if(output) printf("%8d\t%8d\t\t%8d\n", 0, 0, v_n);
    while(rn > 0) {
        reduce_graph_1();
        if(output) printf("%8d\t%8d\t\t%8d(vertex-reduce)\n", best_solution_size, S->size(), remaining_vertex->size());
        reduce_graph_2(remaining_vertex);
        if(output) printf("%8d\t%8d\t\t%8d(subgraph-reduce)\n", best_solution_size, S->size(), remaining_vertex->size());

        if(remaining_vertex->size() < rn) rn = remaining_vertex->size();
        else break;
    }
}

int calculate_upbound_in_BB(int *begin, int *end)//
{
    int i = 0,vd,v,j;
    int *ii;
    int tt = 0;
    for(ii = begin; ii < end; ii++)
    {
        v = *ii;
        neighbor_len[v] = vertex[v].state;
        if(neighbor_len[v] + para_k + neighbor_in_solution[v] > best_solution_size)
            tt++;
        temp_index[v] = 0;
    }
    return tt;
}

void reduceU1(int* &tbegin, int* &tend) // vertex-reduction in BB
{
    int *ii, v, i, u;
	std::queue<int> remove_que;
    for(ii = tbegin; ii < tend; ii++)//
        temp_index[*ii] = 1;

    for(ii = tbegin; ii < tend; ii++)
    {
        v = *ii;
        vertex[v].state = 0;
        for(i = 0; i < vertex[v].degree; i++)
            if(temp_index[vertex[v].neighbor[i]] > 0)
                vertex[v].state++;//
        if(vertex[v].state + neighbor_in_solution[v] + para_k <= best_solution_size)//
        {
            remove_que.push(v);
            temp_index[v] = 2;
        }
    }

    while(!remove_que.empty())//
    {
        v = remove_que.front();
        remove_que.pop();
        temp_index[v] = 0;
        for(i = 0; i < vertex[v].degree; i++)
        {
            u = vertex[v].neighbor[i];
            if(temp_index[u] != 1)
                continue;
            vertex[u].state--;
            if(vertex[u].state + neighbor_in_solution[u] + para_k <= best_solution_size)
            {
                temp_index[u] = 0;
                remove_que.push(u);
            }
        }
    }

    for(ii = tbegin; ii < tend; ii++)//
    {
        v = *ii;
        if(!temp_index[v])
        {
            delete_i(ii, tend);
            ii--;
            temp_remove[temp_remove_size++] = v;
        }
    }
}

void reduceV(int* &begin, int* &end, int v)// v-reduction
{
    int i,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int ssize;
    int *ii;
	std::queue<int> remove_que;
    ssize = S->size();
    kv = para_k - (ssize - neighbor_in_solution[v]);//
    rv = best_solution_size - kv - ssize;//
    rd = vertex[v].degree;
    neighbor_len[v] = 0;
    for(j = 0; j < rd; j++)//
    {
        u1 = vertex[v].neighbor[j];
        if(!temp_index[u1])
            continue;
        neighbor[u1][0] = v;
        neighbor_len[u1] = 1;
        neighbor[v][neighbor_len[v]++] = u1;
    }
    for(j = 0; j < rd; j++)//
    {
        u1 = vertex[v].neighbor[j];
        if(!temp_index[u1])
            continue;
        for(k = j + 1; k < rd; k++)
        {
            u2 = vertex[v].neighbor[k];
            if(!temp_index[u2])
                continue;
            if(edge_is(u1, u2))
            {
                neighbor[u1][neighbor_len[u1]++] = u2;
                neighbor[u2][neighbor_len[u2]++] = u1;
            }
        }
        ku = para_k - 1 - (ssize - neighbor_in_solution[u1]);
        temp_array[u1] = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
        if(temp_array[u1] < rv)//
        {
            temp_array[u1] = 0;
            remove_que.push(u1);
            //cout << "%";
        }
    }
    temp_array[v] = neighbor_len[v];

    while(!remove_que.empty())//
    {
        u1 = remove_que.front();
        remove_que.pop();
        temp_array[v]--;//
        temp_index[u1] = 0;
        for(j = 0; j < vertex[u1].degree; j++)
        {
            u2 = vertex[u1].neighbor[j];//
            if(temp_index[u2])
                vertex[u2].state--;
        }
        for(j = 1; j < neighbor_len[u1]; j++)//
        {
            u2 = neighbor[u1][j];
            //vertex[u2].state--;
            if(!temp_array[u2])
                continue;
            temp_array[u2]--;
            if(temp_array[u2] > temp_array[v])
                temp_array[u2] = temp_array[v];
            if(temp_array[u2] < rv)
            {
                temp_array[u2] = 0;
                remove_que.push(u2);
            }
        }
    }

    for(ii = begin; ii < end; ii++)//
    {
        u1 = *ii;
        if(!temp_index[u1])
        {
            delete_i(ii, end);
            temp_remove[temp_remove_size++] = u1;
            ii--;
        }
    }
}

int is_S;
int nn= 0;
void update_candidate(int* &tbegin, int* &tend, int mark, int rev)//
{
    int tt1,tt2 = 0,tt3 = 0,tt4;
    int v, i, j;
    int *ii;
    if(mark)//
    {
        reduceU1(tbegin, tend);
        reduce_graph_in_BB(tbegin, tend, rev);
    }
    else//
    {
        for(ii = tbegin; ii < tend; ii++)
            temp_mark[*ii] = 0;
        is_S = 0;
        for(i = S->begin(); i < S->size(); i++)//
        {
            v = S->element_at(i);
            if(neighbor_in_solution[v] + para_k == S->size())
            {
                for(j = 0; j < vertex[v].degree; j++)
                        temp_mark[vertex[v].neighbor[j]]++;
                is_S++;
            }
        }
        for(ii = tbegin; ii < tend; ii++)
        {
            v = *ii;
            if((is_S && temp_mark[v] != is_S) || neighbor_in_solution[v] + para_k <= S->size())
            {
                delete_i(ii, tend);
                temp_remove[temp_remove_size++] = v;
                //temp_index[v] = 0;
                ii--;
            }
        }
        // the above is k-reduction

        reduceU1(tbegin, tend);
        reduceV(tbegin, tend, rev);
    }
}

int ct =  0;//

int BB(int* &tbegin, int* &tend, int *tv, int nt, bool fix_start_vertex)//
{
    int u, i, upper, lower, maxd, v, vd, d;
    int *maxv;
    int *ii;
    int *begin, *end;
    if(tbegin < tend)
    {
        ct++;
        nt++;
        v = *tv;
        delete_i(tv, tend);//
        S->insert_element(v);//
        for(i = 0; i < vertex[v].degree; i++)//
            neighbor_in_solution[vertex[v].neighbor[i]]++;
        begin = tbegin;
        end = tend;
        temp_remove_size = 0;
        update_candidate(tbegin, tend, 0, v);
        upper = calculate_upbound_in_BB(tbegin, tend) + S->size();

        if(end - tend != temp_remove_size)
        {
            printf("Error in add BB: %ld %d  ", end- begin, temp_remove_size);
        }

        begin = tend;
        end = begin;
        for(i = 0; i < temp_remove_size; i++)
            insert_v(end, temp_remove[i]);

        if(upper > best_solution_size)
        {
            maxd = -1;
			int minkv=para_k;
            //for(ii = begin; ii < end; ii++)//
            for(ii = tbegin; ii < tend; ii++)
            {
                u = *ii;
                //d = neighbor_len[u] + para_k - 1 - (S->size() - neighbor_in_solution[u]);//
                d = neighbor_len[u];//
                int kkk=para_k - 1 - (S->size() - neighbor_in_solution[u]);
				if(kkk==0 )
				{
					if(minkv==0 && d>maxd || minkv>0)
					{
					   maxv = ii;
					   minkv= kkk;
					   maxd = d;
					}

				}

				else if(minkv>0 && d > maxd)
                {

                    maxd = d;
                    maxv = ii;
					minkv=kkk;
                }

				/*if(d > maxd)
				{
					maxd = d;
                    maxv = ii;

				}*/

            }
            //cout << "a" << -(tbegin - tend) << " ";
            lower = BB(tbegin, tend, maxv, nt, fix_start_vertex);
            if(lower > best_solution_size)//
            {
                assert(lower == S->size());
                update_best_solution();
                printf_solution();
            }
        }

        for(ii = begin; ii < end; ii++)
            insert_v(tend, *ii);

        S->delete_element(v);//

        for(i = 0; i < vertex[v].degree; i++)//
            neighbor_in_solution[vertex[v].neighbor[i]]--;

        temp_remove_size = 0;
        begin = tbegin;
        end = tend;
        update_candidate(tbegin, tend, 1, v);
        upper = calculate_upbound_in_BB(tbegin, tend) + S->size();

        if(end - tend != temp_remove_size)
        {
            printf("Error in remove BB: %ld %d  ", end - begin, temp_remove_size);
        }

        begin = tend;
        end = begin;
        for(i = 0; i < temp_remove_size; i++)
            insert_v(end, temp_remove[i]);
        if(upper > best_solution_size&&(!fix_start_vertex || nt > 1))
        {
            maxd = -1;
			int minkv=para_k;
            for(ii = tbegin; ii < tend; ii++)//
            {
                u = *ii;
                //d = neighbor_len[u] + para_k - 1 - (S->size() - neighbor_in_solution[u]);//
                d = neighbor_len[u];//
                int kkk=para_k - 1 - (S->size() - neighbor_in_solution[u]);
				if(kkk==0 )
				{
					if(minkv==0 && d>maxd || minkv>0)
					{
						maxv = ii;
						minkv= kkk;
						maxd = d;
					}

				}

				else if(minkv>0 && d > maxd)
                {

                    maxd = d;
                    maxv = ii;
					minkv=kkk;
                }

				/*if(d > maxd)
				{
					maxd = d;
                    maxv = ii;

				}*/
            }
            //cout << "r" << -(tbegin - tend) << " ";
            lower = BB(tbegin, tend, maxv, nt, fix_start_vertex);
            if(lower > best_solution_size)
            {
                assert(lower == S->size());
                update_best_solution();
                printf_solution();
            }
        }
        for(ii = begin; ii < end; ii++)
            insert_v(tend, *ii);
        insert_v(tend, v);
        return best_solution_size;
    }
    else
    {
        //cout << nt << " ";
        return S->size();
    }

}

void search(bool fix_start_vertex = false) {
    int v, maxd = -1;
    int *begin, *end;
    int *maxv;
    begin = end = t_U;
    int en = 0;
    for(int i = remaining_vertex->begin(); i < remaining_vertex->size(); i++) {
        v = remaining_vertex->element_at(i);
        //U->insert_element(v);
        neighbor_in_solution[v] = 0;
        insert_v(end, v);//
        temp_index[v] = 0;
        en += vertex[v].degree; // only used for collecting statistics, can be deleted
        if(!fix_start_vertex&&vertex[v].degree > maxd) {
            maxd = vertex[v].degree;
            maxv = end - 1;
        }
		else if(fix_start_vertex&&v == 0) maxv = end-1;
    }
    S->clear();//
    BB(begin, end, maxv, 0, fix_start_vertex);
}

void kPlex(int K, std::vector<unsigned int> &kplex, bool must_include_0) {
    para_k = K;
    best_solution_size = kplex.size();
    reduce_graph();
    if(!must_include_0||remaining_vertex->is_in_array(0)) search(must_include_0);
    if(best_solution_size > kplex.size()) {
        kplex.clear();
        for(int i = 0;i < best_solution_size;i ++) kplex.push_back(best_solution[i]);
    }
}

int main(int argc, char *argv[])
{
	if(argc < 3) {
		printf("Usage: [1]exe [2]dir [3]k [4 optinal]time_limit\n");
		return 0;
	}
	printf("**** kplex build at %s %s ***\n", __TIME__, __DATE__);
    char* File_Name;
    File_Name = argv[1];//
    para_k = atoi(argv[2]);//
    printf("%s: \n", File_Name);
    readGraph(File_Name);//
    //readGraph_binary(File_Name);//
    read_time = get_utime();
    printf("----------------------------readGraph finish----------------------------\n");
    reduce_graph();//
    init_time = get_utime() - read_time;
    //printf_solution();
    search();//
    printf("read time: %f; init time: %f; search time: %f\n", read_time, init_time, get_utime() - read_time - init_time);
    printf("total time excluding read time: %f\n", get_utime() - read_time);
    //cout << ct << endl;
    //cout << nn << endl;
    printf("------------------------------------------------------------------------\n");
    return 0;
}
} ;
