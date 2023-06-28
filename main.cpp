#include "Graph.h"
#include "Utility.h"
#include "Timer.h"
#include "popl.hpp"

using namespace std;
using namespace popl;

void print_usage() {
	printf("Example usage: ./kPlexS -g path_to_graph -a exact -k 3 -o\n");
}

int main(int argc, char *argv[]) {
#ifndef NDEBUG
	printf("**** kPlexS (Debug) build at %s %s ***\n", __TIME__, __DATE__);
	printf("!!! You may want to define NDEBUG in Utility.h to get better performance!\n");
#else
	printf("**** kPlexS (Release) build at %s %s ***\n", __TIME__, __DATE__);
#endif

	bool output = false;
	bool binary_input = false;

	OptionParser op("Allowed options");
	auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
	auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
	auto alg_option = op.add<Value<string>>("a", "alg", "\'algorithm name\' (exact | verify)");
	auto k_option = op.add<Value<int>>("k", "k", "\'the value of k for k-plex\'");
	op.add<Switch>("o", "output", "\'write the kplex to ./kplex.txt\'", &output);
	op.add<Switch>("b", "binary", "\'read the input graph from binary files b_adj.bin and b_degree.bin\'", &binary_input);

	op.parse(argc, argv);

	if(help_option->is_set()||argc <= 1) {
		cout << op << endl;
		if(argc <= 1) {
			print_usage();
			return 0;
		}
	}
	if(!graph_option->is_set()) {
		printf("!!! Path to input graph file is not provided! Exit !!!\n");
		return 0;
	}
	if(!k_option->is_set()) {
		printf("!!! k is not provided! Exit !!!\n");
		return 0;
	}

	string alg = "exact";
	if(alg_option->is_set()) alg = alg_option->value();

	Graph *graph = new Graph(graph_option->value().c_str(), k_option->value());
	if(binary_input) graph->read_graph_binary();
	else graph->read_graph();

#ifndef NDEBUG
	printf("\t*** Finished reading graph\n");
#endif

	if(strcmp(alg.c_str(), "exact") == 0) graph->kPlex_exact();
	else if(strcmp(alg.c_str(), "verify") == 0) graph->verify_kplex();
	else print_usage();

	if(output) graph->output_one_kplex();

	delete graph;

	printf("\n");
	return 0;
}
