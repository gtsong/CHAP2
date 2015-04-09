#ifndef GRAPH_BIPAR_H
#define GRAPH_BIPAR_H

struct bipar_node {
	int degree;
	int label; // an event index of list ops
	struct I reg;
	int sp_code;
	struct edge_ent *adj_list;
};

struct edge_ent {
	int endpoint;
	int pid;
	int len;
	int algn_id;
};

struct matrix_ent {
	int weight;
	int tmp_val; // if tmp_val == -1, then either no edge or already picked
	int algn_id;
	bool is_in;
	int val;
	bool mark;
	bool is_required;
	int sign;
};

#endif /* GRAPH_BIPAR_H */
