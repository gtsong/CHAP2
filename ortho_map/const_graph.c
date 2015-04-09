#include "main.h"
#include "const_graph.h"
#include "regions.h"
#include "util_gen.h"
#include "util.h"
#include "util_i.h"
#include "util_ops.h"
#include "apply_ops.h"
#include "read_algn.h"
#include "read_maf.h"

extern int debug_mode;

//void const_graph(int num_algns, struct DotList *algns, FILE *f, struct ops_list *ops, int num_ops, int size) // the initial alignments
void const_graph(int num_algns, struct DotList *algns, FILE *f)
{
	int i = 0, j = 0, k = 0;
	int num_edges = 0, num_nodes = 0;
	struct slist *st;
	int cur_id = 0, tmp_id = 0;
	struct I cur_reg, tmp_reg;
	int lo, hi;
	int num_degree = 1;
	int num_left_nodes = 0;
	int num_right_nodes = 0;
	int entry = 1;
	int cur_entry = 0;
	int cur_weight = 0, old_weight = 0;
	struct bipar_node *graph;
	struct matrix_ent *init_matrix;
	int round = 0;
// an alignment maps to an edge
	for( i = 0; i < num_algns; i++ ) {
		if( (algns[i].x.upper - algns[i].xr_diff - algns[i].x.lower - algns[i].xl_diff) < DEL_TH ) {
			algns[i].sign = DELETED;
		}
		else if( (algns[i].y.upper - algns[i].yr_diff - algns[i].y.lower - algns[i].yl_diff) < DEL_TH ) {
			algns[i].sign = DELETED;
		}
		else if( (algns[i].sign != DELETED) && (algns[i].sp_id == PAIR) ) {
			num_edges++;
		}
	}

	if( num_edges <= 0 ) {}
	else {
		graph = (struct bipar_node *) ckalloc((2 * num_edges) * sizeof(struct bipar_node) );
		initialize_graph(graph, 2*num_edges);
		st = (struct slist *) ckalloc(num_algns * sizeof(struct slist) );

		num_left_nodes = 0;
		for ( round = 0; round <= 1; round++ ) 
		{
//			sort_init_algns(st, algns, num_algns, SELF1);
			initialize_slist(st, 0, num_algns);
			for( i = 0; i < num_algns; i++ ) st[i].id = i;
			sort_init_algns(st, algns, num_algns, round); // SELF1 if round=0, SELF2 if round=1

			k = 1;
			num_degree = 1;
			if( round == 0 ) {
				cur_reg = assign_I(algns[st[0].id].x.lower + algns[st[0].id].xl_diff, algns[st[0].id].x.upper - algns[st[0].id].xr_diff);
			}
			else if( round == 1 ) {
				cur_reg = assign_I(algns[st[0].id].y.lower + algns[st[0].id].yl_diff, algns[st[0].id].y.upper - algns[st[0].id].yr_diff);
			}

			while( (k < num_algns) && (((cur_reg.upper-cur_reg.lower) < DEL_TH) || (algns[st[k-1].id].sign == DELETED) || (algns[st[k-1].id].sp_id != PAIR)) ) {
				if( round == 0 ) {
					cur_reg = assign_I(algns[st[k].id].x.lower + algns[st[k].id].xl_diff, algns[st[k].id].x.upper - algns[st[k].id].xr_diff);
				}
				else if( round == 1 ) {
					cur_reg = assign_I(algns[st[k].id].y.lower + algns[st[k].id].yl_diff, algns[st[k].id].y.upper - algns[st[k].id].yr_diff);
				}
				k++;
			}

			if( k >= num_algns ) { // note that st[k-1].id is for cur_reg
				if(((cur_reg.upper-cur_reg.lower) < DEL_TH) || (algns[st[k-1].id].sign == DELETED) || (algns[st[k-1].id].sp_id != PAIR)) {
					fatal("number of edges should be zero in const_graph.c\n");
				}
			}

			j = num_left_nodes;
			for( i = k; i < num_algns; i++ ) {
				cur_id = st[i].id;

				if( round == 0 ) {
					tmp_reg = assign_I(algns[cur_id].x.lower + algns[cur_id].xl_diff, algns[cur_id].x.upper - algns[cur_id].xr_diff);
				}
				else if( round == 1 ) {
					tmp_reg = assign_I(algns[cur_id].y.lower + algns[cur_id].yl_diff, algns[cur_id].y.upper - algns[cur_id].yr_diff);
				}

				while( ((i+1) < num_algns) && (((tmp_reg.upper-tmp_reg.lower) < DEL_TH) || (algns[cur_id].sign == DELETED) || (algns[cur_id].sp_id != PAIR)) ) {
					cur_id = st[i+1].id;
					if( round == 0 ) {
						tmp_reg = assign_I(algns[cur_id].x.lower + algns[cur_id].xl_diff, algns[cur_id].x.upper - algns[cur_id].xr_diff);
					}
					else if( round == 1 ) {
						tmp_reg = assign_I(algns[cur_id].y.lower + algns[cur_id].yl_diff, algns[cur_id].y.upper - algns[cur_id].yr_diff);
					}
					i++;
				}

				if( (i > (num_algns-1)) || ((tmp_reg.upper-tmp_reg.lower) < DEL_TH) || (algns[st[i].id].sign == DELETED) || (algns[st[i].id].sp_id != PAIR)) {}
				else if( ( algns[cur_id].sign != DELETED ) && ( algns[cur_id].sp_id == PAIR ) ) 
				{
					if( ( strict_almost_equal(tmp_reg, cur_reg) == true ) || ( f_loose_subset(tmp_reg, cur_reg, STRICT) == true ) || ( f_loose_subset(cur_reg, tmp_reg, STRICT) == true)) {
						num_degree++;
						if( cur_reg.lower < tmp_reg.lower ) lo = cur_reg.lower;
						else lo = tmp_reg.lower;
						if( cur_reg.upper > tmp_reg.upper ) hi = cur_reg.upper;
						else hi = tmp_reg.upper;

						if( lo > hi ) {
							fatalf("empty interval: %d-%d", lo, hi);
						}
						cur_reg = assign_I(lo, hi);
					}
					else {
						graph[j].degree = num_degree;
						graph[j].adj_list = (struct edge_ent *) ckalloc(num_degree * (sizeof(struct edge_ent)));
						initialize_edge_ent(graph[j].adj_list, num_degree);
						j++;
						num_degree = 1;
						cur_reg = assign_I(tmp_reg.lower, tmp_reg.upper);
					}
				}
				else {
					fatalf("unexpected case: [%d-%d] and [%d-%d]\n", cur_reg.lower, cur_reg.upper, tmp_reg.lower, tmp_reg.upper);
				}
			}

			graph[j].degree = num_degree;
			graph[j].adj_list = (struct edge_ent *) ckalloc(num_degree * (sizeof(struct edge_ent)));
			initialize_edge_ent(graph[j].adj_list, num_degree);
	
			k = 1;
			cur_id = st[0].id;
			if( round == 0 ) {
				cur_reg = assign_I(algns[st[0].id].x.lower + algns[st[0].id].xl_diff, algns[st[0].id].x.upper - algns[st[0].id].xr_diff);
			}
			else if( round == 1 ) {
				cur_reg = assign_I(algns[st[0].id].y.lower + algns[st[0].id].yl_diff, algns[st[0].id].y.upper - algns[st[0].id].yr_diff);
			}

			while( (k < num_algns) && (((cur_reg.upper-cur_reg.lower) < DEL_TH) || (algns[st[k-1].id].sign == DELETED) || (algns[st[k-1].id].sp_id != PAIR)) ) {
				if( round == 0 ) {
					cur_reg = assign_I(algns[st[k].id].x.lower + algns[st[k].id].xl_diff, algns[st[k].id].x.upper - algns[st[k].id].xr_diff);
				}
				else if( round == 1 ) {
					cur_reg = assign_I(algns[st[k].id].y.lower + algns[st[k].id].yl_diff, algns[st[k].id].y.upper - algns[st[k].id].yr_diff);
				}
				cur_id = st[k].id;
				k++;
			}

			if( k >= num_algns ) { // note that st[k-1].id is for cur_reg
				if(((cur_reg.upper-cur_reg.lower) < DEL_TH) || (algns[st[k-1].id].sign == DELETED) || (algns[st[k-1].id].sp_id != PAIR)) {
					fatal("number of edges should be zero in const_graph.c\n");
				}
			}

			j = num_left_nodes;
			entry = 1;
			graph[j].reg = assign_I(cur_reg.lower, cur_reg.upper);
			graph[j].adj_list[0].pid = algns[cur_id].identity;
			graph[j].adj_list[0].len = width(algns[cur_id].x)-algns[cur_id].xl_diff-algns[cur_id].xr_diff;
			graph[j].adj_list[0].algn_id = cur_id;
			graph[j].adj_list[0].endpoint = -1;
			graph[j].sp_code = round+1; // SP_1=1 if round == 0 and SP_2=2 else

			for( i = k; i < num_algns; i++ ) {
				cur_id = st[i].id;

				if( round == 0 ) {
					tmp_reg = assign_I(algns[cur_id].x.lower + algns[cur_id].xl_diff, algns[cur_id].x.upper - algns[cur_id].xr_diff);
				}
				else if( round == 1 ) {
					tmp_reg = assign_I(algns[cur_id].y.lower + algns[cur_id].yl_diff, algns[cur_id].y.upper - algns[cur_id].yr_diff);
				}

				while( ((i+1) < num_algns) && (((tmp_reg.upper-tmp_reg.lower) < DEL_TH) || (algns[cur_id].sign == DELETED) || (algns[cur_id].sp_id != PAIR)) ) {
					cur_id = st[i+1].id;
					if( round == 0 ) {
						tmp_reg = assign_I(algns[cur_id].x.lower + algns[cur_id].xl_diff, algns[cur_id].x.upper - algns[cur_id].xr_diff);
					}
					else if( round == 1 ) {
						tmp_reg = assign_I(algns[cur_id].y.lower + algns[cur_id].yl_diff, algns[cur_id].y.upper - algns[cur_id].yr_diff);
					}
					i++;
				}

				if( (i > (num_algns-1)) || ((tmp_reg.upper-tmp_reg.lower) < DEL_TH) || (algns[st[i].id].sign == DELETED) || (algns[st[i].id].sp_id != PAIR) ) {}
				else if( ( algns[cur_id].sign != DELETED ) && ( algns[cur_id].sp_id == PAIR ) ) 
				{
					if( ( strict_almost_equal(tmp_reg, cur_reg) == true ) || ( f_loose_subset(tmp_reg, cur_reg, STRICT) == true ) || ( f_loose_subset(cur_reg, tmp_reg, STRICT) == true)) {
						if( cur_reg.lower < tmp_reg.lower ) lo = cur_reg.lower;
						else lo = tmp_reg.lower;
						if( cur_reg.upper > tmp_reg.upper ) hi = cur_reg.upper;
						else hi = tmp_reg.upper;

						if( lo > hi ) {
							fatalf("empty interval: %d-%d", lo, hi);
						}
						cur_reg = assign_I(lo, hi);
						graph[j].reg = assign_I(lo, hi);
						if( entry >= graph[j].degree ) {
							fatalf("edge list for node %d overflow\n", j);
						}
						graph[j].adj_list[entry].pid = algns[cur_id].identity;
						graph[j].adj_list[entry].len = width(algns[cur_id].x) - algns[cur_id].xl_diff-algns[cur_id].xr_diff;
						graph[j].adj_list[entry].algn_id = cur_id;
						graph[j].adj_list[entry].endpoint = -1;
						entry++;
					}
					else {
						j++;	
						graph[j].reg = assign_I(tmp_reg.lower, tmp_reg.upper);
						graph[j].sp_code = round+1; // SP_1=1 if round == 0 and SP_2=2 else
						graph[j].adj_list[0].pid = algns[cur_id].identity;
						graph[j].adj_list[0].len = width(algns[cur_id].x) - algns[cur_id].xl_diff-algns[cur_id].xr_diff;
						graph[j].adj_list[0].algn_id = cur_id;
						graph[j].adj_list[0].endpoint = -1;
						entry = 1;
						cur_reg = assign_I(tmp_reg.lower, tmp_reg.upper);
					}
				}
				else {
	
				}
			}
			j = j+1;
			
			if( round == 0 ) {
				num_left_nodes = j;
			}
			else {
				num_nodes = j;
				num_right_nodes = num_nodes - num_left_nodes;
			}
		}


		init_matrix = (struct matrix_ent *) ckalloc((num_left_nodes * num_right_nodes) * (sizeof(struct matrix_ent)));
		initialize_matrix(init_matrix, num_left_nodes, num_right_nodes);
// map edges in G to R
		for( i = num_left_nodes; i < num_nodes; i++ ) { 
			for( k = 0; k < graph[i].degree; k++ ) {
				cur_id = graph[i].adj_list[k].algn_id;

				if(cur_id == -1) {
					fatalf("unassigned edge %d th in node %d\n", k, i);
				}
				else if( (algns[cur_id].sign == DELETED) || (algns[cur_id].sp_id != PAIR) )
				{
					fatalf("edge unavailable : [%d, %d] and [%d, %d]\n", algns[cur_id].x.lower, algns[cur_id].x.upper, algns[cur_id].y.lower, algns[cur_id].y.upper);
				}
				else {
					cur_reg = assign_I(algns[cur_id].x.lower + algns[cur_id].xl_diff, algns[cur_id].x.upper - algns[cur_id].xr_diff);

					j = 0;
					tmp_reg = assign_I(graph[j].reg.lower, graph[j].reg.upper);
					while( (j < num_left_nodes) && (strict_almost_equal(cur_reg, tmp_reg) == false) && (f_loose_subset(cur_reg, tmp_reg, STRICT) == false) ) {
						j++;
						if( j < num_left_nodes ) {
							tmp_reg = assign_I(graph[j].reg.lower, graph[j].reg.upper);
						}
					}

					if( j >= num_left_nodes ) {
						fatalf("graph construction error in %d-%d\n", cur_reg.lower, cur_reg.upper);
					}
					else {
						tmp_id = j;
						if( graph[i].adj_list[k].endpoint == -1 ) {
							graph[i].adj_list[k].endpoint = tmp_id;
						}
						else {
							fatalf("endpoint for %d is already assigned for %d-%d\n", graph[i].adj_list[k].endpoint, graph[i].reg.lower, graph[i].reg.upper);
						}

						j = 0; 
						while( (tmp_id < num_left_nodes) && (j < graph[tmp_id].degree) && (graph[tmp_id].adj_list[j].algn_id != cur_id) ) {
							j++;
							if( j >= graph[tmp_id].degree ) {
								if( graph[i].adj_list[k].endpoint == tmp_id ) {
									graph[i].adj_list[k].endpoint = tmp_id + 1;
								}
								tmp_id = tmp_id + 1;
								j = 0;
							}
						}

						if( tmp_id >= num_left_nodes ) {
							graph[i].adj_list[k].endpoint = -1;
							fatalf("graph construction error in %d-%d\n", cur_reg.lower, cur_reg.upper);
						}
						else {
							if( (graph[tmp_id].adj_list[j].endpoint == -1) && (graph[tmp_id].adj_list[j].algn_id != -1)) {
								graph[tmp_id].adj_list[j].endpoint = i;	
								cur_entry = (tmp_id * num_right_nodes) + (i - num_left_nodes);
								old_weight = init_matrix[cur_entry].weight;
								cur_weight = (int)(((float)graph[tmp_id].adj_list[j].len) * ((float)graph[tmp_id].adj_list[j].pid/(float)100));
								if( (old_weight == -1) || (cur_weight > old_weight) ) {
									init_matrix[cur_entry].tmp_val = (int)(((float)(graph[tmp_id].adj_list[j].len - countNs(algns, graph[tmp_id].adj_list[j].algn_id, f, SELF1))) * ((float)graph[tmp_id].adj_list[j].pid/(float)100));
									init_matrix[cur_entry].weight = init_matrix[cur_entry].tmp_val;
									init_matrix[cur_entry].algn_id = graph[tmp_id].adj_list[j].algn_id;	
									init_matrix[cur_entry].sign = algns[graph[tmp_id].adj_list[j].algn_id].init_sign;	
								}
							}
							else {
								fatalf("endpoint for %d is already assigned for %d-%d\n", graph[tmp_id].adj_list[j].endpoint, graph[tmp_id].reg.lower, graph[tmp_id].reg.upper);
							}
						}
					}
				}
			}
		}
	
//		mark_event_src_nodes(graph, num_left_nodes, num_right_nodes, ops, num_ops, size);
// 		init_map_in_greedy(init_matrix, num_left_nodes, num_right_nodes, algns, num_algns, f, graph);
 		init_map_in_dp(init_matrix, num_left_nodes, num_right_nodes, algns, f);
//		adjust_event_src(init_matrix, num_left_nodes, num_rigth_nodes, ops, num_ops, size, algns, num_algns, f, graph);

		if( (num_left_nodes > 0) && (num_right_nodes > 0) ) {
//  		max_bipar_weight_match(init_matrix, num_left_nodes, num_right_nodes, algns, f);
//  		max_bipar_weight_match(init_matrix, num_left_nodes, num_right_nodes);
			update_init_algns_for_ancestral(init_matrix, num_left_nodes, num_right_nodes, algns, num_algns);
		}

 		free_graph(graph, 2*num_edges);
		free(graph);
		free(init_matrix);
		free(st);
	}
}

/*
void adjust_event_src(struct matrix_ent *init_matrix, int num_left, int num_right, struct ops_list *ops, int num_ops, int size, struct DotList *algns, int num_algns, FILE *f, struct bipar_node *graph)
{
	int i = 0, j = 0;
	struct ops_list *ops1, *ops2;
	int num_ops1 = 0, num_ops2 = 0;
	struct I reg;
	int sum_width = 0;

	reg = assign_I(0, 1);

	ops1 = (struct ops_list *) ckalloc(num_ops * sizeof(struct ops_list));
	ops2 = (struct ops_list *) ckalloc(num_ops * sizeof(struct ops_list));

	init_ops(ops1, 0, num_ops);
	init_ops(ops2, 0, num_ops);
	num_ops1 = cal_cur_pos_ops(num_ops, ops, ops1, SELF1, 0);
	num_ops2 = cal_cur_pos_ops(num_ops, ops, ops2, SELF2, size);
	for( i = 0; i < num_ops1; i++ ) {
		if( ( ops1[i].sign == '+' ) || ( ops1[i].sign == '-' ) ) {
			reg = assign_I(ops1[i].srcStart, ops1[i].srcEnd);
			if( width(reg) < (3*ERR_TH) ) {}
			else if( is_on_prev_events(reg, ops1, i+1, num_ops1-1) == true ) {}
			else {
				sum_width = 0;
				for( j = 0; j < num_left; j++ ) {
					if( (proper_overlap(reg, graph[j].reg) == true) && (sum_width < width(intersect(reg, graph[j].reg))) ) 
					{
						sum_width = sum_width + intersect(reg, graph[j].reg);
					}
				}
			}
		}
	}

	for( i = 0; i < num_ops2; i++ ) {
		if( ( ops2[i].sign == '+' ) || ( ops2[i].sign == '-' ) ) {
			reg = assign_I(ops2[i].srcStart, ops2[i].srcEnd);
			if( width(reg) < (3*ERR_TH) ) {}
			else if( is_on_prev_events(reg, ops2, i+1, num_ops2-1) == true ) {}
			else {
				ops2[i].id;
			}
		}
	}
	free(ops1);
	free(ops2);
}
*/

void mark_event_src_nodes(struct bipar_node *graph, int num_left, int num_right, struct ops_list *ops, int num_ops, int size)
{
	int i = 0;
	struct ops_list *ops1, *ops2;
	int num_ops1 = 0, num_ops2 = 0;
	struct I reg;

	reg = assign_I(0, 1);

	ops1 = (struct ops_list *) ckalloc(num_ops * sizeof(struct ops_list));
	ops2 = (struct ops_list *) ckalloc(num_ops * sizeof(struct ops_list));

	init_ops(ops1, 0, num_ops);
	init_ops(ops2, 0, num_ops);
	num_ops1 = cal_cur_pos_ops(num_ops, ops, ops1, SELF1, 0);
	num_ops2 = cal_cur_pos_ops(num_ops, ops, ops2, SELF2, size);
	for( i = 0; i < num_ops1; i++ ) {
		if( ( ops1[i].sign == '+' ) || ( ops1[i].sign == '-' ) ) {
			reg = assign_I(ops1[i].srcStart, ops1[i].srcEnd);
			if( width(reg) < (3*ERR_TH) ) {}
			else if( is_on_prev_events(reg, ops1, i+1, num_ops1-1) == true ) {}
			else {
				mark_label_on_graph(reg, graph, 0, num_left-1, ops1[i].id);
			}
		}
	}

	for( i = 0; i < num_ops2; i++ ) {
		if( ( ops2[i].sign == '+' ) || ( ops2[i].sign == '-' ) ) {
			reg = assign_I(ops2[i].srcStart, ops2[i].srcEnd);
			if( width(reg) < (3*ERR_TH) ) {}
			else if( is_on_prev_events(reg, ops2, i+1, num_ops2-1) == true ) {}
			else {
				mark_label_on_graph(reg, graph, num_left, num_left+num_right-1, ops2[i].id);
			}
		}
	}
	free(ops1);
	free(ops2);
}

void mark_label_on_graph(struct I reg, struct bipar_node *graph, int from, int to, int ops_id)
{
	int i = 0;
	int max_id = -1;
	int max_width = 0;

	for( i = from; i <= to; i++ ) {
		if( (proper_overlap(reg, graph[i].reg) == true) && (max_width < width(intersect(reg, graph[i].reg))) ) 
		{
			max_id = i;
			max_width = width(intersect(reg, graph[i].reg));
		}
	}
	
	if( max_id != -1 ) {
		graph[max_id].label = ops_id;
	}
}

void init_map_in_greedy(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, int num_init_algns, FILE *f, struct bipar_node *graph)
{
	int i = 0;
	int *id1, *id2;
	int cur_max_val = -1;
	int num_entry = 0;
	int entry = 0;
	bool *src_node_l, *src_node_r;
	int prev_id1 = -1, prev_id2 = -1;

	src_node_l = (bool *) ckalloc(num_left * sizeof(bool));
	src_node_r = (bool *) ckalloc(num_right * sizeof(bool));
	id1 = (int *) ckalloc(sizeof(int));
	id2 = (int *) ckalloc(sizeof(int));
	if( num_left > num_right ) num_entry = num_left;
	else num_entry = num_right;

	for( i = 0; i < num_left; i++ ) {
		src_node_l[i] = false;
		if( graph[i].label != -1 ) src_node_l[i] = true;
	}
	for( i = 0; i < num_right; i++ ) {
		src_node_r[i] = false;
		if( graph[num_left+i].label != -1 ) src_node_r[i] = true;
	}

	i = 1;
	*id1 = -1;
	*id2 = -1;
	cur_max_val = pick_max_element(init_matrix, num_left, num_right, id1, id2, src_node_l, src_node_r, prev_id1, prev_id2);
	while( (cur_max_val != -1) && (i < num_entry) ) {
		if( ((*id1) != -1) && ((*id2) != -1) ) {
			if( cur_max_val != -1 ) {
				entry = (*id1) * num_right + (*id2);
				init_matrix[entry].is_in = true;
				src_node_l[*id1] = false;
				src_node_r[*id2] = false;
				graph[*id1].label = -1;
				graph[num_left+(*id2)].label = -1;
			}
		}

		if( init_matrix[entry].algn_id == -1 ) {
			fatalf("alignment not assigned to edge (%d,%d) in the matrix\n", *id1, *id2);
		}
		else {
			update_matrix_tmp_val(graph, init_matrix, num_left, num_right, *id1, *id2, init_algns, num_init_algns, f);
		}
		prev_id1 = *id1;
		prev_id2 = *id2;
		cur_max_val = pick_max_element(init_matrix, num_left, num_right, id1, id2, src_node_l, src_node_r, prev_id1, prev_id2);
		i++;
	}

	free(id1);
	free(id2);
	free(src_node_r);
	free(src_node_l);
}

void update_matrix_tmp_val(struct bipar_node *graph, struct matrix_ent *init_matrix, int num_left, int num_right, int l_entry, int r_entry, struct DotList *init_algns, int num_init_algns, FILE *f)
{
	int b = 0, e = 0;	
	int i = 0, j = 0, k = 0;
	int cur_id = 0, init_cur_id = 0;
	struct I cur_reg, tmp_reg1, tmp_reg2;
	char S1[BIG], T1[BIG];
	int beg1 = 0, end1 = 0, len1 = 0;
	struct b_list *a_info;
	int cur_len = 0, end_pos = 0;
	int *t_b;
	float pid;
	int cur_entry = 0;
	int b1 = 0, e1 = 0, b2 = 0, e2 = 0;
	int entry = 0;
	
	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	t_b = (int *) ckalloc(sizeof(int));
	*t_b = 0;
	cur_reg = assign_I(0,1);
	tmp_reg1 = assign_I(0,1);
	tmp_reg2 = assign_I(0,1);

	if( (l_entry >= num_left) || (r_entry >= num_right) ) {
		fatalf("out of init_matrix: %d,%d\n", l_entry, r_entry);
	}
	entry = l_entry * num_right + r_entry;
	init_cur_id = init_matrix[entry].algn_id;
	init_matrix[entry].tmp_val = -1;

	if( init_cur_id == -1 ) {
		fatalf("unassigned entry in matrix [%d,%d]\n", l_entry, r_entry);
	}
	else if( init_cur_id >= num_init_algns ) {
		fatalf("out of init_algns: %d\n", cur_id);
	}
	else {
		for( k = 0; k < 2; k++ ) {
			if( k == 0 ) {
				cur_reg = assign_I(init_algns[init_cur_id].x.lower+init_algns[init_cur_id].xl_diff, init_algns[init_cur_id].x.upper-init_algns[init_cur_id].xr_diff);
				cur_entry = l_entry;
			}
			else {
				cur_reg = assign_I(init_algns[init_cur_id].y.lower+init_algns[init_cur_id].yl_diff, init_algns[init_cur_id].y.upper-init_algns[init_cur_id].yr_diff);
				cur_entry = r_entry + num_left;
			}
//			graph[cur_entry].reg = assign_I(0,1);
			i = cur_entry - 1;
			while((((k == 0) && (i >= 0)) || ((k == 1) && (i >= num_left))) && (proper_overlap(graph[i].reg, cur_reg) == true)) {
				i--;
			}
			b = i+1;
			i = cur_entry + 1;
			while((((k == 0) && (i < num_left)) || ((k == 1) && (i < (num_left+num_right)))) && (proper_overlap(graph[i].reg, cur_reg) == true)) {
				i++;
			}
			e = i-1;

			if( k == 0 ) {
				b1 = b;
				e1 = e;
				b2 = 0;
				e2 = num_right-1;
			}
			else {
				b1 = 0;
				e1 = num_left-1;
				b2 = b - num_left;
				e2 = e - num_left;
			}
	
			for( i = b1; i <= e1; i++ ) { // in species1
				for( j = b2; j <= e2; j++ ) { // in species2
					entry = i * num_right + j;
					cur_id = init_matrix[entry].algn_id;
					if( (i == l_entry) && (j == r_entry) ) {} // it is saved
/*					else if( (i == l_entry) || (j == r_entry) ) {
						init_matrix[entry].tmp_val = -1;
					} */
					else if( (init_matrix[entry].tmp_val != -1) && (cur_id != -1) ) {
						if( k == 0 ) {
							tmp_reg1 = assign_I(init_algns[cur_id].x.lower+init_algns[cur_id].xl_diff, init_algns[cur_id].x.upper-init_algns[cur_id].xr_diff);
						}
						else {
							tmp_reg1 = assign_I(init_algns[cur_id].y.lower+init_algns[cur_id].yl_diff, init_algns[cur_id].y.upper-init_algns[cur_id].yr_diff);
						}
	
						if( (strict_almost_equal(tmp_reg1, cur_reg) == true) || (f_loose_subset(tmp_reg1, cur_reg, STRICT) == true) || (f_loose_subset(cur_reg, tmp_reg1, STRICT) == true)) {
							init_matrix[entry].tmp_val = -1;	
						}
						else if( proper_overlap(tmp_reg1, cur_reg) == true ) {
							beg1 = init_algns[cur_id].xl_diff + init_algns[cur_id].xl_offset;	
							len1 = init_algns[cur_id].x.upper - init_algns[cur_id].xr_diff - init_algns[cur_id].xr_offset - init_algns[cur_id].x.lower - init_algns[cur_id].xl_offset;
							end1 = find_xloc_one(init_algns[cur_id], f, len1, NO_GAP_INC);
							get_nth_algn(S1, T1, init_algns[cur_id].fid, beg1, f, a_info, REG);

							if( S1[strlen(S1)-1] == '\n' ) end_pos = strlen(S1)-1;
							else end_pos = strlen(S1);

							if( end1 < end_pos ) {
								S1[end1] = '\0';
								end_pos = end1-1;
							}

							if( tmp_reg1.lower < cur_reg.lower ) {
								tmp_reg2 = assign_I(tmp_reg1.lower, cur_reg.lower);
							}
							else {
								tmp_reg2 = assign_I(cur_reg.upper, tmp_reg1.upper);
							}

							if( init_algns[cur_id].init_sign == 0 ) {
								cur_len = count_ncol(tmp_reg1, tmp_reg2, S1, end_pos, t_b);
							}
							else if( (k == 1) && (init_algns[cur_id].init_sign == 1) ) {
								cur_len = count_ncol_rev(tmp_reg1, tmp_reg2, T1, end_pos, t_b);
							}
							else if( (k == 0) && (init_algns[cur_id].init_sign == 1) ) {
								cur_len = count_ncol(tmp_reg1, tmp_reg2, T1, end_pos, t_b);
							}

							pid = cal_pid_maf_beg(S1, T1, *t_b, cur_len);
							init_matrix[entry].tmp_val = (int)((float)width(tmp_reg2))*(pid/(float)100);
						}
					}
				}
			}
		}
	}

	free(t_b);
	free(a_info);
}

void initialize_graph(struct bipar_node *graph, int num_nodes)
{
	int i = 0;

	for( i = 0; i < num_nodes; i++ ) {
		graph[i].degree = 0;
		graph[i].label = -1;
  	graph[i].reg = assign_I(0,1);
		graph[i].sp_code = -1;
		graph[i].adj_list = NULL;
	}
}

void initialize_edge_ent(struct edge_ent *edge, int num_edges)
{
	int i = 0;

	for( i = 0; i < num_edges; i++ ) {
  	edge[i].endpoint = -1;
  	edge[i].pid = 0;
  	edge[i].len = 0;
  	edge[i].algn_id = -1;
	}
}

void initialize_matrix(struct matrix_ent * init_matrix, int num_left, int num_right)
{
	int i = 0, j = 0;
	int cur_entry = 0;

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			cur_entry = i * num_right + j;
			init_matrix[cur_entry].weight = -1;
			init_matrix[cur_entry].tmp_val = -1;
			init_matrix[cur_entry].algn_id = -1;
			init_matrix[cur_entry].is_in = false;
			init_matrix[cur_entry].mark = false;
			init_matrix[cur_entry].val = -1;
			init_matrix[cur_entry].sign = DELETED;
			init_matrix[cur_entry].is_required = false;
		}
	}
}

void pick_max_unvisited(int *l_dist, int *r_dist, int num_left, int num_right, int *id1, int *id2, bool *l_visited, bool *r_visited)
{
	int j = 0;
	int val = 0;

	for( j = 0; j < num_right; j++ ) {
		if( (val < r_dist[j]) && (r_visited[j] == false)) 
		{
			val = r_dist[j];
			*id1 = -1;
			*id2 = j;
		}
	}

	for( j = 0; j < num_left; j++ ) {
		if( (val < l_dist[j]) && (l_visited[j] == false) ) {
			val = l_dist[j];
			*id1 = j;
			*id2 = -1;
		}
	}

	if( val == 0 ) {
		*id1 = -1;
		*id2 = -1;
	}
}

int pick_max_element(struct matrix_ent *init_matrix, int num_left, int num_right, int *id1, int *id2, bool *src_l, bool *src_r, int prev_id1, int prev_id2)
{
	int i = 0, j = 0;
	int val = -1;
	int tmp_id1 = -1, tmp_id2 = -1;
	int entry = 0;
	bool is_all_false = true;
	int count = 0;
	int only_id = -1;

	if( (prev_id1 == -1) && (prev_id2 == -1) ) {
		while( i < num_left ) {
			count = 0;
			for( j = 0; j < num_right; j++ ) {
				entry = i * num_right + j;
				if( init_matrix[entry].tmp_val != -1 ) {
					count++;
					if( count == 1 ) only_id = j;
				}
			}

			if( (count == 1) && (only_id != -1) ) {
				if( (val == -1) || (val < init_matrix[entry].tmp_val) ) {
					tmp_id1 = i;
					tmp_id2 = only_id;
					entry = tmp_id1 * num_right + tmp_id2;
					val = init_matrix[entry].tmp_val;
				}
			}

			only_id = -1;
			i++;
		}
	}

	for( i = 0; i < num_left; i++ ) {
		count = 0;
		for( j = 0; j < num_right; j++ ) {
			entry = i * num_right + j;
			if( init_matrix[entry].tmp_val != -1 ) {
				count++;
			}
		}

		if( count == 0 ) {
			src_l[i] = false;
		}
	}

	for( i = 0; i < num_right; i++ ) {
		count = 0;
		for( j = 0; j < num_left; j++ ) {
			entry = j * num_right + i;
			if( init_matrix[entry].tmp_val != -1 ) {
				count++;
			}
		}

		if( count == 0 ) {
			src_r[i] = false;
		}
	}

	for( i = 0; i < num_left; i++ ) {
		if( src_l[i] == true ) is_all_false = false;
		i++;
	}

	for( i = 0; i < num_right; i++ ) {
		if( src_r[i] == true ) is_all_false = false;
		i++;
	}

	if( val == -1 ) {
		for( i = 0; i < num_left; i++ ) {
			for( j = 0; j < num_right; j++ ) {
				if( (is_all_false == true) || (src_l[i] == true) || (src_r[j] == true) ) 
				{
					entry = i * num_right + j;
					if( (init_matrix[entry].tmp_val != -1) && (init_matrix[entry].is_in == false) && (val < init_matrix[entry].tmp_val) ) {
						val = init_matrix[entry].tmp_val;
						tmp_id1 = i;
						tmp_id2 = j;
					}
				}
			}
		}
	}

	if( val == -1 ) {
		tmp_id1 = -1;
		tmp_id2 = -1;
	}

	*id1 = tmp_id1;
	*id2 = tmp_id2; 

	return(val);
}

//void max_bipar_weight_match(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f)
void max_bipar_weight_match(struct matrix_ent *init_matrix, int num_left, int num_right)
{
	bool *l_matched, *r_matched;
	int *r_dist, inc_weight;
	int i = 0, j = 0, k = 0;
	int max_lid = -1, max_rid = -1, max_inc = 0;
	int entry = 0;
	bool is_all_mapped = true;
	int *l_prev, *r_prev, *max_l_prev, *max_r_prev;
	int dest = -1;
	int *l_weight, *r_weight;

	r_dist = (int *) ckalloc(sizeof(int) * num_right);

	l_matched = (bool *) ckalloc(num_left * sizeof(bool));
	r_matched = (bool *) ckalloc(num_right * sizeof(bool));
	l_prev = (int *) ckalloc(num_left * sizeof(int));
	r_prev = (int *) ckalloc(num_right * sizeof(int));
	max_l_prev = (int *) ckalloc(num_left * sizeof(int));
	max_r_prev = (int *) ckalloc(num_right * sizeof(int));
	l_weight = (int *) ckalloc(num_left * sizeof(int));
	r_weight = (int *) ckalloc(num_right * sizeof(int));

	for( i = 0; i < num_left; i++ ) {
		l_matched[i] = false;
		l_weight[i] = 0;
	}
	for( i = 0; i < num_right; i++ ) {
		r_matched[i] = false;
		r_weight[i] = 0;
	}

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			entry = i * num_right + j;
			if( init_matrix[entry].is_in == true ) {
				l_weight[i] = l_weight[i] + init_matrix[entry].weight;
				r_weight[j] = r_weight[j] + init_matrix[entry].weight;
			}
			else {
			}
		}
	}

	mark_matched(l_matched, num_left, num_right, init_matrix, LEFT);
	mark_matched(r_matched, num_right, num_left, init_matrix, RIGHT);

	for( i = 0; i < num_left; i++ ) {
		if( l_matched[i] == false ) is_all_mapped = false;
	}

	if( is_all_mapped == false ) {
		for( i = 0; i < num_right; i++ ) {
			if( r_matched[i] == false ) is_all_mapped = false;
		}
	}
	
	if( is_all_mapped == true ) max_inc = -1;
	else max_inc = 0;

	k = 0;
	while( (k < (num_left+num_right)) && (max_inc != -1) ) {
		max_lid = -1;
		max_rid = -1;
		max_inc = -1;
		for( i = 0; i < num_left; i++ ) {
			if( l_matched[i] == false ) {
				dest = -1;
				dest = shortest_augmenting_path(i, init_matrix, num_left, num_right, l_prev, r_prev, r_dist, r_matched, r_weight); // this search always begins from L

				if( (dest != -1) && (dest >= 0) && (dest < num_right) && (r_prev[dest] != -1) ) { 
					if( r_matched[dest] == false ) {
						if( r_dist[dest] > 0 ) {
							inc_weight = cal_inc_weight_for_path(init_matrix, num_left, num_right, l_prev, r_prev, dest, i, FIRST_RUN); // calculating weight
							if( inc_weight <= 0 ) {}
							else if( max_inc < inc_weight ) {
								max_inc = inc_weight;
								max_lid = i;
								max_rid = dest;
								for( j = 0; j < num_left; j++ ) {
									max_l_prev[j] = l_prev[j];
								}
								for( j = 0; j < num_right; j++ ) {
									max_r_prev[j] = r_prev[j];
								}
							}
						}
					}
				}
			}
		}

		if( (max_lid == -1) || (max_rid == -1) ) {
			max_inc = -1;
		}	
		else {
			max_inc = cal_inc_weight_for_path(init_matrix, num_left, num_right, max_l_prev, max_r_prev, max_rid, max_lid, SECOND_RUN); // updating init_matrix and init_algns
			mark_matched(l_matched, num_left, num_right, init_matrix, LEFT);
			mark_matched(r_matched, num_right, num_left, init_matrix, RIGHT);
		}

		k++;
	}

	free(max_l_prev);
	free(max_r_prev);
	free(l_prev);
	free(r_prev);
	free(r_dist);
	free(l_matched);
	free(r_matched);
}

int cal_inc_weight_for_path(struct matrix_ent *init_matrix, int num_left, int num_right, int *l_prev, int *r_prev, int id, int source, int run_mode)
{
	int i = 0, j = 0;
	int new_score = 0;
	int mode = RIGHT;
	int cur_prev = 0;
	int entry = -1;

	i = id;
	cur_prev = r_prev[id];
	mode = RIGHT;
	while((j < (num_left + num_right)) && (cur_prev != -1) )	
	{
		if( mode == RIGHT ) {
			entry = cur_prev * num_right + i;
			if( init_matrix[entry].weight == -1 ) {
				fatalf("unexpected unaligned nodes: %d vs %d\n", cur_prev, i);
			}
			else {
				new_score = new_score + init_matrix[entry].weight;
			}

			i = cur_prev;
			cur_prev = l_prev[i];
			mode = LEFT;
		}
		else if( mode == LEFT ) {
			entry = i * num_right + cur_prev;
			if( init_matrix[entry].weight == -1 ) {
				fatalf("unexpected unaligned nodes: %d vs %d\n", cur_prev, i);
			}
			else {
				new_score = new_score - init_matrix[entry].weight;
			}

			i = cur_prev;
			cur_prev = r_prev[i];
			if( cur_prev == -1 ) {
				if( i == source ) {}
				else {
					fatalf("path from %d to %d has a problem\n", id, source);
				}
			}
			mode = RIGHT;
		}

		if( run_mode == SECOND_RUN ) {
			if( (cur_prev == -1) || (i == -1) ) {}
			else if( init_matrix[entry].is_in == false ) {
				init_matrix[entry].is_in = true;
			}
			else {
				init_matrix[entry].is_in = false;
			}
		}
		j++;
	}	

	return(new_score);
}

int cal_inc_weight(struct matrix_ent *init_matrix, int num_left, int num_right, int *l_prev, int *r_prev, int id, int source, struct DotList *init_algns, FILE *f)
{
	int i = 0, j = 0;
	int cur_prev = 0;
	int mode = RIGHT;
	int res = 0;
	int num_match = 0, num_unmatch = 0;
	int cur_id = 0;
	int cur_score = 0, new_score = 0;
	struct DotList *ortho_candi;
	int num_ortho = 0;
	int entry = 0;

	i = id;
	cur_prev = r_prev[i];
	mode = LEFT;
	while((j < (num_left + num_right)) && (cur_prev != -1) && (!((mode == LEFT) &&(cur_prev == source))))	
	{
		if( mode == LEFT ) {
			entry = cur_prev * num_right + i;
			if( init_matrix[entry].tmp_val == -1 ) num_unmatch++;
			else num_match++;

			i = cur_prev;
			cur_prev = l_prev[i]; // the node in L connected to i in the augmenting path
			mode = RIGHT;
		}
		else if( mode == RIGHT ) {
			entry = i * num_right + cur_prev;
			if( init_matrix[entry].tmp_val == -1 ) num_unmatch++;
			else num_match++;

			i = cur_prev;
			cur_prev = r_prev[i];
			mode = LEFT;
		}
		j++;
	}

	if( (mode == LEFT) && (cur_prev == source) ) {
		cur_score = 0;
	}
	else {
//		fatalf("inappropriate path: %d-%d\n", cur_prev, mode);
		if( debug_mode == true ) {
			printf("Warning: inappropriate path: %d-%d\n", cur_prev, mode);
		}
		cur_score = -1;
	}

	if( cur_score == -1 ) {
		res = 0;
	}
	else {  
		for( i = 0; i < num_left; i++ ) {
			for( j = 0; j < num_right; j++ ) {
				entry = i * num_right + j;
				if( init_matrix[entry].tmp_val > 0 ) {
					cur_score = cur_score + init_matrix[entry].tmp_val;
					init_matrix[entry].mark = true;
					init_matrix[entry].val = init_matrix[entry].tmp_val;
				}
			}
		}

		num_match = 0;
		num_unmatch = 0;
		i = id;
		cur_prev = r_prev[i];
		mode = LEFT;
		new_score = cur_score;
		while((j < (num_left + num_right)) && (cur_prev != -1) )	
		{
			if( mode == LEFT ) {
				entry = cur_prev * num_right + i;
				if( init_matrix[entry].tmp_val == -1 ) {
					init_matrix[entry].mark = true;
					num_unmatch++;
				}
				else {
					init_matrix[entry].mark = false;
				init_matrix[entry].val = -1;
				new_score = new_score - init_matrix[entry].tmp_val;
				num_match++;
				}

				i = cur_prev;
				cur_prev = l_prev[i];
				mode = RIGHT;
			}
			else if( mode == RIGHT ) {
				entry = i * num_right + cur_prev;
				if( init_matrix[entry].tmp_val == -1 ) {
					init_matrix[entry].mark = true;
					num_unmatch++;
				}
				else {
					init_matrix[entry].mark = false;
					init_matrix[entry].val = -1;
					new_score = new_score - init_matrix[entry].tmp_val;
					num_match++;
				}

				i = cur_prev;
				cur_prev = r_prev[i];
				mode = LEFT;
			}
			j++;
		}

		for( i = 0; i < num_left; i++ ) {
			for( j = 0; j < num_right; j++ ) {
				entry = i * num_right + j;
				if( init_matrix[entry].mark == true ) num_ortho++;
			}
		}

		ortho_candi = (struct DotList *) ckalloc(sizeof(struct DotList) * (num_ortho));
		num_ortho = 0;
		for( i = 0; i < num_left; i++ ) {
			for( j = 0; j < num_right; j++ ) {
				entry = i * num_right + j;
				if( (init_matrix[entry].mark == true) && (init_matrix[entry].val > 0 ) ) {
					cur_id = init_matrix[entry].algn_id;
					assign_algn(ortho_candi, num_ortho, init_algns[cur_id]); 
					num_ortho++;
				}
			}
		}
	
		for( i = 0; i < num_left; i++ ) {
			for( j = 0; j < num_right; j++ ) {
				entry = i * num_right + j;
				if( (init_matrix[entry].mark == true) && (init_matrix[entry].val == -1) ) {
					update_init_matrix_val(init_matrix, i, j, num_right, init_algns, ortho_candi, num_ortho, f);
					if( init_matrix[entry].val != -1 ) {
						new_score = new_score + init_matrix[entry].val;
						assign_algn(ortho_candi, num_ortho, init_algns[init_matrix[entry].algn_id]);
						num_ortho++;
					}
				}
			}
		}		
		
		if( (new_score - cur_score) > 0 ) {
			res = new_score - cur_score;
		}
		else {
			res = 0;
		}
	
		free(ortho_candi);
	}

	return(res);
}

void mark_matched(bool *matched, int num1, int num2, struct matrix_ent *init_matrix, int flag)
{
	int i = 0, j = 0;
	bool cur_status = false;
	int entry = 0;

	for( i = 0; i < num1; i++ ) {
		j = 0;
		cur_status = false;
		while( (j < num2) && (cur_status == false) ) {
			if( flag == LEFT ) {
				entry = i * num2 + j;
				if( init_matrix[entry].is_in == true ) cur_status = true;
			}
			else {
				entry = j * num1 + i;
				if( init_matrix[entry].is_in == true ) cur_status = true;
			}
			j++;
		}
		if( cur_status == true ) {
			matched[i] = true;
		}
		else {
			matched[i] = false;
		}
	}
}


int shortest_augmenting_path(int v, struct matrix_ent *init_matrix, int num_left, int num_right, int *l_prev, int *r_prev, int *r_dist, bool *r_matched, int *r_weight)  
{
	int i = 0, j = 0,  k = 0; 
	bool is_failed = false, all_r_visited = false;
	bool *l_visited, *r_visited;
	int *l_dist;
	int mode = LEFT;
	int alt = -1;
	int entry = 0;
	int dest = -1;
	int count = 0;
	bool is_reached = false;

	l_visited = (bool *) ckalloc(sizeof(bool) * num_left);
	r_visited = (bool *) ckalloc(sizeof(bool) * num_right);
	l_dist = (int *) ckalloc(sizeof(int) * num_left);

	for( i = 0; i < num_right; i++ ) {
		r_visited[i] = false;
		r_prev[i] = -1;
		r_dist[i] = 0;
	}

	for( i = 0; i < num_left; i++ ) {
		l_visited[i] = false;
		l_prev[i] = -1;
		l_dist[i] = 0; // sum of weights of a maximum weighted augmenting path from a source to i in L
	}

	k = 0;
	while((is_reached == false) && (j != -1) && (k < (num_right+num_left)) && (is_failed == false) && (all_r_visited == false))
	{
		if( k == 0 ) {
			j = v;			
			mode = LEFT;
			l_visited[j] = true;
		}
/*
		else {
			pick_max_unvisited(l_dist, r_dist, num_left, num_right, id1, id2, l_visited, r_visited);

			if( ((*id1) == -1) && ((*id2) == -1) ) {
				j = -1;
			}
			else if( (*id1) == -1 ) {
				j = *id1;
				l_visited[j] = true; // remove j from Q
				mode = LEFT;
				if( l_dist[j] == 0 ) j = -1; 
			}
			else if( (*id2) == -1 ) {
				j = *id2;
				r_visited[j] = true; // remove j from Q
				mode = RIGHT;
				if( r_dist[j] == 0 ) j = -1; 
			}
			else {
				fatalf("two nodes picked : %d and %d\n", *id1, *id2);
			}
		}
*/

		if( j != -1 ) { 
			if( mode == LEFT ) {
				count = 0;
				dest = -1;
				for( i = 0; i < num_right; i++ ) {
					entry = j * num_right + i;
					if( ((k > 0) || ((k == 0) && (r_matched[i] == true))) && (init_matrix[entry].is_in == false) && (init_matrix[entry].weight != -1) && (r_visited[i] != true) ) {
						alt = init_matrix[entry].weight - r_weight[i];	// recalculate "alt"

						if( count == 0 ) {
							r_dist[i] = alt;
							r_prev[i] = j;
							dest = i;
						}
						else if( alt > r_dist[i] ) {
							r_dist[i] = alt;
							r_prev[i] = j;
							dest = i;
						}

						count++;
					}
				}

				if( (k > 0) && (dest != -1) && (r_matched[dest] == false) ) {
					is_reached = true;
				}

				if( (dest == -1) || (count == 0) ) {
					is_failed = true;
				}
				else {
					r_visited[dest] = true;
					mode = RIGHT;
				}
			}
			else if( mode == RIGHT ) {
				count = 0;
				dest = -1;
				for( i = 0; i < num_left; i++ ) {
					entry = i * num_right + j;
					if( (init_matrix[entry].is_in == true) && (init_matrix[entry].weight != -1) && (l_visited[i] != true) ) {
						alt = init_matrix[entry].weight - r_weight[j];	

						if( count == 0 ) {
							l_dist[i] = alt;
							l_prev[i] = j;
							dest = i;
						}
						else if( alt > l_dist[i] ) {
							l_dist[i] = alt;
							l_prev[i] = j;
							dest = i;
						}
						count++;
					}
				}

				if( (dest == -1) || (count == 0) ) {
					is_failed = true;
				}
				else {
					l_visited[dest] = true;
					mode = LEFT;
				}
			}

			j = dest;
		}
	
		all_r_visited = true;
		for( i = 0; i < num_right; i++ ) {
			if( r_visited[i] ==  false ) all_r_visited = false;
		}

		if( all_r_visited == false ) {
			all_r_visited = true;
			for( i = 0; i < num_left; i++ ) {
				if( l_visited[i] ==  false ) all_r_visited = false;
			}
		}

		k++;
	}

	free(l_visited);
	free(r_visited);
	free(l_dist);

	return(dest);
}

int update_init_matrix_val(struct matrix_ent *init_matrix, int l_entry, int r_entry, int num_right, struct DotList *init_algns, struct DotList *candi_algns, int num_algns, FILE *f)
{
	int b = 0, e = 0;	
	int i = 0, k = 0;
	int cur_id = 0;
	struct I cur_reg, tmp_reg;
	char S1[BIG], T1[BIG];
	struct b_list *a_info;
	float pid;
	int b1 = 0, e1 = 0, b2 = 0, e2 = 0;
	struct slist *st;
	int res = 0;
	int len = 0;
	int entry = 0;
	
	st = (struct slist *) ckalloc(sizeof(struct slist) * num_algns);
	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));
	cur_reg = assign_I(0,1);
	tmp_reg = assign_I(0,1);

	entry = l_entry * num_right + r_entry;
	cur_id = init_matrix[entry].algn_id;
	for( k = 0; k < 2; k++ ) {
		initialize_slist(st, 0, num_algns);
		for( i = 0; i < num_algns; i++ ) st[i].id = i;
		if( k == 0 ) {
			cur_reg = assign_I(init_algns[cur_id].x.lower, init_algns[cur_id].x.upper);
			sort_init_algns(st, candi_algns, num_algns, SELF1);
		}
		else {
			cur_reg = assign_I(init_algns[cur_id].y.lower, init_algns[cur_id].y.upper);
			sort_init_algns(st, candi_algns, num_algns, SELF2);
		}

		i = 0;
		if( k == 0 ) {
			tmp_reg = assign_I(candi_algns[st[i].id].x.lower, candi_algns[st[i].id].x.upper);
		}
		else {
			tmp_reg = assign_I(candi_algns[st[i].id].y.lower, candi_algns[st[i].id].y.upper);
		}

		while( (i < num_algns) && (cur_reg.upper > tmp_reg.lower) && ((b1+e1) < width(cur_reg)) && ((b2+e2) < width(cur_reg)) ) {
			if( (f_loose_subset(cur_reg, tmp_reg, STRICT) == true) || (f_loose_subset(tmp_reg, cur_reg, STRICT) == true) ) {
				if( k == 0 ) {
					b1 = width(cur_reg);
					e1 = width(cur_reg);
				}
				else {
					b2 = width(cur_reg);
					e2 = width(cur_reg);
				}
			}
			else if( proper_overlap(cur_reg, tmp_reg) == true ) {
				if( (cur_reg.lower >= tmp_reg.lower) && (cur_reg.lower <= tmp_reg.upper) ) 
				{
					if( (k == 0) && (tmp_reg.upper - cur_reg.lower) > b1  ) {
						b1 = tmp_reg.upper - cur_reg.lower;
					}
					else if( (k == 1) && (tmp_reg.upper - cur_reg.lower) > b2  ) {
						b2 = tmp_reg.upper - cur_reg.lower;
					}
				}
				else {
					if( (k == 0) && ((tmp_reg.lower - cur_reg.upper) > e1)  ) {
						e1 = tmp_reg.lower - cur_reg.upper;
					}
					else if( (k == 1) && ((tmp_reg.lower - cur_reg.upper) > e2)  ) {
						e2 = tmp_reg.lower - cur_reg.upper;
					}
				}
			}

			i++;
			if( i >= num_algns ) {}
			else {
				if( k == 0 ) {
					tmp_reg = assign_I(candi_algns[st[i].id].x.lower, candi_algns[st[i].id].x.upper);
				}
				else {
					tmp_reg = assign_I(candi_algns[st[i].id].y.lower, candi_algns[st[i].id].y.upper);
				}
			}
		}
	}

	len = width(cur_reg) - b1 - e1;
	if( len > ( width(cur_reg) - b2 - e2 ) ) {
		len = width(cur_reg) - b2 - e2;
	}

	if( len <= 0 ) {
		res = 0;
	}
	else if( ((b1 + e1) >= (width(init_algns[cur_id].x) - DEL_TH)) || ((b2 + e2) >= (width(init_algns[cur_id].y) - DEL_TH))) {
		res = 0;
	}
	else {
		if( init_algns[cur_id].sign == 0 ) {
			b = find_yloc_one(init_algns[cur_id], f, b2, GAP_INC_IN_Y);			
			b = find_xloc_one(init_algns[cur_id], f, b, GAP_INC) - init_algns[cur_id].x.lower;
		}
		else if( init_algns[cur_id].sign == 1 ) {
			b = find_yloc_one(init_algns[cur_id], f, e2, GAP_INC_IN_Y);			
			b = find_xloc_one(init_algns[cur_id], f, b, GAP_INC) - init_algns[cur_id].x.lower;
		}
		if( b < b1 ) b = b1;

		e1 = find_xloc_one(init_algns[cur_id], f, width(init_algns[cur_id].x)-e1, NO_GAP_INC);
		if( init_algns[cur_id].sign == 0 ) {
			e = find_yloc_one(init_algns[cur_id], f, width(init_algns[cur_id].y)-e2, GAP_INC_IN_Y);			
		}
		else if( init_algns[cur_id].sign == 1 ) {
			e = find_yloc_one(init_algns[cur_id], f, e1, GAP_INC_IN_Y);			
		}

		if( b < b1 ) b = b1; // beginning of the alignment
		if( e > e1 ) e = e1; // number of columns

		get_nth_algn(S1, T1, init_algns[cur_id].fid, b, f, a_info, REG);
		pid = cal_pid_maf_beg(S1, T1, b, e);
		res = (int)((pid/((float)100)) * (float)(len));
	}
	init_matrix[entry].val = res;
	
	free(st);
	return(res);
}

void free_graph(struct bipar_node *graph, int num_nodes)
{
	int i = 0;

	for( i = 0; i < num_nodes; i++ ) {
		free(graph[i].adj_list);			
	}
}

void update_init_algns_for_ancestral(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, int num_init_algns)
{
	int i = 0, j = 0;
	int cur_id = 0;
	int entry = 0;

	for( i = 0; i < num_init_algns; i++ ) {
//		init_algns[i].init_sign = init_algns[i].sign;
		init_algns[i].lock = init_algns[i].sign;
		init_algns[i].sign = DELETED;
	}

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			entry = i * num_right + j;
			if( init_matrix[entry].is_in == true ) {
				cur_id = init_matrix[entry].algn_id;	
				init_algns[cur_id].sign = init_algns[cur_id].lock;
			}
		}
	}
}

void init_map_in_dp(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f)
{
	int **scores, **s1, **s2;
	struct point **paths, **p1, **p2; 
	int i = 0;
	struct point **t1, **t2;
	int *id1, *id2;
	int b = -1, e = -1, prev_b = -1;
	int entry = 0;

	id1 = (int *) ckalloc(sizeof(int));
	id2 = (int *) ckalloc(sizeof(int));
	t1 = (struct point **) ckalloc(num_left * sizeof(struct point *));
	t2 = (struct point **) ckalloc(num_left * sizeof(struct point *));
	for( i = 0; i < num_left; i++ ) {
		t1[i] = (struct point *) ckalloc(num_right * sizeof(struct point));
	} // the closest non-zero index
	for( i = 0; i < num_left; i++ ) {
		t2[i] = (struct point *) ckalloc(num_right * sizeof(struct point));
	} // the closest non-zero index

	scores = (int **) ckalloc(num_left * sizeof(int *));
	s1 = (int **) ckalloc(num_left * sizeof(int *));
	s2 = (int **) ckalloc(num_left * sizeof(int *));
	paths = (struct point **) ckalloc(num_left * sizeof(struct point *));
	p1 = (struct point **) ckalloc(num_left * sizeof(struct point *));
	p2 = (struct point **) ckalloc(num_left * sizeof(struct point *));

	for( i = 0; i < num_left; i++ ) {
		scores[i] = (int *) ckalloc(num_right * sizeof(int));
	}

	for( i = 0; i < num_left; i++ ) {
		s1[i] = (int *) ckalloc(num_right * sizeof(int));
	}

	for( i = 0; i < num_left; i++ ) {
		s2[i] = (int *) ckalloc(num_right * sizeof(int));
	}

	for( i = 0; i < num_left; i++ ) {
		paths[i] = (struct point *) ckalloc(num_right * sizeof(struct point));
	}

	for( i = 0; i < num_left; i++ ) {
		p1[i] = (struct point *) ckalloc(num_right * sizeof(struct point));
	}

	for( i = 0; i < num_left; i++ ) {
		p2[i] = (struct point *) ckalloc(num_right * sizeof(struct point));
	}

	*id1 = -1;
	*id2 = -1;
	cal_scores(s1, p1, t1, init_matrix, num_left, num_right, init_algns, f);
	cal_scores_rev(s2, p2, t2, init_matrix, num_left, num_right, init_algns, f);
	combine_scores(s1, s2, p1, p2, t1, t2, scores, paths, num_left, num_right, id1, id2);


	if( ((*id1) != -1) && ((*id2) != -1) ) {
		b = *id1;
		e = *id2;
		i = 0;
		while( (i < (num_left+num_right)) &&  (b != -1) && (e != -1) && (b >= 0) && (b < num_left) && (e >= 0) && (e < num_right) ) {
			entry = b * num_right + e;
			init_matrix[entry].is_in = true;
			prev_b = b; 
			b = paths[b][e].x;
			e = paths[prev_b][e].y;
			i++;
		}
	}

	check_unmapped_nodes(init_matrix, num_left, num_right, init_algns);
	for( i = 0; i < num_left; i++ ) {
		free(scores[i]);
		free(paths[i]);
		free(p1[i]);
		free(p2[i]);
		free(s1[i]);
		free(s2[i]);
		free(t1[i]);
		free(t2[i]);
	}
	free(t1);
	free(t2);
	free(s1);
	free(s2);
	free(p1);
	free(p2);
	free(scores);
	free(paths);
	free(id1);
	free(id2);
}

void check_unmapped_nodes(struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns)
{
	bool *l_matched, *r_matched;
	int *l_weight, *r_weight;
	struct I *l_regs, *r_regs;
	int i = 0, j = 0;
	int entry = 0;
	struct I reg1, reg2;
	int id = -1;
	int lo = 0, hi = 0;

	l_matched = (bool *) ckalloc(sizeof(bool) * num_left);
	r_matched = (bool *) ckalloc(sizeof(bool) * num_right);
	l_weight = (int *) ckalloc(num_left * sizeof(int));
	r_weight = (int *) ckalloc(num_right * sizeof(int));
	l_regs = (struct I *) ckalloc(num_left * sizeof(struct I));
	r_regs = (struct I *) ckalloc(num_right * sizeof(struct I));

	reg1 = assign_I(0, 1);
	reg2 = assign_I(0, 1);
	for( i = 0; i < num_left; i++ ) {
		l_matched[i] = false;
		l_weight[i] = 0;
		l_regs[i] = assign_I(-1,0);
	}
	for( i = 0; i < num_right; i++ ) {
		r_matched[i] = false;
		r_weight[i] = 0;
		r_regs[i] = assign_I(-1,0);
	}

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			entry = i * num_right + j;
			if( init_matrix[entry].is_in == true ) {
				l_weight[i] = l_weight[i] + init_matrix[entry].weight;
				r_weight[j] = r_weight[j] + init_matrix[entry].weight;

				id = init_matrix[entry].algn_id;
				reg1 = assign_I(init_algns[id].x.lower+init_algns[id].xl_diff, init_algns[id].x.upper-init_algns[id].xr_diff);
				reg2 = assign_I(init_algns[id].y.lower+init_algns[id].yl_diff, init_algns[id].y.upper-init_algns[id].yr_diff);

				if( l_regs[i].lower == -1 ) {
					l_regs[i] = assign_I(reg1.lower, reg1.upper);
				}
				else {
					if(l_regs[i].lower < reg1.lower) lo = l_regs[i].lower;
					else lo = reg1.lower;

					if(l_regs[i].upper > reg1.upper) hi = l_regs[i].upper;
					else hi = reg1.upper;

					l_regs[i] = assign_I(lo, hi);
				}

				if( r_regs[j].lower == -1 ) {
					r_regs[j] = assign_I(reg2.lower, reg2.upper);
				}
				else {
					if(r_regs[j].lower < reg2.lower) lo = r_regs[j].lower;
					else lo = reg2.lower;

					if(r_regs[j].upper > reg2.upper) hi = r_regs[j].upper;
					else hi = reg2.upper;

					r_regs[j] = assign_I(lo, hi);
				}
			}
			else {
			}
		}
	}

	mark_matched(l_matched, num_left, num_right, init_matrix, LEFT);
	mark_matched(r_matched, num_right, num_left, init_matrix, RIGHT);

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			entry = i * num_right + j;
			if( (((l_matched[i] == false) && (r_matched[j] == true)) || ((l_matched[i] == true) && (r_matched[j] == false))) && (init_matrix[entry].weight != -1) && (init_matrix[entry].weight > l_weight[i]) && (init_matrix[entry].weight > r_weight[j])) 
			{
				id = init_matrix[entry].algn_id;
				reg1 = assign_I(init_algns[id].x.lower+init_algns[id].xl_diff, init_algns[id].x.upper-init_algns[id].xr_diff);
				reg2 = assign_I(init_algns[id].y.lower+init_algns[id].yl_diff, init_algns[id].y.upper-init_algns[id].yr_diff);
				if( l_matched[i] == true ) {
					if( (strict_almost_equal(reg1, l_regs[i]) == true) || (strict_almost_equal(l_regs[i], reg1) == true) ) {}
					else if( (f_loose_subset(l_regs[i], reg1, STRICT) == true) && (width(l_regs[i]) <= (int)(0.8*((float)width(reg1))) ) ) {
						r_matched[j] = true;
						r_weight[j] = r_weight[j] + init_matrix[entry].weight;
						l_weight[i] = l_weight[i] + init_matrix[entry].weight;

						if(l_regs[i].lower < reg1.lower) lo = l_regs[i].lower;
						else lo = reg1.lower;

						if(l_regs[i].upper > reg1.upper) hi = l_regs[i].upper;
						else hi = reg1.upper;

						l_regs[i] = assign_I(lo, hi);
						r_regs[j] = assign_I(reg2.lower, reg2.upper);
						init_matrix[entry].is_in = true;
					}
				}
				else if( r_matched[j] == true ) {
					if( (strict_almost_equal(reg2, r_regs[j]) == true) || (strict_almost_equal(r_regs[j], reg2) == true) ) {}
					else if( (f_loose_subset(r_regs[j], reg2, STRICT) == true) && (width(r_regs[j]) <= (int)(0.8*((float)width(reg2))) ) ) {
						l_matched[i] = true;
						r_weight[j] = r_weight[j] + init_matrix[entry].weight;
						l_weight[i] = l_weight[i] + init_matrix[entry].weight;

						if(r_regs[j].lower < reg2.lower) lo = r_regs[j].lower;
						else lo = reg2.lower;

						if(r_regs[j].upper > reg2.upper) hi = r_regs[j].upper;
						else hi = reg2.upper;

						r_regs[j] = assign_I(lo, hi);
						l_regs[i] = assign_I(reg1.lower, reg1.upper);
						init_matrix[entry].is_in = true;
					}
				}
			}
		}
	}

	free(r_regs);
	free(l_regs);
	free(r_matched);
	free(l_matched);
	free(r_weight);
	free(l_weight);
}

int find_max_score(int **s1, struct point **p1, struct point **t1, int b1, int e1, int b2, int e2, int *id1, int *id2, int *end_id1, int *end_id2, int mode)	
{
	int i = 0, j = 0;
	int max_score = 0;
	int max_id1 = -1, max_id2 = -1;
	int max_end_id1 = -1, max_end_id2 = -1;
	int b = -1, e = -1;
	int prev_b = -1, prev_e = -1;
	int tmp_val = 0;

	j = e2;
	for( i = b1; i <= e1; i++ ) {
		prev_b = -1;
		prev_e = -1;
		tmp_val = 0;
		if( (t1[i][j].x < b1) || (t1[i][j].x > e1) || (t1[i][j].y < b2) || (t1[i][j].y > e2 ) ) {}
		else {
			b = t1[i][j].x;
			e = t1[i][j].y;
			tmp_val = s1[b][e];

			while( (b >= b1) && (b <= e1) && (e >= b2) && (e <= e2) ) {
				prev_b = b;
				prev_e = e;
				b = p1[b][e].x;
				e = p1[prev_b][e].y;
			}

			if( (b <= -1) || (e <= -1) ) {}
			else {
				if( s1[b][e] != 0 ) {
					tmp_val = tmp_val - s1[b][e];
				}
			}

			if( tmp_val == 0 ) {}
			else if( max_score < tmp_val ) {
				max_score = tmp_val; 
				max_id1 = t1[i][j].x;
				max_id2 = t1[i][j].y;
				max_end_id1 = prev_b;
				max_end_id2 = prev_e;
			}
		}
	}

	if( mode == FIRST_RUN ) {
		i = e1;
	}
	else if( mode == SECOND_RUN ) {
		i = b1;
	}
	for( j = b2; j <= e2; j++ ) {
		prev_b = -1;
		prev_e = -1;
		tmp_val = 0;
		if( (t1[i][j].x < b1) || (t1[i][j].x > e1) || (t1[i][j].y < b2) || (t1[i][j].y > e2 ) ) {}
		else {
			b = t1[i][j].x;
			e = t1[i][j].y;
			tmp_val = s1[b][e];

			while( (b >= b1) && (b <= e1) && (e >= b2) && (e <= e2) ) {
				prev_b = b;
				prev_e = e;
				b = p1[b][e].x;
				e = p1[prev_b][e].y;
			}

			if( (b <= -1) || (e <= -1) ) {}
			else {
				if( s1[b][e] != 0 ) {
					tmp_val = tmp_val - s1[b][e];
				}
			}

			if( tmp_val == 0 ) {}
			else if( max_score < tmp_val ) {
				max_score = tmp_val; 
				max_id1 = t1[i][j].x;
				max_id2 = t1[i][j].y;
				max_end_id1 = prev_b;
				max_end_id2 = prev_e;
			}
		}
	}

	*id1 = max_id1;
	*id2 = max_id2;
	*end_id1 = max_end_id1;
	*end_id2 = max_end_id2;
	return(max_score);
}

void combine_scores(int **s1, int **s2, struct point **p1, struct point **p2, struct point **t1, struct point **t2, int **scores, struct point **paths, int num_left, int num_right, int *res_id1, int *res_id2)
{
	int score1 = 0, score2 = 0;
	int *id1, *id2, *rev_id1, *rev_id2;
	int *end_id1, *end_id2, *rev_end_id1, *rev_end_id2;
	int max_id1 = -1, max_id2 = -1, max_end_id1 = -1, max_end_id2 = -1;
	int prev_id1 = -1, prev_id2 = -1;
	int start_id1 = -1, start_id2 = -1;
	bool *is_filled1, *is_filled2;
	bool is_done = false, is_done_pos = false, is_done_rev = false;
	int b1 = 0, e1 = num_left-1, b2 = 0, e2 = num_right-1;
	int old_b1 = 0, old_e1 = 0, old_b2 = 0, old_e2 = 0;
	int i = 0, j = 0;
	bool is_rev = false;
	int temp1 = 0, temp2 = 0;

	id1 = (int *) ckalloc(sizeof(int));
	id2 = (int *) ckalloc(sizeof(int));
	rev_id1 = (int *) ckalloc(sizeof(int));
	rev_id2 = (int *) ckalloc(sizeof(int));
	end_id1 = (int *) ckalloc(sizeof(int));
	end_id2 = (int *) ckalloc(sizeof(int));
	rev_end_id1 = (int *) ckalloc(sizeof(int));
	rev_end_id2 = (int *) ckalloc(sizeof(int));
	is_filled1 = (bool *) ckalloc(num_left * sizeof(bool));
	is_filled2 = (bool *) ckalloc(num_right * sizeof(bool));

	*id1 = -1;
	*id2 = -1;
	*rev_id1 = -1;
	*rev_id2 = -1;
	*end_id1 = -1;
	*end_id2 = -1;
	*rev_end_id1 = -1;
	*rev_end_id2 = -1;

	init_bool_list(is_filled1, num_left);
	init_bool_list(is_filled2, num_right);

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			if( s1[i][j] != 0) {
				scores[i][j] = s1[i][j];
				paths[i][j].x = p1[i][j].x;
				paths[i][j].y = p1[i][j].y;
			}
			else if (s2[i][j] != 0) 
			{
				scores[i][j] = s2[i][j];
				paths[i][j].x = p2[i][j].x;
				paths[i][j].y = p2[i][j].y;
			}
			else {
				scores[i][j] = 0;
				paths[i][j].x = -1;
				paths[i][j].y = -1;
			}
		}
	}

	i = 0;
	while( (i < (num_right + num_left)) && (is_done == false) ) {
		score1 = 0;
		score2 = 0;
		max_id1 = -1;
		max_id2 = -1;
//		prev_id1 = max_end_id1;
//		prev_id2 = max_end_id2;
		max_end_id1 = -1;
		max_end_id2 = -1;
		is_rev = false;
		if( is_done_pos == false ) {
			score1 = find_max_score(s1, p1, t1, b1, e1, b2, e2, id1, id2, end_id1, end_id2, FIRST_RUN);	
			if( (score1 == 0) || ((*id1) == -1) || ((*id2) == -1) ) {
				is_done_pos = true;
			}	
		}

		if( is_done_rev == false ) {
			score2 = find_max_score(s2, p2, t2, b1, e1, b2, e2, rev_id1, rev_id2, rev_end_id1, rev_end_id2, SECOND_RUN);	
			if( (score2 == 0) || ((*rev_id1) == -1) || ((*rev_id2) == -1) ) {
				is_done_rev = true;
			}
		}

		if( (score1 == 0) && (score2 == 0) ) {
			is_done = true;
		}
		else if( score1 == 0 ) {
			is_done_pos = true;
			is_rev = true;
			max_id1 = *rev_id1;
			max_id2 = *rev_id2;
			max_end_id1 = *rev_end_id1;
			max_end_id2 = *rev_end_id2;
		}
		else if( score2 == 0 ) {
			max_id1 = *id1;
			max_id2 = *id2;
			max_end_id1 = *end_id1;
			max_end_id2 = *end_id2;
			is_done_rev = true;
		}
		else if( score1 >= score2 ) {
			max_id1 = *id1;
			max_id2 = *id2;
			max_end_id1 = *end_id1;
			max_end_id2 = *end_id2;
		}
		else {
			is_rev = true;
			max_id1 = *rev_id1;
			max_id2 = *rev_id2;
			max_end_id1 = *rev_end_id1;
			max_end_id2 = *rev_end_id2;
		}

		if( ((start_id1 == -1) || (start_id2 == -1)) && (max_id1 != -1) && (max_id2 != -1) ) {
			start_id1 = max_id1;
			start_id2 = max_id2;
		}
		else if( (start_id1 != -1) && (start_id2 != -1) && (max_id1 != -1)  && (max_id2 != -1) && (max_end_id1 != -1) && (max_end_id2 != -1)) {
			if((max_id2 > start_id2 ) && (max_end_id2 >= start_id2)) {
				temp1 = max_id1;
				temp2 = max_id2;
				max_id1 = start_id1;
				max_id2 = start_id2;
				start_id1 = temp1;
				start_id2 = temp2;

				temp1 = max_end_id1;
				temp2 = max_end_id2;
				max_end_id1 = prev_id1;
				max_end_id2 = prev_id2;
				prev_id1 = temp1;
				prev_id2 = temp2;
			}
		}

		if( (is_done == true) || (max_id1 == -1) || (max_id2 == -1) ) {
		}
		else {
			update_matrix_paths(scores, paths, start_id1, start_id2, prev_id1, prev_id2, max_id1, max_id2, max_end_id1, max_end_id2, is_filled1, is_filled2, num_left, num_right);
		}

		if( is_done == true ) {}
		else if( (is_done_pos == true) && (is_done_rev == true)  ) {
			is_done = true;
		}
		else {
			old_b1 = b1;
			old_b2 = b2;
			old_e1 = e1;
			old_e2 = e2;
			is_done = find_next_bound(paths, start_id1, start_id2, is_filled1, is_filled2, num_left, num_right, &b1, &e1, &b2, &e2, &prev_id1, &prev_id2);

			if( (old_b1 == b1) && (old_b2 == b2) && (old_e1 == e1) && (old_e2 == e2) ) {
				is_done = true;
			}
		}

		i++;
	}

	*res_id1 = start_id1;
	*res_id2 = start_id2;
	free(is_filled1);
	free(is_filled2);
	free(end_id1);
	free(end_id2);
	free(rev_end_id1);
	free(rev_end_id2);
	free(rev_id1);
	free(rev_id2);
	free(id1);
	free(id2);
}

bool find_next_bound(struct point **paths, int start_id1, int start_id2, bool *is_filled1, bool *is_filled2, int num_left, int num_right, int *b1, int *e1, int *b2, int *e2, int *prev_id1, int *prev_id2)
{
	int i = 0;
	int from = 0, old_from = -1, old_to = -1, to = 0;
	int max_id1 = 0, max_id2 = 0;
	bool res = false;

	i = 0;
	while(i < num_left) {
		while( (i < num_left ) && (is_filled1[i] == true) ) i++;
		
		if( i < num_left ) from = i;
		else from = num_left-1;

		while( (i < num_left ) && (is_filled1[i] == false) ) i++;
		if( i < num_left ) to = (i-1);
		else to = num_left-1;

		if( (to - from + 1) > (max_id2 - max_id1 + 1) ) {
			max_id1 = from;
			max_id2 = to;
		}
		i++;	
	}

	if( ( max_id2 - max_id1 + 1 ) < 2 ) {
		res = true;
	}
	else {
		*b1 = max_id1;
		*e1 = max_id2;

		max_id1 = 0;
		max_id2 = 0;
		i = 0;
		while(i < num_right) {
			while( (i < num_right ) && (is_filled2[i] == true) ) i++;
			
			if( i < num_right ) from = i;
			else from = num_right-1;

			while( (i < num_right ) && (is_filled2[i] == false) ) i++;
			if( i < num_right ) to = (i-1);
			else to = num_right-1;

			if( (to - from + 1) > (max_id2 - max_id1 + 1) ) {
				max_id1 = from;
				max_id2 = to;
			}
			i++;	
		}

		if( ( max_id2 - max_id1 + 1 ) < 2 ) {
			res = true;
		}
		else {
			*b2 = max_id1;
			*e2 = max_id2;
		}
	}

	if( res == false ) {
		from = start_id1;
		to = start_id2;
		i = 0;
		while( (i <= (num_left+num_right) ) && (from != -1) && (to != -1) && (is_still_in_the_same_part(start_id1, start_id2, from, to, *b1, *e1, *b2, *e2) == true)) {
			*prev_id1 = from;
			*prev_id2 = to;
			old_from = from;
			old_to = to;
			from = paths[old_from][old_to].x;
			to = paths[old_from][old_to].y;

			if( (from != -1) && ( to != -1) && (is_already_visited(old_from, old_to, paths, num_left, num_right, start_id1, start_id2, i+1) == true) )
			{
				paths[old_from][old_to].x = -1;
				paths[old_from][old_to].y = -1;
				from = -1;
				to = -1;
			}

			if( (from == old_from) && (to == old_to) ) {
				from = -1;
				to = -1;
			}
			i++;
		}
	}
	else {
		*prev_id1 = -1;
		*prev_id2 = -1;
	}
	
	return(res);
}

bool is_already_visited(int b, int e, struct point **paths, int num_left, int num_right, int s1, int s2, int cur_run)
{
	bool res = false;
	int old_t1 = s1, old_t2 = s2;
	int t1 = s1, t2 = s2;
	int i = 0;
	int cur_b = -1, cur_e = -1;

	cur_b = paths[b][e].x;
	cur_e = paths[b][e].y;
	while( (i <= (num_left+num_right) && (t1 != cur_b) && (t2 != cur_e) ) ) {
		old_t1 = t1;
		old_t2 = t2;
		t1 = paths[old_t1][old_t2].x;
		t2 = paths[old_t1][old_t2].y;
		i++;
	}

	if( (old_t1 == b) && (old_t2 == e) ) {
		if( i != cur_run ) {
			res = true;
		}
	}
	else {
		res = true;
	}

	if( i > (num_left+num_right) ) res = true;

	return(res);
}

void update_matrix_paths(int **s, struct point **p, int start_id1, int start_id2, int prev_id1, int prev_id2, int b1, int e1, int b2, int e2, bool *is_filled1, bool *is_filled2, int num_left, int num_right)
{
	int i = 0, j = 0;
	int temp_b = -1, temp_e = -1, b = -1, e = -1;
	int tmp_val = 0;
	int last_id1 = -1, last_id2 = -1;
	int last_score = 0;

	if( (b1 == -1) || (e1 == -1) ) {
		fatalf("unexpected negative value %d %d\n", b1, e1);
	}

	b = b1;
	e = e1;
	tmp_val = s[b][e];
	while( (i < num_left) && (i < num_right) && (b != -1) && (e != -1) && (b >= 0) && (e >= 0 ) && (b < num_left) && (e < num_right) && ((b != b2) || (e != e2))) 
	{
		is_filled1[b] = true;
		is_filled2[e] = true;
		temp_b = b;
		temp_e = e;
		b = p[b][e].x;
		e = p[temp_b][e].y;

		if( (temp_b != -1) && (temp_e != -1) && (temp_b >= 0) && (temp_e >= 0 ) && (temp_b < num_left) && (temp_e < num_right) && (b != -1) && (e != -1) && (b >= 0) && (e >= 0 ) && (b < num_left) && (e < num_right) && ((b != b2) || (e != e2)) ) {
			if( abs(b-temp_b) <= 2 ) {
				if( temp_e > e ) {
					for( j = e; j < temp_e; j++ ) {
						is_filled2[j] = true;
					}
				}
				else {
					for( j = temp_e; j < e; j++ ) {
						is_filled2[j] = true;
					}
				}
			}

			if( abs(e-temp_e) <= 2 ) {
				if( temp_b > b ) {
					for( j = b; j < temp_b; j++ ) {
						is_filled1[j] = true;
					}
				}
				else {
					for( j = temp_b; j < b; j++ ) {
						is_filled1[j] = true;
					}
				}
			}
		}
		i++;
	}

	if( (b == -1) || (e == -1) ) {}
	else if((temp_b < 0) || (temp_e < 0 ) || (temp_b >= num_left) || (temp_e > num_right) || (b < 0) || (e < 0 ) || (b >= num_left) || (e >= num_right)) {}
	else if( (temp_b == -1) || (temp_e == -1) || (p[b][e].x == -1) || (p[b][e].y == -1) ) {
		is_filled1[b] = true;
		is_filled2[e] = true;
	}
	else {
		if( abs(b-temp_b) <= 2 ) {
			if( temp_e > e ) {
				for( j = e; j < temp_e; j++ ) {
					is_filled2[j] = true;
				}
			}
			else {
				for( j = temp_e; j < e; j++ ) {
					is_filled2[j] = true;
				}
			}
		}

		if( abs(e-temp_e) <= 2 ) {
			if( temp_b > b ) {
				for( j = b; j < temp_b; j++ ) {
					is_filled1[j] = true;
				}
			}
			else {
				for( j = temp_b; j < b; j++ ) {
					is_filled1[j] = true;
				}
			}
		}
		is_filled1[b] = true;
		is_filled2[e] = true;
		tmp_val = tmp_val - s[p[b][e].x][p[b][e].y];	
	}

	if( (prev_id1 == -1) || (prev_id2 == -1) ) {
		if( (start_id1 == b1) && (start_id2 == e1) ) {
		}
		else {
			fatalf("%d and %d are supposed to be identical to %d and %d\n", start_id1, start_id2, b1, e1);
		}
	}
	else {
		b = start_id1;
		e = start_id2;
		while( (i < (num_left + num_right)) && (b != -1) && (e != -1) && (b >= 0) && (e >= 0 ) && (b < num_left) && (e < num_right) && ((b != prev_id1) || (e != prev_id2))) 
		{
			is_filled1[b] = true;
			is_filled2[e] = true;
			s[b][e] = s[b][e] + tmp_val;
			temp_b = b;
			b = p[b][e].x;
			e = p[temp_b][e].y;
			i++;
		}

		if( (b == -1) || (e == -1) || (b != prev_id1) || (e != prev_id2) ) {
			fatalf("maximal path doesn't pass [%d,%d]\n", prev_id1, prev_id2);
		}
		else {
			is_filled1[b] = true;
			is_filled2[e] = true;
			s[b][e] = s[b][e] + tmp_val;

			if( (p[prev_id1][prev_id2].x != -1) && (p[prev_id1][prev_id2].y != -1) ) {
				last_id1 = p[prev_id1][prev_id2].x;
				last_id2 = p[prev_id1][prev_id2].y;
				last_score = s[last_id1][last_id2];

				i = 0;
				b = b1;
				e = e1;
				tmp_val = s[b][e];
				while( (i < (num_left+num_right)) && (b != -1) && (e != -1) && (b >= 0) && (e >= 0 ) && (b < num_left) && (e < num_right) && ((b != b2) || (e != e2))) 
				{
					s[b][e] = s[b][e] + last_score;
					temp_b = b;
					b = p[b][e].x;
					e = p[temp_b][e].y;
					i++;
				}

				if( (b == -1) || (e == -1) || (b != b2) || (e != e2) ) {
					fatalf("maximal path doesn't pass [%d,%d]\n", b, e);
				}
				else {
					s[b][e] = s[b][e] + last_score;
					p[b][e].x = last_id1;
					p[b][e].y = last_id2;
				}
			}

			p[prev_id1][prev_id2].x = b1;
			p[prev_id1][prev_id2].y = e1;
		}	
	}
}

void cal_scores_rev(int **s2, struct point **p2, struct point **temp, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f)
{
	int i = 0, j = 0, k = 0;
	int num = 0;
	int b = 0, e = 0;
	int entry = 0;

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			entry = i * num_right + j;
			temp[i][j].x = -1;
			temp[i][j].y = -1;
			if( init_matrix[entry].sign != 1 ) {
				s2[i][j] = 0;
			}
			else if( init_matrix[entry].tmp_val == -1 ) { 
				s2[i][j] = 0;
			}
			else {
				s2[i][j] = init_matrix[entry].tmp_val;
				temp[i][j].x = i;
				temp[i][j].y = j;
			}

			p2[i][j].x = -1;
			p2[i][j].y = -1;
		}
	}

	num = num_left + num_right - 1;
	for( i = 0; i < num; i++ )
	{
		if ( (num_left - 1 - i) < 0 ) b = 0;
		else b = num_right - 1 - i;

		if ( (b+i) > (num_left-1) ) e = num_left-1;
		else e = b+i;

		for( j = b; j <= e; j++ )
		{
			k = j - (num_left - i - 1);
			if( ((j < 0)) || (j > (num_left-1)) || (k < 0) || (k > (num_right-1)) ) {}
			else {
				update_score_in_matrix_rev(s2, p2, temp, j, k, init_matrix, num_left, num_right, init_algns, f);		
			}
		}
	}
}

void cal_scores(int **s1, struct point **p1, struct point **temp, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f)
{
	int i = 0, j = 0, k = 0;
	int num = 0;
	int b = 0, e = 0;
	int entry = 0;

	for( i = 0; i < num_left; i++ ) {
		for( j = 0; j < num_right; j++ ) {
			entry = i * num_right + j;
			temp[i][j].x = -1;
			temp[i][j].y = -1;
			if( init_matrix[entry].sign != 0 ) {
				s1[i][j] = 0;
			}
			else if( init_matrix[entry].tmp_val == -1 ) { 
				s1[i][j] = 0;
			}
			else {
				s1[i][j] = init_matrix[entry].tmp_val;
				temp[i][j].x = i;
				temp[i][j].y = j;
			}

			p1[i][j].x = -1;
			p1[i][j].y = -1;
		}
	}

	num = num_left + num_right - 1;
	for( i = 0; i < num; i++ )
	{
		if ( (i - num_right + 1) < 0 ) b = 0;
		else b = i-num_right + 1;

		if ( (b+i) > (num_left-1) ) e = num_left-1;
		else e = b+i;

		for( j = b; j <= e; j++ )
		{
			k = i - j;
			if( (j < 0) || (j > (num_left-1)) || (k < 0) || (k > (num_right-1)) ) {}
			else {
				update_score_in_matrix(s1, p1, temp, j, k, init_matrix, num_left, num_right, init_algns, f);		
			}
		}
	}
}

void update_score_in_matrix_rev(int **s1, struct point **p1, struct point **temp, int id1, int id2, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f)
{
	int i = 0, j = 0;
	int max_val = 0;
	int max_id1 = -1, max_id2 = -1;
	int algn_id1 = -1, algn_id2 = -1;
	int entry = 0;
	int cur_entry = 0;
	float gap_penalty = (float)0;
	int b = 0, e = 0;
	int tmp_val = 0;
	int temp_entry = 0;
	int temp_id = -1;
	
	cur_entry = id1 * num_right + id2;
	algn_id2 = init_matrix[cur_entry].algn_id;

	if( algn_id2 == -1 ) {
		max_val = 0;
	}
	else if( s1[id1][id2] != 0 ) {
		max_val = s1[id1][id2];
	}

	for( j = (id2 - 1); j <= id2; j++ ) {
		if( (j < 0) || (j >= num_right) ) {
		}
		else {
			for( i = (id1+1); i < num_left ; i++ ) {
				temp_entry = i * num_right + j;
				temp_id = init_matrix[temp_entry].algn_id;
				if( (temp[i][j].x == -1) || (temp[i][j].y == -1) ) {}
				else if( (j == id2) && (temp_id == -1) ) {}
				else {
					b = temp[i][j].x;
					e = temp[i][j].y;

					if( (b <= id1) || (e > id2) ) {
						fatalf("%d should be larger than %d and %d smaller than %d\n", b, id1, e, id2);
					}

					if( (algn_id2 != -1) && (temp_id != -1) && (j == id2) && ((f_loose_overlap(init_algns[temp_id].m_x, init_algns[algn_id2].m_x, LOOSE) == true) || (f_loose_overlap(init_algns[temp_id].m_y, init_algns[algn_id2].m_y, LOOSE) == true)) ) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (temp_id != -1) && (j == id2) && ((init_algns[temp_id].m_y.lower >= init_algns[algn_id2].m_y.lower) && (init_algns[temp_id].m_y.upper >= init_algns[algn_id2].m_y.upper) )) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (s1[id1][id2] != 0) ) {
						entry = b * num_right + e;
						algn_id1 = init_matrix[entry].algn_id;	
						gap_penalty = 1.00 - (0.01 * (float)(abs(id1-1-b)+abs(id2-1-e)));
						tmp_val = (int)(gap_penalty * (float)s1[b][e]) + cal_non_merged_score(init_algns, algn_id1, algn_id2, f, s1[id1][id2]);
					}
					else {
						tmp_val = s1[b][e];
					}
	
					if( tmp_val > max_val ) {
						max_val = tmp_val;
						max_id1 = b;
						max_id2 = e;
					}
				}
			}
		}
	}

	for( i = (id1 + 1); i >= id1; i-- ) {
		if( (i < 0) || (i >= num_left) ) {
		}
		else {
			for( j = 0; j < (id2-1); j++ ) {
				temp_entry = i * num_right + j;
				temp_id = init_matrix[temp_entry].algn_id;
				if( (temp[i][j].x == -1) || (temp[i][j].y == -1) ) {}
				else if( (i == id1) && (temp_id == -1) ) {}
				else {
					b = temp[i][j].x;
					e = temp[i][j].y;

					if( (b < id1) || (e >= id2) ) {
						fatalf("%d should be larger than %d and %d smaller than %d\n", b, id1, e, id2);
					}
	
					if( (algn_id2 != -1) && (temp_id != -1) && (i == id1) && ((f_loose_overlap(init_algns[temp_id].m_x, init_algns[algn_id2].m_x, LOOSE) == true) || (f_loose_overlap(init_algns[temp_id].m_y, init_algns[algn_id2].m_y, LOOSE) == true)) ) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (temp_id != -1) && (i == id1) && ((init_algns[temp_id].m_y.lower >= init_algns[algn_id2].m_y.lower) && (init_algns[temp_id].m_y.upper >= init_algns[algn_id2].m_y.upper) )) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (s1[id1][id2] != 0) ) {
						entry = b * num_right + e;
						algn_id1 = init_matrix[entry].algn_id;	
						gap_penalty = 1.00 - (0.01 * (float)(abs(id1-1-b)+abs(id2-1-e)));
						tmp_val = (int)(gap_penalty * (float)s1[b][e]) + cal_non_merged_score(init_algns, algn_id1, algn_id2, f, s1[id1][id2]);
					}
					else {
						tmp_val = s1[b][e];
					}

					if( tmp_val > max_val ) {
						max_val = tmp_val;
						max_id1 = b;
						max_id2 = e;
					}
				}
			}
		}
	}

	if( (algn_id2 == -1) || (s1[id1][id2] == 0) ) {
		temp[id1][id2].x = max_id1;
		temp[id1][id2].y = max_id2;
		p1[id1][id2].x = -1;
		p1[id1][id2].y = -1;
	}
	else {
		temp[id1][id2].x = id1;
		temp[id1][id2].y = id2;
		p1[id1][id2].x = max_id1;
		p1[id1][id2].y = max_id2;
		if( max_val != 0 ) {
			s1[id1][id2] = max_val;
		}
	}
}

void update_score_in_matrix(int **s1, struct point **p1, struct point **temp, int id1, int id2, struct matrix_ent *init_matrix, int num_left, int num_right, struct DotList *init_algns, FILE *f)
{
	int i = 0, j = 0;
	int max_val = 0;
	int max_id1 = -1, max_id2 = -1;
	int algn_id1 = -1, algn_id2 = -1;
	int entry = 0;
	int cur_entry = 0;
	float gap_penalty = (float)0;
	int b = 0, e = 0;
	int tmp_val = 0;
	int temp_entry = 0;
	int temp_id = -1;
	
	cur_entry = id1 * num_right + id2;
	algn_id2 = init_matrix[cur_entry].algn_id;

	if( algn_id2 == -1 ) {
		max_val = 0;
	}
	else if( s1[id1][id2] != 0 ) {
		max_val = s1[id1][id2];
	}

	for( j = (id2 - 1); j <= id2; j++ ) {
		if( (j < 0) || (j >= num_right) ) {
		}
		else {
			for( i = 0; i < id1; i++ ) {
				temp_entry = i * num_right + j;
				temp_id = init_matrix[temp_entry].algn_id;
				if( (temp[i][j].x == -1) || (temp[i][j].y == -1) ) {}
				else if( (j == id2) && (temp_id == -1) ) {}
				else {
					b = temp[i][j].x;
					e = temp[i][j].y;

					if( (b >= id1) || (e > id2) ) {
						fatalf("%d should be smaller than %d and %d than %d\n", b, id1, e, id2);
					}

					if( (algn_id2 != -1) && (temp_id != -1) && (j == id2) && ((f_loose_overlap(init_algns[temp_id].m_x, init_algns[algn_id2].m_x, LOOSE) == true) || (f_loose_overlap(init_algns[temp_id].m_y, init_algns[algn_id2].m_y, LOOSE) == true)) ) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (temp_id != -1) && (j == id2) && ((init_algns[temp_id].m_y.lower >= init_algns[algn_id2].m_y.lower) && (init_algns[temp_id].m_y.upper >= init_algns[algn_id2].m_y.upper) )) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (s1[id1][id2] != 0) ) {
						entry = b * num_right + e;
						algn_id1 = init_matrix[entry].algn_id;	
						gap_penalty = 1.00 - (0.01 * (float)(abs(id1-1-b)+abs(id2-1-e)));
						tmp_val = (int)(gap_penalty * (float)s1[b][e]) + cal_non_merged_score(init_algns, algn_id1, algn_id2, f, s1[id1][id2]);
					}
					else {
						tmp_val = s1[b][e];
					}

					if( tmp_val > max_val ) {
						max_val = tmp_val;
						max_id1 = b;
						max_id2 = e;
					}
				}
			}
		}	
	}

	for( i = (id1 - 1); i <= id1; i++ ) {
		if( (i < 0) || (i >= num_left) ) {
		}
		else {
			for( j = 0; j < (id2-1); j++ ) {
				temp_entry = i * num_right + j;
				temp_id = init_matrix[temp_entry].algn_id;
				if( (temp[i][j].x == -1) || (temp[i][j].y == -1) ) {}
				else if( (i == id1) && (temp_id == -1) ) {}
				else {
					b = temp[i][j].x;
					e = temp[i][j].y;
	
					if( (b > id1) || (e >= id2) ) {
						fatalf("%d should be smaller than %d and %d than %d\n", b, id1, e, id2);
					}

					if( (algn_id2 != -1) && (temp_id != -1) && (i == id1) && ((f_loose_overlap(init_algns[temp_id].m_x, init_algns[algn_id2].m_x, LOOSE) == true) || (f_loose_overlap(init_algns[temp_id].m_y, init_algns[algn_id2].m_y, LOOSE) == true)) ) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (temp_id != -1) && (j == id2) && ((init_algns[temp_id].m_y.lower >= init_algns[algn_id2].m_y.lower) && (init_algns[temp_id].m_y.upper >= init_algns[algn_id2].m_y.upper) )) {
						tmp_val = 0;
					}
					else if( (algn_id2 != -1) && (s1[id1][id2] != 0) ) {
						entry = b * num_right + e;
						algn_id1 = init_matrix[entry].algn_id;	
						gap_penalty = 1.00 - (0.01 * (float)(abs(id1-1-b)+abs(id2-1-e)));
						tmp_val = (int)(gap_penalty * (float)s1[b][e]) + cal_non_merged_score(init_algns, algn_id1, algn_id2, f, s1[id1][id2]);
					}
					else {
						tmp_val = s1[b][e];
					}

					if( tmp_val > max_val ) {
						max_val = tmp_val;
						max_id1 = b;
						max_id2 = e;
					}
				}
			}
		}
	}

	if( (algn_id2 == -1) || (s1[id1][id2] == 0) ) {
		temp[id1][id2].x = max_id1;
		temp[id1][id2].y = max_id2;
		p1[id1][id2].x = -1;
		p1[id1][id2].y = -1;
	}
	else {
		temp[id1][id2].x = id1;
		temp[id1][id2].y = id2;
		p1[id1][id2].x = max_id1;
		p1[id1][id2].y = max_id2;
		if( max_val != 0 ) {
			s1[id1][id2] = max_val;
		}
	}
}

int cal_non_merged_score(struct DotList *init_algns, int id1, int id2, FILE *f, int score)
{
	int tmp_val = 0;
	struct I tmp_reg1, tmp_reg2;
	struct I cur_reg1, cur_reg2;
	int beg1 = 0, end1 = 0, len1 = 0;
	int cur_len = 0;
	int end_pos = 0;
	int *t_b;
	float pid = (float)0;
	bool is_x = true;
	char S1[BIG], T1[BIG];
	struct b_list *a_info;

	strcpy(S1, "");
	strcpy(T1, "");

	t_b = (int *) ckalloc(sizeof(int));
	a_info = (struct b_list *) ckalloc(sizeof(struct b_list));

	cur_reg1 = assign_I(0, 1);
	tmp_reg1 = assign_I(0, 1);
	cur_reg2 = assign_I(0, 1);
	tmp_reg2 = assign_I(0, 1);
	
	cur_reg1 = assign_I(init_algns[id1].x.lower+init_algns[id1].xl_diff, init_algns[id1].x.upper-init_algns[id1].xr_diff);
	tmp_reg1 = assign_I(init_algns[id2].x.lower+init_algns[id2].xl_diff, init_algns[id2].x.upper-init_algns[id2].xr_diff);

	cur_reg2 = assign_I(init_algns[id1].y.lower+init_algns[id1].yl_diff, init_algns[id1].y.upper-init_algns[id1].yr_diff);
	tmp_reg2 = assign_I(init_algns[id2].y.lower+init_algns[id2].yl_diff, init_algns[id2].y.upper-init_algns[id2].yr_diff);

	if( (proper_overlap(tmp_reg1, cur_reg1) == true) && (proper_overlap(tmp_reg2, cur_reg2) == true) ) {
		if( width(intersect(tmp_reg1, cur_reg1)) < width(intersect(tmp_reg2, cur_reg2)) ) {
			cur_reg1 = assign_I(init_algns[id1].y.lower+init_algns[id1].yl_diff, init_algns[id1].y.upper-init_algns[id1].yr_diff);
			tmp_reg1 = assign_I(init_algns[id2].y.lower+init_algns[id2].yl_diff, init_algns[id2].y.upper-init_algns[id2].yr_diff);
			is_x = false;
		}
	}
	else if(proper_overlap(tmp_reg2, cur_reg2) == true) {
		cur_reg1 = assign_I(init_algns[id1].y.lower+init_algns[id1].yl_diff, init_algns[id1].y.upper-init_algns[id1].yr_diff);
		tmp_reg1 = assign_I(init_algns[id2].y.lower+init_algns[id2].yl_diff, init_algns[id2].y.upper-init_algns[id2].yr_diff);
		is_x = false;
	}

	if( subset(tmp_reg1, cur_reg1) == true ) {
		tmp_val = 0;
	}
	else if( proper_overlap(tmp_reg1, cur_reg1) == true ) {
		beg1 = init_algns[id2].xl_diff + init_algns[id2].xl_offset;	
		len1 = init_algns[id2].x.upper - init_algns[id2].xr_diff - init_algns[id2].xr_offset - init_algns[id2].x.lower - init_algns[id2].xl_offset;
		end1 = find_xloc_one(init_algns[id2], f, len1, NO_GAP_INC);

		get_nth_algn(S1, T1, init_algns[id2].fid, beg1, f, a_info, REG);
		if( S1[strlen(S1)-1] == '\n' ) end_pos = strlen(S1)-1;
		else end_pos = strlen(S1);

		if( end1 < end_pos ) {
			S1[end1] = '\0';
			end_pos = end1-1;
		}

		if( subset(cur_reg1, tmp_reg1) == true ) {
			pid = (float)init_algns[id2].identity;
			tmp_val = (int)((float)(width(tmp_reg1)-width(cur_reg1)))*(pid/(float)100);
		}
		else {
			if( tmp_reg1.lower < cur_reg1.lower ) {
				tmp_reg2 = assign_I(tmp_reg1.lower, cur_reg1.lower);
			}
			else {
				tmp_reg2 = assign_I(cur_reg1.upper, tmp_reg1.upper);
			}

			if( is_x == true ) {
				cur_len = count_ncol(tmp_reg1, tmp_reg2, S1, end_pos, t_b);
			}
			else if( init_algns[id1].init_sign == 0 ) {
				cur_len = count_ncol(tmp_reg1, tmp_reg2, T1, end_pos, t_b);
			}
			else if( init_algns[id2].init_sign == 1 ) {
				cur_len = count_ncol_rev(tmp_reg1, tmp_reg2, T1, end_pos, t_b);
			}
			else {
				fatalf("an expected case: %d %d\n", id1, id2);
			}

			pid = cal_pid_maf_beg(S1, T1, *t_b, cur_len);

			if( is_x == true ) {
				tmp_val = (int)((float)(width(tmp_reg2)-countN_in_seq(S1, *t_b, cur_len)))*(pid/(float)100);
			}
			else {
				tmp_val = (int)((float)(width(tmp_reg2)-countN_in_seq(T1, *t_b, cur_len)))*(pid/(float)100);
			}
		}
	}
	else {
		tmp_val = score;
	}

	free(t_b);
	free(a_info);
	return(tmp_val);
}

bool is_still_in_the_same_part(int sid1, int sid2, int from, int to, int b1, int e1, int b2, int e2)
{
	bool res = false;

	if( (sid1 == -1) || (sid2 == -1) || (from == -1) || (to == -1) || (b1 == -1) || (e1 == -1) || (b2 == -1) || (e2 == -1) ) {
		res = false;
	}
	else if( (sid1 > e1) && (sid2 > e2) ) {
		if( (from > e1) && (to > e2) ) res = true;
	}
	else if( (sid1 > e1) && (sid2 < b2) ) {
		if( (from > e1) && (to < b2) ) res = true;
	}
	else if( (sid1 < b1) && (sid2 < b2) ) {
		if( (from < b1) && (to < b2) ) res = true;
	}
	else if( (sid1 < b1) && (sid2 > e2) ) {
		if( (from < b1) && (to > e2) ) res = true;
	}
	else {
		res = false;
	}

	return(res);
}
