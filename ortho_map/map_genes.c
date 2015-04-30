#include "main.h"
#include "map_genes.h"
#include "util.h"
#include "util_i.h"
#include "util_input.h"
#include "util_gen.h"
#include "id_ortho_conv.h"
#include "regions.h"

#define MAP_ERR_TH 305

void map_genes(struct DotList *algns, int num_algns, struct exons_list *exons1, int num_exons1, struct exons_list *genes1, int num_genes1, struct exons_list *genes2, int num_genes2, FILE *f)
{
	int i = 0, j = 0, k = 0;
	struct I cur_reg, temp_reg;
	int s_loc = 0, e_loc = 0;
	struct slist *sorted;
	int *res_b, *res_e;
	int count = 0;

	res_b = (int *) ckalloc(sizeof(int));
	res_e = (int *) ckalloc(sizeof(int));
	sorted = (struct slist *) ckalloc(num_algns * sizeof(struct slist));
	for( i = 0; i < num_algns; i++ ) sorted[i].id = i;
	sort_init_algns(sorted, algns, num_algns, INIT_SELF2);

	cur_reg = assign_I(0, 1);
	temp_reg = assign_I(0, 1);

	for( i = 0; i < num_genes2; i++ ) {
		cur_reg = assign_I(genes2[i].reg.lower, genes2[i].reg.upper);	
		s_loc = search_range_b(sorted, algns, num_algns, cur_reg.lower, SELF2);
		e_loc = search_range_e(sorted, algns, num_algns, cur_reg.upper, SELF2);
		count = 0;

		for( j = s_loc; j <= e_loc; j++ ) {
			k = sorted[j].id;
			temp_reg = assign_I(0, 1);
	    *res_b = -1;
 		  *res_e = -1;
    	if((proper_overlap(algns[k].y, cur_reg) == true) && (width(intersect(algns[k].y, cur_reg)) >= DEL_TH)) {
				temp_reg = intersect(algns[k].y, cur_reg);
        find_overlapping_ends(cur_reg, algns[k].y, SELF2, algns, k, f, res_b, res_e);
        if(((*res_b) != -1) && ((*res_e) != -1) && ((*res_b) < (*res_e) )) {
          temp_reg = assign_I(*res_b, *res_e);

					genes2[i].val = get_gene_index(temp_reg, genes1, num_genes1, exons1, num_exons1);
					if( genes2[i].val != -1 ) {
						if( count == 0 ) {
							if( strstr(genes2[i].name, "_ps") != NULL ) {
								printf("\t\\%d", genes2[i].val+1);
							}
							else {
								printf("\t%d", genes2[i].val+1);
							}
						}
						else {
							printf(",%d", genes2[i].val+1);
						}
						count++;
					}
				}	
			}
		}
		if( count > 0 ) {
		}
		else {
			if( strstr(genes2[i].name, "_ps") != NULL ) {
				printf("\t\\%d", -1);
			}
			else {
				printf("\t%d", -1);
			}
		}
	}

	free(sorted);
	free(res_b);
	free(res_e);
}

int get_gene_index(struct I reg, struct exons_list *genes, int num_genes, struct exons_list *exons, int num_exons)
{
  int i = 0, j;
	int max_width = 0;
	bool is_in = false;

	j = -1;
  for( i = 0; i < num_genes; i++ ) {
    if( proper_overlap(genes[i].reg, reg) == true) {
			if( max_width < (width(intersect(genes[i].reg, reg)))) {
				max_width = width(intersect(genes[i].reg, reg));
				j = genes[i].fid;
			}
    }
  }

	if( j != -1 ) {
		if( strstr(genes[j].name, "_ps") != NULL ) {
			if( max_width >= (int)(0.1 * (float)width(genes[j].reg)) ) {
				return(j);
			}
			else return(-1);
		}
		else {
			for( i = 0; i < num_exons; i++ ) {
				if( exons[i].fid == j ) {
					if( (subset(exons[i].reg, reg) == true) || ((proper_overlap(exons[i].reg, reg) == true) && (width(intersect(exons[i].reg, reg)) > (int)(0.8*(float)width(exons[i].reg)) )) ) {
						is_in = true;
					}
				}
			}

			if( is_in == true ) {
  			return(j);
			}
			else {
				return(-1);
			}
		}
	}
	else {
		return(-1);
	}
}

void map_genes_partition(struct DotList *algns, int num_algns, struct exons_list *exons1, int num_exons1, struct exons_list *genes1, int num_genes1, struct exons_list *genes2, int num_genes2, FILE *f)
{
	int i = 0, j = 0, k = 0, n = 0;
	struct I cur_reg, temp_reg, cmp_reg;
	int lo = 0, hi = 0;
	int s_loc = 0, e_loc = 0;
	struct slist *sorted;
	int *res_b, *res_e;
	int count = 0;
	struct int_map *parts;
	int index = 0;
	int val = -1, prev_val = -1, prev_index = -1;
	int *bp;
	int num_bp = 0;
	int num_par = 0;
	bool is_ps_marked = false;
	int num_commas = 0;

	parts = (struct int_map *) ckalloc(num_algns * sizeof(struct int_map));
	bp = (int *) ckalloc((2*num_algns) * sizeof(int));
	res_b = (int *) ckalloc(sizeof(int));
	res_e = (int *) ckalloc(sizeof(int));
	sorted = (struct slist *) ckalloc(num_algns * sizeof(struct slist));
	for( i = 0; i < num_algns; i++ ) {
		sorted[i].id = i;
		parts[i].reg = assign_I(0,1);
		parts[i].bp = assign_I(0,1);
		parts[i].index = -1;
		bp[2*i] = -1;
		bp[2*i+1] = -1;
	}
	sort_init_algns(sorted, algns, num_algns, INIT_SELF2);

	cur_reg = assign_I(0, 1);
	temp_reg = assign_I(0, 1);
	cmp_reg = assign_I(0, 1);

	for( i = 0; i < num_genes2; i++ ) {
		index = 0;
		num_bp = 0;
		cur_reg = assign_I(genes2[i].reg.lower, genes2[i].reg.upper);	
		s_loc = search_range_b(sorted, algns, num_algns, cur_reg.lower, SELF2);
		e_loc = search_range_e(sorted, algns, num_algns, cur_reg.upper, SELF2);
		count = 0;

		val = -1;
		prev_val = -1;

		bp[0] = cur_reg.lower;
		bp[1] = cur_reg.upper;
		num_bp++;
		for( j = s_loc; j <= e_loc; j++ ) {
			k = sorted[j].id;
			temp_reg = assign_I(0, 1);
	    *res_b = -1;
 		  *res_e = -1;
    	if((proper_overlap(algns[k].y, cur_reg) == true) && (width(intersect(algns[k].y, cur_reg)) >= MAP_ERR_TH)) {
				temp_reg = intersect(algns[k].y, cur_reg);
        find_overlapping_ends(cur_reg, algns[k].y, SELF2, algns, k, f, res_b, res_e);
        if(((*res_b) != -1) && ((*res_e) != -1) && ((*res_b) < (*res_e) )) {
          temp_reg = assign_I(*res_b, *res_e);

					prev_val = val;
					val = get_gene_index(temp_reg, genes1, num_genes1, exons1, num_exons1);
					if( val != -1 ) {
						if( index > 0 ) {
							prev_index = -1;
							for( n = 0; n < index; n++ ) {
								if( val == parts[n].index ) {
									prev_index = n;
								}
							}
						}

						if( (index > 0 ) && (prev_index != -1) && ((algns[k].y.lower - parts[prev_index].reg.upper) <= 3*MAP_ERR_TH )) {
							if( algns[k].y.lower < parts[prev_index].reg.lower ) {
								lo = algns[k].y.lower;
							}
							else {
								lo = parts[prev_index].reg.lower;
							}

							if( algns[k].y.upper > parts[prev_index].reg.upper ) {
								hi = algns[k].y.upper;
							}
							else {
								hi = parts[prev_index].reg.upper;
							}
							cmp_reg = assign_I(lo, hi);
							parts[prev_index].reg = intersect(cmp_reg, cur_reg);
							prev_val = parts[prev_index].bp.lower;
							bp[prev_val] = parts[prev_index].reg.lower;
							bp[prev_val+1] = parts[prev_index].reg.upper;
						}
						else {
							parts[index].reg = intersect(algns[k].y, cur_reg);
							parts[index].index = val;
							bp[2*num_bp] = parts[index].reg.lower;
							parts[index].bp.lower = 2*num_bp;
							bp[2*num_bp+1] = parts[index].reg.upper;
							parts[index].bp.upper = 2*num_bp+1;
							index++;
							num_bp++;
						}
					}
				}
			}
		}
		
		if( num_bp > 0 ) {
			quick_sort_num_list_inc(bp, 0, 2*num_bp-1);
		
			num_par = 0;
			for(n = 0; n < (2*num_bp-1); n++) {
				if( (bp[n+1]-bp[n]) > MAP_ERR_TH ) {
					temp_reg = assign_I(bp[n], bp[n+1]-1);

					for( j = 0; j < index; j++ ) {
						if( subset(temp_reg, parts[j].reg) == true ) {
							num_par++;
						}
					}
				}
			}
		}

		is_ps_marked = false;
		if( num_par > 0 ) {
			count = 0;
			for(n = 0; n < (2*num_bp-1); n++) {
				if( (bp[n+1]-bp[n]) > MAP_ERR_TH ) {
					temp_reg = assign_I(bp[n], bp[n+1]-1);

					k = 0;
					num_commas = 0;
					for( j = 0; j < index; j++ ) {
						if( subset(temp_reg, parts[j].reg) == true ) {
							if( (count == 0) && (num_commas == 0) ) {
								if( strstr(genes2[i].name, "_ps") != NULL ) {
									printf("\t\\%d", parts[j].index+1);
									is_ps_marked = true;
								}
								else {
									printf("\t%d", parts[j].index+1);
									is_ps_marked = true;
								}
								num_commas++;
							}
							else {
								if( num_commas == 0 ) {
									printf("%d", parts[j].index+1);
								}
								else {
									printf(",%d", parts[j].index+1);
								}
								num_commas++;
							}
							k++;
						}
					}	

					if( k > 0 ) {
						printf(";");
						count++;
					}
					else {
						if( width(cur_reg) <= 3*MAP_ERR_TH ) {}
						else {
							if( is_ps_marked == false ) {
								if( strstr(genes2[i].name, "_ps") != NULL ) {
									printf("\t\\%d;", -1);
									is_ps_marked = true;
								}
								else {
									printf("\t%d;", -1);
									is_ps_marked = true;
								}
							}
							else {
								printf("-1;");
							}
							count++;
						}
					}
				}
			}
		}

		if( is_ps_marked == false ) {
			if( strstr(genes2[i].name, "_ps") != NULL ) {
				printf("\t\\%d;", -1);
			}
			else {
				printf("\t%d;", -1);
			}
		}
	}

	free(bp);
	free(parts);
	free(sorted);
	free(res_b);
	free(res_e);
}
