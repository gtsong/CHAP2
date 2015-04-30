#include "main.h"
#include "find_ps.h"
#include "util.h"
#include "util_i.h"
#include "util_input.h"
#include "util_gen.h"
#include "id_ortho_conv.h"
#include "regions.h"

void find_pseudogenes(struct exons_list *genes1, int *num_genes1, struct exons_list *genes2, int *num_genes2, struct exons_list *exons2, int *num_exons2, struct DotList *algns, int num_algns, FILE *f, struct exons_list **candi_ps, int *num_candi, int *num_alloc)
{
	struct slist *sorted;
	int i = 0, j = 0;
	int from = 0, to = 0;

	if( num_algns > 0 ) {
		sorted = (struct slist *) ckalloc(num_algns * sizeof(struct slist));
		initialize_slist(sorted, 0, num_algns);
		sort_init_algns(sorted, algns, num_algns, INIT_SELF1);
	}

	for( i = 0; i < (*num_exons2); i++ ) {
		from = i;
		while((i < ((*num_exons2)-1)) && (exons2[i].fid == exons2[i+1].fid)) i++;	
		to = i;
		j = exons2[i].fid;

		if( (j >= (*num_genes2)) || (j < 0) ) {
			fatalf("%d out of gene index range [0-%d]\n", j, (*num_genes2)-1); 
		}
		update_ps_candidates(genes2[j], exons2, from, to, genes1, num_genes1, algns, num_algns, f, candi_ps, num_candi, num_alloc, sorted);
	}

	remove_redun_candi(*candi_ps, *num_candi, genes1, *num_genes1);
	free(sorted);
}

void remove_redun_candi(struct exons_list *candi_ps, int num, struct exons_list *genes, int num_genes)
{
	int i = 0, j = 0;
	struct I reg1, reg2;

	reg1 = assign_I(0, 1);
	reg2 = assign_I(0, 1);

	for( i = 0; i < num; i++ ) {
		if( candi_ps[i].sign == 'n' ) {}
		else {
			reg1 = assign_I(candi_ps[i].reg.lower-DEL_TH, candi_ps[i].reg.upper+DEL_TH);
			for( j = (i+1); j < num; j++ ) {
				reg2 = assign_I(candi_ps[j].reg.lower-DEL_TH, candi_ps[j].reg.upper+DEL_TH);
				if ( candi_ps[j].sign == 'n' ) {}
				else if( equal(candi_ps[i].reg, candi_ps[j].reg) == true ) {
					candi_ps[j].sign = 'n';
				}
				else if( subset(candi_ps[i].reg, reg2) == true ) {
					candi_ps[i].sign = 'n';
				}
				else if( subset(candi_ps[j].reg, reg1) == true ) {
					candi_ps[j].sign = 'n';
				}
			}
		}
	}

	for( i = 0; i < num; i++ ) {
		if( candi_ps[i].sign != 'n' ) {
			if( width(candi_ps[i].reg) <= (int)(0.25 * (float)(width(candi_ps[i].cmp_reg))) ) 
			{
				candi_ps[i].sign = 'n';
			}
		}
	}

	for( i = 0; i < num; i++ ) {
		if( candi_ps[i].sign == 'n' ) {}
		else {
			reg1 = assign_I(candi_ps[i].reg.lower-DEL_TH, candi_ps[i].reg.upper+DEL_TH);
			for( j = 0; j < num_genes; j++ ) {
				reg2 = assign_I(genes[j].reg.lower-DEL_TH, genes[j].reg.upper+DEL_TH);
				if ( candi_ps[i].sign == 'n' ) {}
				else if( equal(candi_ps[i].reg, genes[j].reg) == true ) {
					candi_ps[i].sign = 'n';
				}
				else if( subset(candi_ps[i].reg, reg2) == true ) {
					candi_ps[i].sign = 'n';
				}
				else if( subset(genes[j].reg, reg1) == true ) {
					candi_ps[i].sign = 'n';
				}
				else if( (proper_overlap(genes[j].reg, candi_ps[i].reg) == true) && (width(intersect(genes[j].reg, candi_ps[i].reg)) >= (int)(0.25 * (width(genes[j].reg))) ) ) {
					candi_ps[i].sign = 'n';
				}
			}
		}
	}

	for( j = 0; j < num_genes; j++ ) {
		if( (strstr(genes[j].name, "_ps") != NULL) && (genes[j].sign != 'n') ) {
			reg1 = assign_I(genes[j].reg.lower-DEL_TH, genes[j].reg.upper+DEL_TH);
			for( i = (j+1); i < num_genes; i++ ) {
				reg2 = assign_I(genes[i].reg.lower-DEL_TH, genes[i].reg.upper+DEL_TH);
				if ( genes[i].sign == 'n' ) {}
				else if( equal(genes[j].reg, genes[i].reg) == true ) {
					genes[j].sign = 'n';
				}
				else if( subset(genes[j].reg, reg2) == true ) {
					genes[j].sign = 'n';
				}
				else if( subset(genes[i].reg, reg1) == true ) {
					genes[j].sign = 'n';
				}
				else if( (proper_overlap(genes[i].reg, genes[j].reg) == true) && (width(intersect(genes[i].reg, genes[j].reg)) >= (int)(0.25 * (width(genes[i].reg))) ) ) {
					genes[j].sign = 'n';
				}
			}
		}
	}
}

void update_ps_candidates(struct exons_list cur_gene, struct exons_list * exons2, int from, int to, struct exons_list *genes1, int *num_genes1, struct DotList *algns, int num_algns, FILE *f, struct exons_list **candi_ps, int *num_candi, int *num_alloc, struct slist *sorted)
{
	int i = 0, j = 0;
	int s_loc = 0, e_loc = 0;
	struct I temp_reg;
	int *res_b, *res_e;
	int temp_num = 0;
	struct I cur_reg;
	char cur_sign = '<';

	cur_reg = assign_I(cur_gene.reg.lower, cur_gene.reg.upper);

	temp_num = *num_candi;

	res_b = (int *) ckalloc(sizeof(int));
	res_e = (int *) ckalloc(sizeof(int));

	*res_b = -1;
	*res_e = -1;
	temp_reg = assign_I(0,1);
	s_loc = search_range_b(sorted, algns, num_algns, cur_reg.lower, SELF1);
	e_loc = search_range_e(sorted, algns, num_algns, cur_reg.upper, SELF1);

	for( i = s_loc; i <= e_loc; i++ ) {
		j = sorted[i].id;
		temp_reg = assign_I(0,1);
		*res_b = -1;
		*res_e = -1;
		if(proper_overlap(algns[j].x, cur_reg) == true) {
			if( algns[j].sign == 0 ) {
				cur_sign = cur_gene.sign;
			}
			else { // algns[j].sign == 1
				if( cur_gene.sign == '<' ) {
					cur_sign = '>';
				}
				else {
					cur_sign = '<';
				}
			}

			temp_reg = intersect(algns[j].x, cur_reg);
			if( is_exons_included(temp_reg, cur_sign, exons2, from, to) == true ) {
				find_overlapping_ends(cur_reg, algns[j].x, SELF1, algns, j, f, res_b, res_e);
				if(((*res_b) != -1) && ((*res_e) != -1) && ((*res_b) < (*res_e) )) { 
					temp_reg = assign_I(*res_b, *res_e);	

					if( is_genes_included(temp_reg, cur_sign, genes1, 0, (*num_genes1)-1) == true ) 
					{
					}
					else if( is_ps_overlapped(temp_reg, cur_reg, cur_sign,  *candi_ps, temp_num) == true ) 
					{
					}
					else {
						if( temp_num >= (*num_alloc) ) {
							*candi_ps = ckrealloc(*candi_ps, ((*num_alloc) + ALLOC_UNIT) * sizeof(struct exons_list)); 	
							initialize_exons_list(*candi_ps, *num_alloc, *num_alloc + ALLOC_UNIT);
							*num_alloc = *num_alloc + ALLOC_UNIT;
						}

						(*candi_ps)[temp_num].reg = assign_I(temp_reg.lower, temp_reg.upper);
						(*candi_ps)[temp_num].cmp_reg = assign_I(cur_reg.lower, cur_reg.upper);
						(*candi_ps)[temp_num].sign = cur_sign;
						(*candi_ps)[temp_num].ctg_id = algns[j].ctg_id2;
						(*candi_ps)[temp_num].sp_id = SELF1;
						(*candi_ps)[temp_num].val = 0;
						(*candi_ps)[temp_num].fid = -1;
						strcpy((*candi_ps)[temp_num].name, cur_gene.name);
						temp_num++;
					}
				}
			}
		}
	}

	*num_candi = temp_num;
	free(res_b);
	free(res_e);
}

bool is_ps_overlapped(struct I temp_reg, struct I cmp_reg, char cur_sign, struct exons_list *candi_ps, int num)
{
	int i = 0;
	int j = -1;
	int max_width = 0;
	int temp_width = 0;
	int lo = 0, hi = 0;
	bool res = false;
	struct I reg;
	bool is_end = false;

	reg = assign_I(0,1);

	for( i = 0; i < num; i++ ) {
		if( (cur_sign == candi_ps[i].sign) && (proper_overlap(temp_reg, candi_ps[i].reg) == true) ) {
			candi_ps[i].cmp_reg = assign_I(cmp_reg.lower, cmp_reg.upper);
			temp_width = width(intersect(temp_reg, candi_ps[i].reg));
			if( temp_width > max_width ) {
				j = i;
				max_width = temp_width;
			}
		}
	}

	i = 0;
	while( (i < num) && (is_end == false) ) {
		if( (cur_sign == candi_ps[i].sign ) && (equal(cmp_reg, candi_ps[i].cmp_reg) == true )) {
			if( temp_reg.lower < candi_ps[i].reg.lower ) lo = temp_reg.lower;
			else lo = candi_ps[i].reg.lower;

			if( temp_reg.upper < candi_ps[i].reg.upper ) hi = candi_ps[i].reg.upper;
			else hi = temp_reg.upper;

			reg = assign_I(lo, hi);

			if( (width(reg) >= (int)(0.9*(float)width(cmp_reg))) && (width(reg) <= (int)(1.1*(float)width(cmp_reg))) ) {
				is_end = true;
				j = i;
			}
		}
		i++;
	}

	if( j != -1 ) {
		if( candi_ps[j].reg.lower < temp_reg.lower ) {
			lo = candi_ps[j].reg.lower;
		}
		else lo = temp_reg.lower;

		if( candi_ps[j].reg.upper > temp_reg.upper ) {
			hi = candi_ps[j].reg.upper;
		}
		else hi = temp_reg.upper;

		candi_ps[j].reg = assign_I(lo, hi);
		res = true;
	}

	return(res);
}

bool is_exons_included(struct I temp_reg, char cur_sign, struct exons_list *exons, int from, int to)
{
	int i = 0;
	bool res = false;

	for( i = from; i <= to; i++ ) {
		if( (cur_sign == exons[i].sign) && (proper_overlap(exons[i].reg, temp_reg) == true) && (width(intersect(exons[i].reg, temp_reg)) >= (int)(0.25*(float)(width(exons[i].reg))))) {
			res = true;
		}
	}

	return(res);
}

bool is_genes_included(struct I temp_reg, char cur_sign, struct exons_list *genes, int from, int to)
{
	int i = 0;
	bool res = false;

	for( i = from; i <= to; i++ ) {
		if( (cur_sign == genes[i].sign) && ((strict_subset(genes[i].reg, temp_reg) == true) || (strict_subset(temp_reg, genes[i].reg) == true))) {
			if( strstr(genes[i].name, "_ps") != NULL ) {
				if( strict_subset(genes[i].reg, temp_reg) == true ) {
					genes[i].reg = assign_I(temp_reg.lower, temp_reg.upper);
				}
			}
			res = true;
		}
	}
	return(res);
}
