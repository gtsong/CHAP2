#include "main.h"
#include "analysis.h"
#include "read_algn.h"
#include "util.h"
#include "regions.h"
#include "util_i.h"
#include "util_gen.h"

//void get_numbers(struct DotList *content_ortho, int num_content, struct DotList *position_ortho, int num_position, struct exons_list *exons, int num_exons, FILE *f, FILE *g, int size1, int size2)
void get_numbers(struct DotList *content_ortho, int num_content, struct DotList *position_ortho, int num_position, FILE *f, FILE *g, int size1, int size2)
{
	struct slist *mark1, *mark2, *common;
	int *num1, *num2;
	int num_algns_content = 0, num_algns_position = 0, num_common_algns = 0;
	int num_bases_content = 0, num_bases_position = 0, num_common_bases = 0;
	int num_mbases_content = 0, num_mbases_position = 0, num_common_mbases = 0;
	int i = 0, j = 0;

	mark1 = (struct slist *) ckalloc(num_content * sizeof(struct slist));
	common = (struct slist *) ckalloc(num_content * sizeof(struct slist));
	mark2 = (struct slist *) ckalloc(num_position * sizeof(struct slist));
	num1 = (int *) ckalloc(sizeof(int));
	num2 = (int *) ckalloc(sizeof(int));

	initialize_slist(mark1, 0, num_content);
	initialize_slist(mark2, 0, num_position);
	for( i = 0; i < num_content; i++ ) {
		mark1[i].val = 0; // the number of nucleotides in a local alignment
		mark1[i].val_red = 0; // the number of nucleotides in common
		mark1[i].sp_state = 0; // the number of aligned nucleotides 
		mark1[i].add_sp_state = 0; // the number of matched nucleotides
		mark1[i].id = content_ortho[i].indiv_fid;
		mark1[i].is_x = false; // is this alignment common in both of ortho_content and ortho_position
	}
	
	for( i = 0; i < num_position; i++ ) {
		mark2[i].val = 0;
		mark2[i].val_red = 0;
		mark2[i].sp_state = 0;
		mark2[i].add_sp_state = 0;
		mark2[i].id = position_ortho[i].indiv_fid;
		mark2[i].is_x = false;
	}

	for( i = 0; i < num_content; i++ ) {
		common[i].val = 0;
		common[i].val_red = 0;
		common[i].sp_state = 0;
		common[i].add_sp_state = 0;
		common[i].id = content_ortho[i].indiv_fid;
		common[i].is_x = false;
	}

	for( i = 0; i < num_content; i++ ) {
		*num1 = 0;
		*num2 = 0;
		mark1[i].add_sp_state = count_nucs(content_ortho[i], f, width(content_ortho[i].x), num1, num2);
		mark1[i].val = *num1;
		mark1[i].sp_state = *num2;
	}

	for( i = 0; i < num_position; i++ ) {
		*num1 = 0;
		*num2 = 0;
		mark2[i].add_sp_state = count_nucs(position_ortho[i], g, width(position_ortho[i].x), num1, num2);
		mark2[i].val = *num1;
		mark2[i].sp_state = *num2;
	}

	compute_common_algns(content_ortho, num_content, position_ortho, num_position, f, common);

	for( i = 0; i < num_content; i++ ) {
		if( mark1[i].is_x == false ) {
			for( j = (i+1); j < num_content; j++ ) {
				if( mark1[i].id == mark1[j].id ) {
					mark1[j].is_x = true;
				}
			}
			num_algns_content++; // the number of alignments
		}
		num_bases_content = num_bases_content + mark1[i].sp_state;
		num_mbases_content = num_mbases_content + mark1[i].add_sp_state;
	}

	for( i = 0; i < num_position; i++ ) {
		if( mark2[i].is_x == false ) {
			for( j = (i+1); j < num_position; j++ ) {
				if( mark2[i].id == mark2[j].id ) {
					mark2[j].is_x = true;
				}
			}
			num_algns_position++;
		}
		num_bases_position = num_bases_position + mark2[i].sp_state;
		num_mbases_position = num_mbases_position + mark2[i].add_sp_state;
	}

	for( i = 0; i < num_content; i++ ) {
		if( common[i].is_x == false ) {
			for( j = (i+1); j < num_content; j++ ) {
				if( common[i].id == common[j].id ) {
					common[j].is_x = true;
					common[i].sp_state = common[i].sp_state + common[j].sp_state; // the number of common bases
					common[i].add_sp_state = common[i].add_sp_state + common[j].sp_state; // the number of common bases
				}
			}
			if( common[i].val >= ERR_SM_TH ) num_common_algns++; // the number of common alignments
			num_common_bases = num_common_bases + common[i].sp_state;
			num_common_mbases = num_common_mbases + common[i].add_sp_state;
		}
	}

	printf("%d %d %d %d %d %d %d %d %d %d %d ", num_algns_content, num_algns_position, num_common_algns, num_bases_content, num_bases_position, num_common_bases, num_mbases_content, num_mbases_position, num_common_mbases, size1, size2);

	free(common);
	free(mark1);
	free(mark2);
	free(num1);
	free(num2);
}

void compute_common_algns(struct DotList *content_ortho, int num_content, struct DotList *position_ortho, int num_position, FILE *f, struct slist *common)
{
	int i = 0, j = 0;
	struct I cmp, tmp;
	int *num1, *num2;
	int temp_num1 = 0, temp_num2 = 0, temp_num3 = 0;
	int b = 0, e = 0;

	num1 = (int *) ckalloc(sizeof(int));
	num2 = (int *) ckalloc(sizeof(int));
	cmp = assign_I(0, 1);
	tmp = assign_I(0, 1);

	for( i = 0; i < num_content; i++ ) {
		cmp = assign_I(content_ortho[i].x.lower, content_ortho[i].x.upper);
		*num1 = 0;
		*num2 = 0;
		for( j = 0; j < num_position; j++ ) {
			tmp = assign_I(position_ortho[j].x.lower, position_ortho[j].x.upper);
			temp_num3 = 0;
			temp_num1 = 0;
			temp_num2 = 0;
			if( (content_ortho[i].indiv_fid == position_ortho[j].indiv_fid) && (proper_overlap(cmp, tmp) == true) ) {
				if( cmp.lower < tmp.lower ) {
					b = tmp.lower;
				}
				else {
					b = cmp.lower;
				}
				
				if( cmp.upper < tmp.upper ) {
					e = cmp.upper;
				}
				else {
					e = tmp.upper;
				}

				if( e > content_ortho[i].x.lower ) {
					temp_num3 = count_nucs(content_ortho[i], f, e-content_ortho[i].x.lower, num1, num2);
					temp_num1 = *num1;
					temp_num2 = *num2;
				}
				else {
					temp_num3 = 0;
					temp_num1 = 0;
					temp_num2 = 0;
				}

				if( b > content_ortho[i].x.lower ) {
					temp_num3 = temp_num3 - count_nucs(content_ortho[i], f, b-content_ortho[i].x.lower, num1, num2);
					temp_num1 = temp_num1 - (*num1);
					temp_num2 = temp_num2 - (*num2);
				}
				else {
				}

				common[i].val = common[i].val + temp_num1;
				common[i].sp_state = common[i].sp_state + temp_num2;
				common[i].add_sp_state = common[i].add_sp_state + temp_num3;
			}
		}
	}

	free(num1);
	free(num2);
}
