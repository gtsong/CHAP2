/* GC_history - reconstruct the evolutionary history of gene conversion. */
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "maf.h"
#include "multi_util.h"
#include "contigs_op.h"

#define OVERLAP_RATIO 0.75

struct TreeNode {
	struct TreeNode *left, *right, *parent;
	char speciesname[100];
	int conversion_num, index;
};

struct GC_Array {
	char chr1[100], chr2[100], orient, identity[100], pvalue[100], outgroup[100], outgroup_name[2][100], outgroup_orient[2], orthologs_block[2][100];
	int index, b1, e1, b2, e2, length, GC_len, min_beg1, min_end1, min_beg2, min_end2, direction, outgroup_start[2], outgroup_end[2], event_index, triplet_status;
	struct GC_Array *next, *next1, *equivalent_event;
};

struct Orthologs {
	int beg, end;
	struct Orthologs *next;
};

struct TreeNode *stree;
struct TreeNode *species_nodes[100];
struct GC_Array *GC_array[100], *GC_array1[100];
struct GC_Array *GC_array_tail;
struct Orthologs *orthologs[2];
char ancestors[100][100][100][100];
int ancestors_num[100][100];
int ancestors_GC_index[100][100];
int species_num = 0;
char species_name[100][100];
int orthologs_beg[2][1000], orthologs_end[2][1000], orthologs_num[2];
int paralogs_overlap_index;
int H1_beg, H1_end, H2_beg, H2_end, H3_beg, H3_end;
float identities[3];	// 0: H1-H2;	1: H1-H3;	2: H2-H3
int global_index = 0;

char* gettoken(char *s, int *n) {
	char *t;

    t = s + *n;
    if(isalnum(s[*n])) {
		while(isalnum(s[*n]) || s[*n] == '_' || s[*n] == '-') 
			(*n)++;
		if(s[*n] == ':')
			*n += 8;
		return t;
    }
    if((s[*n] == '(')  || (s[*n] == ')' ) || (s[*n] == ',')) {
		(*n)++;
		if(s[*n] == ':')
			(*n) += 8;
		return t;
    }
	else {
		exit(-1);
	}
}

struct TreeNode *parseTree(char *s, int *n) {
    char *t;
	struct TreeNode *a, *b, *c;
	int i;
	
	t = gettoken(s, n);
	c = (void *)malloc(sizeof(struct TreeNode));
	c->parent = NULL;
	c->conversion_num = 0;

    if (t[0]=='(') {
        a = parseTree(s,n);
        t = gettoken(s,n);
        b = parseTree(s,n);
        t = gettoken(s,n);

		a->parent = c;
		b->parent = c;
		c->left = a;
		c->right = b;
		
        return c;
    }
    else {
		c->left = c->right = NULL;
		
		for(i=0; i<*n - (t - s); i++)
			c->speciesname[i] = t[i];
		c->speciesname[i] = '\0';
			
		return c;
	}
}


void extract_species_tree(char filename[100]) {
	FILE *s_file;
	char str[1000];
	char buffer[1000];
	int index, newtree;

	s_file = fopen(filename, "r");
	if(s_file == NULL) {
		printf("can't open species_tree\n");
		exit(-1);
	}
	index = 0;
	newtree = 1;
	while(1) {
		if(fscanf(s_file, "%s", buffer) == EOF) {
			index = 2;
			break;
		}
		if(buffer[strlen(buffer) - 1] == ';') { 
			buffer[strlen(buffer) - 1] = '\0';
			index = 1;
		}
		if(newtree == 1) {
			strcpy(str, buffer);
			newtree = 0;
		}
		else
			strcat(str, buffer);
		if(index == 1)
			break;
	}
		
	index = 0;
	stree = parseTree(str, &index);
	//traceTree(stree);

	fclose(s_file);

	return;
}

void trace_subtree(struct TreeNode *tree, int species, int ancestor) {
	if(tree->left == NULL) {
		strcpy(ancestors[species][ancestor][ancestors_num[species][ancestor]], tree->speciesname);
		ancestors_num[species][ancestor]++;
	}
	else {
		trace_subtree(tree->left, species, ancestor);
		trace_subtree(tree->right, species, ancestor);
	}

	return;
}

void extract_ancestors(struct TreeNode *tree, int species) {
	int ancestor = 0;

	while(tree->parent != NULL) {
		//printf("%d\n", ancestor);
		ancestors_num[species][ancestor] = 0;
		if(tree == tree->parent->left)
			trace_subtree(tree->parent->right, species, ancestor);
		else
			trace_subtree(tree->parent->left, species, ancestor);
		ancestor++;
		tree = tree->parent;
	}
	ancestors_num[species][ancestor] = 0;

	return;
}

void extract_all_ancestors(struct TreeNode *tree) {
	if(tree->left == NULL) {
		//printf("%s\n", tree->speciesname);
		extract_ancestors(tree, species_num);
		strcpy(species_name[species_num], tree->speciesname);
		species_nodes[species_num] = tree;
		species_num++;
	}
	else {
		extract_all_ancestors(tree->left);
		extract_all_ancestors(tree->right);
	}

	return;
}

int find_species(char species[100]) {
	int i;
	
	for(i=0; i<species_num; i++)
		if(strcmp(species_name[i], species) == 0)
			return i;

	printf("Can't find species : %s\n", species);
	exit(-1);
}

int find_ancestor(int species_index, char outgroup[100]) {
	int i, j; 
	
	for(i=0; ancestors_num[species_index][i]>0; i++)
		for(j=0; j<ancestors_num[species_index][i]; j++)
			if(strcmp(ancestors[species_index][i][j], outgroup) == 0)
				return i;
				
	printf("Can't find outgroup : %s\n", outgroup);
	exit(-1);
}

void GC_array_add(int species_index, int index, char chr1[100], int b1, int e1, char chr2[100], int b2, int e2, char orient, int length, char identity[100], int GC_len, char pvalue[100], int min_beg1, int min_end1, int min_beg2, int min_end2, int direction, char outgroup[100], char outgroup_name1[100], int outgroup_start1, int outgroup_end1, char outgroup_orient1, char outgroup_name2[100], int outgroup_start2, int outgroup_end2, char outgroup_orient2, char orthologs_block1[100], char orthologs_block2[100], int triplet_status) {
	if(GC_array_tail == NULL) {
		GC_array[species_index] = (void *)malloc(sizeof(struct GC_Array));
		GC_array_tail = GC_array1[species_index] = GC_array[species_index];
	}
	else {
		GC_array_tail->next = GC_array_tail->next1 = (void *)malloc(sizeof(struct GC_Array));
		GC_array_tail = GC_array_tail->next;
	}
	GC_array_tail->index = index;
	strcpy(GC_array_tail->chr1, chr1);
	GC_array_tail->b1 = b1;
	GC_array_tail->e1 = e1;
	strcpy(GC_array_tail->chr2, chr2);
	GC_array_tail->b2 = b2;
	GC_array_tail->e2 = e2;
	GC_array_tail->orient = orient;
	GC_array_tail->length = length;
	strcpy(GC_array_tail->identity, identity);
	GC_array_tail->GC_len = GC_len;
	strcpy(GC_array_tail->pvalue, pvalue);
	GC_array_tail->min_beg1 = min_beg1;
	GC_array_tail->min_end1 = min_end1;
	GC_array_tail->min_beg2 = min_beg2;
	GC_array_tail->min_end2 = min_end2;
	GC_array_tail->direction = direction;
	strcpy(GC_array_tail->outgroup, outgroup);
	strcpy(GC_array_tail->outgroup_name[0], outgroup_name1);
	GC_array_tail->outgroup_start[0] = outgroup_start1;
	GC_array_tail->outgroup_end[0] = outgroup_end1;
	GC_array_tail->outgroup_orient[0] = outgroup_orient1;
	strcpy(GC_array_tail->outgroup_name[1], outgroup_name2);
	GC_array_tail->outgroup_start[1] = outgroup_start2;
	GC_array_tail->outgroup_end[1] = outgroup_end2;
	GC_array_tail->outgroup_orient[1] = outgroup_orient2;
	strcpy(GC_array_tail->orthologs_block[0], orthologs_block1);
	strcpy(GC_array_tail->orthologs_block[1], orthologs_block2);
	GC_array_tail->event_index = -1;
	GC_array_tail->next = NULL;
	GC_array_tail->next1 = NULL;
	GC_array_tail->equivalent_event = NULL;
	GC_array_tail->triplet_status = triplet_status;
	
	return;
}

void remove_multiple_conversion_redundancy(char filename[100]) {
	FILE *fp;
	int index, index_next, b1, e1, b2, e2, length, GC_len[100], min_beg1[100], min_end1[100], min_beg2[100], min_end2[100], direction[100];
	char chr1[100], chr1_next[100], chr2[100], orient, identity[100], pvalue[100][100], outgroup[100][100], buf[1000], outgroup_name[2][100], outgroup_orient[2], orthologs_block[2][100];
	int end_of_file = 0, GC_num, i, species_index, old_species_index = -1, ancestor_index, later_ancestors_GC, remain_ancestors_GC, overlap_beg, overlap_end, outgroup_start[2], outgroup_end[2], triplet_status;
	int remain_GC[100];
	struct GC_Array *GC_array_index_start, *GC_array_equivalent_ptr[100];
	char species[100], contigs[100];

	strcpy(species, "");
	strcpy(contigs, "");
	
	fp = fopen(filename, "r");
	if(!fgets(buf, 5000, fp))
		return;

	while(!end_of_file) {
		GC_num = 0;
		GC_array_index_start = GC_array_tail;
		while(1) {
			sscanf(buf, "%d %s %d %d %s %d %d %c %d %s %d %s %d %d %d %d %d %s %d %d %c %s %d %d %c %s %s %d", &index, chr1, &b1, &e1, chr2, &b2, &e2, &orient, &length, identity, &GC_len[GC_num], pvalue[GC_num], &min_beg1[GC_num], &min_end1[GC_num], &min_beg2[GC_num], &min_end2[GC_num], &direction[GC_num], outgroup_name[0], &outgroup_start[0], &outgroup_end[0], &outgroup_orient[0], outgroup_name[1], &outgroup_start[1], &outgroup_end[1], &outgroup_orient[1], orthologs_block[0], orthologs_block[1], &triplet_status);
			concat_ctg_name(chr1, species, contigs);
			species_index = find_species(species);
			if(species_index != old_species_index) { // new species
				GC_array[species_index] = GC_array_tail = GC_array_index_start = NULL;
				old_species_index = species_index;
			}
			if(strcmp(outgroup_name[0], "NAN") == 0) {
				concat_ctg_name(outgroup_name[1], species, contigs);
				strcpy(outgroup[GC_num], species);
			}
			else {
				concat_ctg_name(outgroup_name[0], species, contigs);
				strcpy(outgroup[GC_num], outgroup_name[0]);
			}

			GC_array_add(species_index, index, chr1, b1, e1, chr2, b2, e2, orient, length, identity, GC_len[GC_num], pvalue[GC_num], min_beg1[GC_num], min_end1[GC_num], min_beg2[GC_num], min_end2[GC_num], direction[GC_num], outgroup[GC_num], outgroup_name[0], outgroup_start[0], outgroup_end[0], outgroup_orient[0], outgroup_name[1], outgroup_start[1], outgroup_end[1], outgroup_orient[1], orthologs_block[0], orthologs_block[1], triplet_status);
			GC_array_equivalent_ptr[GC_num] = GC_array_tail;
			GC_num++;
			if(!fgets(buf, 5000, fp)) {
				end_of_file = 1;
				break;
			}
			sscanf(buf, "%d %s", &index_next, chr1_next);
			if(index_next != index || strcmp(chr1_next, chr1) != 0)
				break;
		}
		
		
		for(i=0; ancestors_num[species_index][i]>0; i++)
			ancestors_GC_index[species_index][i] = -1;		
			
		// find ancestors which detect conversion
		for(i=0; i<GC_num; i++) {
			remain_GC[i] = 0;
			ancestor_index = find_ancestor(species_index, outgroup[i]);
			if(ancestors_GC_index[species_index][ancestor_index] < 0) {
				ancestors_GC_index[species_index][ancestor_index] = i;
			}
			else {
				if(atof(pvalue[i]) < atof(pvalue[ancestors_GC_index[species_index][ancestor_index]])) {
					ancestors_GC_index[species_index][ancestor_index] = i;
				}
			}
		}
		
		// determine when the conversion occurred
		later_ancestors_GC = remain_ancestors_GC = -1;
		for(i=0; ancestors_num[species_index][i]>0; i++) {
			if(ancestors_GC_index[species_index][i] >= 0) {
				if(later_ancestors_GC < 0) {
					remain_GC[ancestors_GC_index[species_index][i]] = -1;
					remain_ancestors_GC = ancestors_GC_index[species_index][i];
				}
				else {
					//check ratio of overlapped regions
					if(min_beg1[ancestors_GC_index[species_index][i]] > min_beg1[later_ancestors_GC])
						overlap_beg = min_beg1[ancestors_GC_index[species_index][i]];
					else
						overlap_beg = min_beg1[later_ancestors_GC];
					if(min_end1[ancestors_GC_index[species_index][i]] < min_end1[later_ancestors_GC])
						overlap_end = min_end1[ancestors_GC_index[species_index][i]];
					else
						overlap_end = min_end1[later_ancestors_GC];
					if((float)(overlap_end - overlap_beg + 1) / (float)(min_end1[ancestors_GC_index[species_index][i]] - min_beg1[ancestors_GC_index[species_index][i]] + 1) < OVERLAP_RATIO) {
						remain_GC[ancestors_GC_index[species_index][i]] = -1;
						remain_ancestors_GC = ancestors_GC_index[species_index][i];
					}
					else {					
						remain_GC[ancestors_GC_index[species_index][i]] = remain_ancestors_GC;
					}	
				}
				later_ancestors_GC = ancestors_GC_index[species_index][i];
			}	
		}

		for(i=0; i<GC_num; i++) {
			if(remain_GC[i] == -1) {
				GC_array_equivalent_ptr[i]->equivalent_event = NULL;
			}
			else if(remain_GC[i] == 0) {
				ancestor_index = find_ancestor(species_index, outgroup[i]);
				if(remain_GC[ancestors_GC_index[species_index][ancestor_index]] == -1)
					GC_array_equivalent_ptr[i]->equivalent_event = GC_array_equivalent_ptr[ancestors_GC_index[species_index][ancestor_index]];
				else
					GC_array_equivalent_ptr[i]->equivalent_event = GC_array_equivalent_ptr[remain_GC[ancestors_GC_index[species_index][ancestor_index]]];
			}
			else {
				GC_array_equivalent_ptr[i]->equivalent_event = GC_array_equivalent_ptr[remain_GC[i]];
			}
		}
		for(i=0; i<GC_num; i++) {
			if(remain_GC[i] >= 0) {
				if(GC_array_index_start == NULL) {
					GC_array[species_index] = GC_array[species_index]->next;
				}
				else {
					GC_array_index_start->next = GC_array_index_start->next->next;
					if(GC_array_index_start->next == NULL)
						GC_array_tail = GC_array_index_start;
				}
			}
			else {
				if(GC_array_index_start == NULL) {
					GC_array_index_start = GC_array[species_index];
				}
				else {
					GC_array_index_start = GC_array_index_start->next;
				}
			}
		}
	}
	
	return;
}

void extract_orthologs(char species1[100], char species2[100], int beg, int end, int index) {
	struct mafFile *mf;
	struct mafAli *ali;
	struct mafComp *mc1, *mc2;
	char filename[1000];
	struct Orthologs *orthologs_tail;
	int overlap_beg, overlap_end, col, i;
	
	orthologs[index] = orthologs_tail = NULL;
	sprintf(filename, "ortho.d/many_to_many_ortho.d/%s.%s.maf", species1, species2);
	mf = mafOpen(filename, 0);
	while((ali = mafNext(mf)) != NULL) {
		mc1 = ali->components;
		if(mc1 == NULL) {
			printf("Wrong orthologs file\n");
			exit(1);
		}
		//check overlapping
		if((beg <= mc1->start + mc1->size - 1) && (end >= mc1->start)) {
			//determine overlapped region
			if(beg >= mc1->start)
				overlap_beg = beg;
			else
				overlap_beg = mc1->start;
			if(end <= mc1->start + mc1->size - 1)
				overlap_end = end;
			else
				overlap_end = mc1->start + mc1->size - 1;
				
			//find orthologs
			mc2 = mc1->next;
			if(mc2 == NULL) {
				printf("Wrong orthologs file\n");
				exit(1);
			}
			
			//memory allocation
			if(orthologs_tail == NULL) {
				orthologs[index] = (void *)malloc(sizeof(struct Orthologs));
				orthologs_tail = orthologs[index];
			}
			else {
				orthologs_tail->next = (void *)malloc(sizeof(struct Orthologs));
				orthologs_tail = orthologs_tail->next;
			}
			
			orthologs_tail->next = NULL;
			
			col = mafPos2Col_v2(mc1, overlap_beg);
			if(mc2->strand == '+') {
				for(orthologs_tail->beg=mc2->start - 1, i=0; i<=col; i++)
					if(mc2->text[i] != '-')
						orthologs_tail->beg++;
			}
			else {
				for(orthologs_tail->end=mc2->srcSize - mc2->start, i=0; i<=col; i++)
					if(mc2->text[i] != '-')
						orthologs_tail->end--;
			}
			
			col = mafPos2Col_v2(mc1, overlap_end);
			if(mc2->strand == '+') {
				for(orthologs_tail->end=mc2->start - 1, i=0; i<=col; i++)
					if(mc2->text[i] != '-')
						orthologs_tail->end++;
			}
			else {
				for(orthologs_tail->beg=mc2->srcSize - mc2->start, i=0; i<=col; i++)
					if(mc2->text[i] != '-')
						orthologs_tail->beg--;
			}
		}
	}

	return;
}

int check_orthologs(int beg, int end, int index) {
	struct Orthologs *orthologs_tail;
	int overlap_beg, overlap_end, union_beg, union_end;
	
	orthologs_tail = orthologs[index];
	while(orthologs_tail != NULL) {
		if(beg >= orthologs_tail->beg) {
			overlap_beg = beg;
			union_beg = orthologs_tail->beg;
		}
		else {
			overlap_beg = orthologs_tail->beg;
			union_beg = beg;
		}
		if(end <= orthologs_tail->end) {
			overlap_end = end;
			union_end = orthologs_tail->end;
		}
		else {
			overlap_end = orthologs_tail->end;
			union_end = end;
		}
		if((float)(overlap_end - overlap_beg + 1) / (float)(union_end - union_beg + 1) > OVERLAP_RATIO)
			return 1;
		orthologs_tail = orthologs_tail->next;
	}
	
	return 0;
}

int check_outgroup(int species_index, char outgroup1[100], char outgroup2[100]) {
	int i, j;
	
	for(i=0; ancestors_num[species_index][i]>0; i++) {
		for(j=0; j<ancestors_num[species_index][i]; j++)
			if(strcmp(outgroup1, ancestors[species_index][i][j]) == 0) 
				break;
		if(j<ancestors_num[species_index][i])
			for(j=0; j<ancestors_num[species_index][i]; j++)
				if(strcmp(outgroup2, ancestors[species_index][i][j]) == 0) 
					return 1;
	}
	
	return 0;
}

void remove_orthologs_redundancy() {
	int i, j, k, l, m;
	struct GC_Array *GC_array_tail1, *GC_array_tail2, *GC_array_pre;
	
	for(i=0; i<species_num; i++) {
		GC_array_tail1 = GC_array[i];
		while(GC_array_tail1 != NULL) {
			//find ancestor's level for outgroup
			for(j=0; ancestors_num[i][j]>0; j++) {
				for(k=0; k<ancestors_num[i][j]; k++)
					if(strcmp(ancestors[i][j][k], GC_array_tail1->outgroup) == 0)
						break;
				if(k<ancestors_num[i][j])
					break;
			}
			if(ancestors_num[i][j] <= 0) {
				printf("Wrong outgroup: %s\n", GC_array_tail1->outgroup);
				exit(1);
			}
			//remove redundant conversion events in other species
			for(k=i+1; k<species_num; k++) {
				for(l=0; l<j; l++) {
					for(m=0; m<ancestors_num[i][l]; m++) {
						if(strcmp(species_name[k], ancestors[i][l][m]) == 0) {
							extract_orthologs(species_name[i], species_name[k], GC_array_tail1->min_beg1, GC_array_tail1->min_end1, 0);
							extract_orthologs(species_name[i], species_name[k], GC_array_tail1->min_beg2, GC_array_tail1->min_end2, 1);
						
							if(orthologs[0] != NULL && orthologs[1] != NULL) {
								GC_array_tail2 = GC_array[k];
								GC_array_pre = NULL;
								while(GC_array_tail2 != NULL) {					
									//check if the outgrougs are in the same level
									if(check_outgroup(i, GC_array_tail1->outgroup, GC_array_tail2->outgroup) == 1) {
										//check if two conversion regions are orthologous
										if((check_orthologs(GC_array_tail2->min_beg1, GC_array_tail2->min_end1, 0) && check_orthologs(GC_array_tail2->min_beg2, GC_array_tail2->min_end2, 1)) || (check_orthologs(GC_array_tail2->min_beg1, GC_array_tail2->min_end1, 1) && check_orthologs(GC_array_tail2->min_beg2, GC_array_tail2->min_end2, 0))) {
											//remove the second conversion event
											GC_array_tail2->equivalent_event = GC_array_tail1;
											if(GC_array_pre == NULL)
												GC_array[k] = GC_array_tail2->next;
											else
												GC_array_pre->next = GC_array_tail2->next;
											
											GC_array_tail2 = GC_array_tail2->next;
											continue;
										}
									}
								
									GC_array_pre = GC_array_tail2;
									GC_array_tail2 = GC_array_tail2->next;
								}
							}
						}
					}
					if(m<ancestors_num[i][l])
						break;
				}
			}

			GC_array_tail1 = GC_array_tail1->next;
		}
	}

	return;
}

int check_overlapping(int b1, int e1, int b2, int e2, int index) {
	int overlap_beg, overlap_end;
	
	if(e1 < b2 || e2 < b1)
		return 0;
	
	if(b1 > b2)
		overlap_beg = b1;
	else
		overlap_beg = b2;
	if(e1 < e2)
		overlap_end = e1;
	else
		overlap_end = e2;
	
	if((float)(overlap_end - overlap_beg + 1) / (float)(e1 - b1 + 1) > OVERLAP_RATIO && (float)(overlap_end - overlap_beg + 1) / (float)(e2 - b2 + 1) > OVERLAP_RATIO) {
		paralogs_overlap_index = index;
		H3_beg = overlap_beg;
		H3_end = overlap_end;
		return 1;
	}
	else
		return 0;
}

void calculate_similarities(int species) {
	struct mafFile *mf;
	struct mafAli *ali, *ali_temp[10000];
	struct mafComp *mc1, *mc2;
	char filename[1000];
	int i, col1, col2, same, different, ali_temp_num = 0, ali_temp_flag[10000];
	
	//printf("H1_beg = %d; H1_end = %d; H2_beg = %d; H2_end = %d; H3_beg = %d; H3_end = %d\n", H1_beg, H1_end, H2_beg, H2_end, H3_beg, H3_end);
	sprintf(filename, "%s.d/%s.remove_repeats.maf", species_name[species], species_name[species]);
	//sprintf(filename, "ortho.d/self.d/%s.maf", species_name[species]);
	mf = mafOpen(filename, 0);
	identities[0] = identities[1] = identities[2] = 0.0;
	while((ali = mafNext(mf)) != NULL) {
		mc1 = ali->components;
		if(mc1 == NULL) {
			printf("Wrong maf file\n");
			exit(1);
		}
		mc2 = mc1->next;
		if(mc2 == NULL) {
			printf("Wrong maf file\n");
			exit(1);
		}
		
		//ignore lower right alignments
		if((mc2->strand == '+' && mc1->start >= mc2->start) || (mc2->strand == '-' && mc1->start >= mc2->srcSize - mc2->start - mc2->size))
			continue;
		
		if(H3_end >= mc1->start && H3_beg <= mc1->start + mc1->size - 1) {	//	mc1 contains H3
			if(H3_beg < mc1->start || H3_end > mc1->start + mc1->size -1)
				continue;
			if (((mc2->strand == '+') && (H1_end >= mc2->start && H1_beg <= mc2->start + mc2->size - 1)) || ((mc2->strand == '-') && (H1_end >= mc2->srcSize - mc2->start - mc2->size && H1_beg <= mc2->srcSize - mc2->start - 1))) {	// mc2 contains H1
				col1 = mafPos2Col_v2(mc1, H3_beg);
				col2 = mafPos2Col_v2(mc1, H3_end);
				
				//calculate the similarity of H1 and H3
				same = different = 0;
				for(i=col1; i<=col2; i++)
					if(toupper(mc1->text[i]) == toupper(mc2->text[i]))
						same++;
					else
						different++;
				identities[1] = (float)same / (float)(same + different);
				//printf("H1-H3: same = %d\tdifferent = %d\n", same, different);
				/*
				//recalculate the position of H1
				if(mc2->strand == '+') {
					for(H1_beg=mc2->start - 1, i=0; i<=col1; i++)
						if(mc2->text[i] != '-')
							H1_beg++;
					for(H1_end=mc2->start - 1, i=0; i<=col2; i++)
						if(mc2->text[i] != '-')
							H1_end++;
				}
				else {
					for(H1_end=mc2->srcSize - mc2->start, i=0; i<=col1; i++)
						if(mc2->text[i] != '-')
							H1_end--;
					for(H1_beg=mc2->srcSize - mc2->start, i=0; i<=col2; i++)
						if(mc2->text[i] != '-')
							H1_beg--;
				}
				*/
			} 
			else if(((mc2->strand == '+') && (H2_end >= mc2->start && H2_beg <= mc2->start + mc2->size - 1)) || ((mc2->strand == '-') && (H2_end >= mc2->srcSize - mc2->start - mc2->size && H2_beg <= mc2->srcSize - mc2->start - 1))) {	// mc2 contains H2
				col1 = mafPos2Col_v2(mc1, H3_beg);
				col2 = mafPos2Col_v2(mc1, H3_end);
				
				//calculate the similarity of H2 and H3
				same = different = 0;
				for(i=col1; i<=col2; i++)
					if(toupper(mc1->text[i]) == toupper(mc2->text[i]))
						same++;
					else
						different++;
				identities[2] = (float)same / (float)(same + different);
				//printf("H2-H3: same = %d\tdifferent = %d\n", same, different);
				/*
				//recalculate the position of H2
				if(mc2->strand == '+') {
					for(H2_beg=mc2->start - 1, i=0; i<=col1; i++)
						if(mc2->text[i] != '-')
							H2_beg++;
					for(H2_end=mc2->start - 1, i=0; i<=col2; i++)
						if(mc2->text[i] != '-')
							H2_end++;
				}
				else {
					for(H2_end=mc2->srcSize - mc2->start, i=0; i<=col1; i++)
						if(mc2->text[i] != '-')
							H2_end--;
					for(H2_beg=mc2->srcSize - mc2->start, i=0; i<=col2; i++)
						if(mc2->text[i] != '-')
							H2_beg--;
				}
				*/
			}
		}
		else if(((mc2->strand == '+') && (H3_end >= mc2->start && H3_beg <= mc2->start + mc2->size - 1)) || ((mc2->strand == '-') && (H3_end >= mc2->srcSize - mc2->start - mc2->size && H3_beg <= mc2->srcSize - mc2->start - 1))) {	// mc2 contains H3
			if(((mc2->strand == '+') && (H3_beg < mc2->start || H3_end > mc2->start + mc2->size - 1)) || ((mc2->strand == '-') && (mc2->srcSize - H3_end - 1 < mc2->start || mc2->srcSize - H3_beg - 1 > mc2->start + mc2->size - 1)))
					continue;
			
			if(H1_end >= mc1->start && H1_beg <= mc1->start + mc1->size - 1) {	// mc1 contains H1
				if(mc2->strand == '+') {
					col1 = mafPos2Col_v2(mc2, H3_beg);
					col2 = mafPos2Col_v2(mc2, H3_end);
				}
				else {
					col1 = mafPos2Col_v2(mc2, mc2->srcSize - H3_end - 1);
					col2 = mafPos2Col_v2(mc2, mc2->srcSize - H3_beg - 1);
				}
				
				//calculate the similarity of H1 and H3
				same = different = 0;
				for(i=col1; i<=col2; i++)
					if(toupper(mc1->text[i]) == toupper(mc2->text[i]))
						same++;
					else
						different++;
				identities[1] = (float)same / (float)(same + different);
				//printf("H1-H3: same = %d\tdifferent = %d\n", same, different);
				/*
				//recalculate the position of H1
				for(H1_beg=mc1->start - 1, i=0; i<=col1; i++)
					if(mc1->text[i] != '-')
						H1_beg++;
				for(H1_end=mc1->start - 1, i=0; i<=col2; i++)
					if(mc1->text[i] != '-')
						H1_end++;
				*/
			} 
			else if(H2_end >= mc1->start && H2_beg <= mc1->start + mc1->size - 1) {	// mc1 contains H2
				if(mc2->strand == '+') {
					col1 = mafPos2Col_v2(mc2, H3_beg);
					col2 = mafPos2Col_v2(mc2, H3_end);
				}
				else {
					col1 = mafPos2Col_v2(mc2, mc2->srcSize - H3_end - 1);
					col2 = mafPos2Col_v2(mc2, mc2->srcSize - H3_beg - 1);
				}

				//calculate the similarity of H2 and H3
				same = different = 0;
				for(i=col1; i<=col2; i++)
					if(toupper(mc1->text[i]) == toupper(mc2->text[i]))
						same++;
					else
						different++;
				identities[2] = (float)same / (float)(same + different);
				//printf("H2-H3: same = %d\tdifferent = %d\n", same, different);
				/*
				//recalculate the position of H2
				for(H2_beg=mc1->start - 1, i=0; i<=col1; i++)
					if(mc1->text[i] != '-')
						H2_beg++;
				for(H2_end=mc1->start - 1, i=0; i<=col2; i++)
					if(mc1->text[i] != '-')
						H2_end++;
				*/
			} 
		}
		
		//store alignments contain both H1 and H2
		if((H1_end >= mc1->start && H1_beg <= mc1->start + mc1->size - 1) && (((mc2->strand == '+') && (H2_end >= mc2->start && H2_beg <= mc2->start + mc2->size - 1)) || ((mc2->strand == '-') && (H2_end >= mc2->srcSize - mc2->start - mc2->size && H2_beg <= mc2->srcSize - mc2->start - 1)))) {	//mc1 contains H1 and mc2 contains H2
			ali_temp[ali_temp_num] = ali;
			ali_temp_flag[ali_temp_num] = 0;
			ali_temp_num++;
		}
		else if((H2_end >= mc1->start && H2_beg <= mc1->start + mc1->size - 1) && (((mc2->strand == '+') && (H1_end >= mc2->start && H1_beg <= mc2->start + mc2->size - 1)) || ((mc2->strand == '-') && (H1_end >= mc2->srcSize - mc2->start - mc2->size && H1_beg <= mc2->srcSize - mc2->start - 1)))) {	//mc1 contains H2 and mc2 contains H1
			ali_temp[ali_temp_num] = ali;
			ali_temp_flag[ali_temp_num] = 1;
			ali_temp_num++;
		}
	}
	
	int largest_index = 0, largest_size = 0;
	//calculate similarity of H1 and H2
	if(ali_temp_num > 1) {		
		for(i=0; i<ali_temp_num; i++) {
			if(ali_temp[i]->components->size > largest_size) {
				largest_size = ali_temp[i]->components->size;
				largest_index = i;
			}
		}
	}
	else if (ali_temp_num == 0) {
		identities[0] = 0.0;
		return;
	}

	mc1 = ali_temp[largest_index]->components;
	mc2 = mc1->next;
	if(ali_temp_flag[largest_index] == 0) {
		if(H1_beg < mc1->start)
			col1 = mafPos2Col_v2(mc1, mc1->start);
		else
			col1 = mafPos2Col_v2(mc1, H1_beg);
		if(H1_end > mc1->start + mc1->size - 1)
			col2 = mafPos2Col_v2(mc1, mc1->start + mc1->size - 1);
		else
			col2 = mafPos2Col_v2(mc1, H1_end);
	}
	else {
		if(H2_beg < mc1->start)
			col1 = mafPos2Col_v2(mc1, mc1->start);
		else
			col1 = mafPos2Col_v2(mc1, H2_beg);
		if(H2_end > mc1->start + mc1->size - 1)
			col2 = mafPos2Col_v2(mc1, mc1->start + mc1->size - 1);
		else
			col2 = mafPos2Col_v2(mc1, H2_end);
	}
	
	same = different = 0;
	for(i=col1; i<=col2; i++)
		if(toupper(mc1->text[i]) == toupper(mc2->text[i]))
			same++;
		else
			different++;
	identities[0] = (float)same / (float)(same + different);
	//printf("H1-H2: same = %d\tdifferent = %d\n", same, different);
	
	return;
}

void remove_paralogs_redundancy() {
	int i;
	struct GC_Array *GC_array_tail1, *GC_array_tail2, *GC_array_pre1, *GC_array_pre2;

	for(i=0; i<species_num; i++) {
		GC_array_tail1 = GC_array[i];
		GC_array_pre1 = NULL;
		while(GC_array_tail1 != NULL) {
			GC_array_tail2 = GC_array_tail1->next;
			GC_array_pre2 = GC_array_tail1;
			while(GC_array_tail2 != NULL) {
				//check if the outgrougs are in the same level
				if(check_outgroup(i, GC_array_tail1->outgroup, GC_array_tail2->outgroup) == 1) {
					//check overlapping
					paralogs_overlap_index = -1;
					if(check_overlapping(GC_array_tail1->min_beg1, GC_array_tail1->min_end1, GC_array_tail2->min_beg1, GC_array_tail2->min_end1, 0) || check_overlapping(GC_array_tail1->min_beg1, GC_array_tail1->min_end1, GC_array_tail2->min_beg2, GC_array_tail2->min_end2, 1) || check_overlapping(GC_array_tail1->min_beg2, GC_array_tail1->min_end2, GC_array_tail2->min_beg1, GC_array_tail2->min_end1, 2) || check_overlapping(GC_array_tail1->min_beg2, GC_array_tail1->min_end2, GC_array_tail2->min_beg2, GC_array_tail2->min_end2, 3)) {
						if(paralogs_overlap_index < 0) {
							printf("Wrong paralogs_overlap_index: %d\n", paralogs_overlap_index);
							exit(1);
						}
						//calculate similarities
						if(paralogs_overlap_index == 0) {
							H1_beg = GC_array_tail1->min_beg2;
							H1_end = GC_array_tail1->min_end2;
							H2_beg = GC_array_tail2->min_beg2;
							H2_end = GC_array_tail2->min_end2;
						} else if(paralogs_overlap_index == 1) {
							H1_beg = GC_array_tail1->min_beg2;
							H1_end = GC_array_tail1->min_end2;
							H2_beg = GC_array_tail2->min_beg1;
							H2_end = GC_array_tail2->min_end1;
						} else if(paralogs_overlap_index == 2) {
							H1_beg = GC_array_tail1->min_beg1;
							H1_end = GC_array_tail1->min_end1;
							H2_beg = GC_array_tail2->min_beg2;
							H2_end = GC_array_tail2->min_end2;
						} else {
							H1_beg = GC_array_tail1->min_beg1;
							H1_end = GC_array_tail1->min_end1;
							H2_beg = GC_array_tail2->min_beg1;
							H2_end = GC_array_tail2->min_end1;
						}
						calculate_similarities(i);
						//printf("H1-H2: %f\tH1-H3: %f\tH2-H3: %f\n\n", identities[0], identities[1], identities[2]);
						//remove redundencies
						if(identities[0] > 0.0) {
							if(identities[2] <= identities[0]) {
								//remove conversion between H2-H3
								GC_array_tail2->equivalent_event = GC_array_tail1;
								GC_array_pre2->next = GC_array_tail2->next;
								GC_array_tail2 = GC_array_tail2->next;
								continue;
							}
							else if(identities[1] <= identities[0]) {
								//remove conversion between H1-H3
								GC_array_tail1->equivalent_event = GC_array_tail2;
								if(GC_array_pre1 == NULL)
									GC_array[i] = GC_array_tail1->next;
								else
									GC_array_pre1->next = GC_array_tail1->next;
						
								GC_array_tail1 = GC_array_tail1->next;					
								break;
							}
						}
					}
				}
				
				GC_array_pre2 = GC_array_tail2;
				GC_array_tail2 = GC_array_tail2->next;
			}
			
			GC_array_pre1 = GC_array_tail1;
			GC_array_tail1 = GC_array_tail1->next;
		}
	}
	
	return;
}

void show_tree(struct TreeNode *tree, FILE *out_fp) {
	if(tree->left == NULL) {
		fprintf(out_fp, "%s[%d]", tree->speciesname, tree->conversion_num);
	}
	else {
		fprintf(out_fp, "(");
		show_tree(tree->left, out_fp);
		fprintf(out_fp, ",");
		show_tree(tree->right, out_fp);
		fprintf(out_fp, ")[%d]", tree->conversion_num);
	}

	return;
}

void show_tree_with_index(struct TreeNode *tree, FILE *out_fp) {
	if(tree->left == NULL) {
		global_index++;
		tree->index = global_index;
		fprintf(out_fp, "%s[%d]", tree->speciesname, global_index);
	}
	else {
		fprintf(out_fp, "(");
		show_tree_with_index(tree->left, out_fp);
		fprintf(out_fp, ",");
		show_tree_with_index(tree->right, out_fp);
		global_index++;
		tree->index = global_index;
		fprintf(out_fp, ")[%d]", global_index);
	}

	return;
}

void show_converion_in_tree() {
	int i, j, outgroup_index;
	struct TreeNode *species_node, *conversion_node;
	FILE *out_fp;

	for(i=0; i<species_num; i++) {
		species_node = species_nodes[i];
		GC_array_tail = GC_array[i];
		while(GC_array_tail != NULL) {
			outgroup_index = find_ancestor(i, GC_array_tail->outgroup);
			conversion_node = species_node;
			for(j=0; j<outgroup_index; j++)
				conversion_node = conversion_node->parent;
			conversion_node->conversion_num++;
			GC_array_tail = GC_array_tail->next;
		}
	}
	
	out_fp = fopen("conversion_in_tree.txt", "w");
	show_tree(stree, out_fp);
	fclose(out_fp);
	
	out_fp = fopen("species_tree_with_index.txt", "w");
	show_tree_with_index(stree, out_fp);
	fprintf(out_fp, "\n");
	fclose(out_fp);
}

int find_tree_index(char species[100], char outgroup[100]) {
	int i, species_index, ancestor_index;
	struct TreeNode *species_node;
	
	species_index = find_species(species);
	species_node = species_nodes[species_index];
	ancestor_index = find_ancestor(species_index, outgroup);
	for(i=0; i<ancestor_index; i++) {
		if(species_node->parent == NULL) {
			printf("Wrong outgroup: %s\n", outgroup);
			exit(1);
		}
		else {
			species_node = species_node->parent;
		}
	}
	
	return species_node->index;
}

int main(int argc, char *argv[]) { 
	int i, event_num = 0;
	FILE *fp;
	struct GC_Array *GC_array_equivalent_ptr;

	if ( argc != 3) {
		printf("GC_history species_tree conversion_file\n");
		return 1;
	}
	extract_species_tree(argv[1]);
	extract_all_ancestors(stree);
	remove_multiple_conversion_redundancy(argv[2]);
	remove_orthologs_redundancy();
	remove_paralogs_redundancy();
	show_converion_in_tree();
	for(i=0; i<species_num; i++) {
		GC_array_tail = GC_array[i];
		while(GC_array_tail != NULL) {
			event_num++;
			GC_array_tail->event_index = event_num;
			printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%3.2f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n", GC_array_tail->event_index, GC_array_tail->chr1, GC_array_tail->b1, GC_array_tail->e1, GC_array_tail->chr2, GC_array_tail->b2, GC_array_tail->e2, GC_array_tail->orient, GC_array_tail->length, atof(GC_array_tail->identity), GC_array_tail->GC_len, GC_array_tail->pvalue, GC_array_tail->min_beg1, GC_array_tail->min_end1, GC_array_tail->min_beg2, GC_array_tail->min_end2, GC_array_tail->direction, GC_array_tail->outgroup);			
			GC_array_tail = GC_array_tail->next;
		}
	}
	
	fp = fopen("all_1.gc", "w");
	for(i=0; i<species_num; i++) {
		GC_array_tail = GC_array1[i];
		while(GC_array_tail != NULL) {
			GC_array_equivalent_ptr = GC_array_tail;
			while(GC_array_equivalent_ptr->equivalent_event != NULL)
				GC_array_equivalent_ptr = GC_array_equivalent_ptr->equivalent_event;
			fprintf(fp, "%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%3.2f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%d\t%d\t%s\t%s\t%d\n", GC_array_tail->index, GC_array_tail->chr1, GC_array_tail->b1, GC_array_tail->e1, GC_array_tail->chr2, GC_array_tail->b2, GC_array_tail->e2, GC_array_tail->orient, GC_array_tail->length, atof(GC_array_tail->identity), GC_array_tail->GC_len, GC_array_tail->pvalue, GC_array_tail->min_beg1, GC_array_tail->min_end1, GC_array_tail->min_beg2, GC_array_tail->min_end2, GC_array_tail->direction, GC_array_tail->outgroup_name[0], GC_array_tail->outgroup_start[0], GC_array_tail->outgroup_end[0], GC_array_tail->outgroup_orient[0], GC_array_tail->outgroup_name[1], GC_array_tail->outgroup_start[1], GC_array_tail->outgroup_end[1], GC_array_tail->outgroup_orient[1], GC_array_equivalent_ptr->event_index, find_tree_index(GC_array_equivalent_ptr->chr1, GC_array_equivalent_ptr->outgroup), GC_array_tail->orthologs_block[0], GC_array_tail->orthologs_block[1], GC_array_tail->triplet_status);			
			GC_array_tail = GC_array_tail->next1;
		}
	}
	fclose(fp);

	return 0;
}
