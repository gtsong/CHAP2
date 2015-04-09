// interval: [a,b]
#include "main.h"
#include "regions.h"
#include "util.h"
#include "util_i.h"
#include "util_gen.h"
#include "read_maf.h"
#include "tree_op.h"

int debug_mode;
char S[BIG], T[BIG];

int sort_merge_intervals(struct I *regs, int num_regs, struct I *new_regs);
int count_bp(struct I *regs, int num_regs);
int count_overlap_bp(struct I *reg1, int num_reg1, struct I *reg2, int num_reg2);
int find_status(struct cv_list cur_cv, struct cv_list *cv, int num_cv);
bool is_pseudo_gene(char *name);
void init_g_list(struct g_list *genes, int b, int e);

int main(int argc, char **argv)
{
	struct DotList *algns;
	int *num_algns;
	int *size1, *size2;
	FILE *cv_f;
	FILE *f;
 	struct cv_list *cv; // the conversion events detected by Chih-Hao's program
	int num_cv;
 	struct cv_list *init_cv; // the conversion events detected by Chih-Hao's program
	int num_init_cv;
	char buf[1000];
  char tree_line[1000], filename[1000], gene_name[100];
  struct sp_list *sp_id;
	int num_sp = 0, temp_num = 0;
	struct I src, dst;
	struct I cmp1, cmp2;
	int i = 0, j = 0, k = 0, l = 0;
	int b = 0, e = 1;
	int count = 0;
	char species[100], species2[100];
	bool is_identical = false;
	bool is_identical_other = false;
	bool is_in = false;
	bool is_ps = false;
	int num_paralogs = 0;
	int num_pairs_cr1 = 0;
	int num_pairs_cr2 = 0;
	int num_pairs_cr = 0;
	int num_cv_cr1 = 0;
	int num_cv_cr2 = 0;
	struct I *cv_reg;
	struct I *cv_reg_cr1;
	struct I *cv_reg_cr2;
	struct I *new_cv_reg;
	struct I *new_cv_reg_cr1;
	struct I *new_cv_reg_cr2;
	int num1 = 0, num2 = 0;
	int num_cv_reg;
	int num_cv_reg_cr1;
	int num_cv_reg_cr2;
	struct I *dup_reg;
	struct I *new_dup_reg;
	int num_dup_reg = 0;
	int num_exons = 0;
	struct I *exons;
	struct I *org_exons;
	int num_bp_dup = 0;
	int num_bp_cr1 = 0;
	int num_bp_cr2 = 0;
	int num_bp_cv = 0;
	int num_bp_exons = 0;
	int num_bp_exons_cr1 = 0;
	int num_bp_exons_cr2 = 0;
	int num_bp_exons_cv = 0;
	struct g_list *genes;
	int num_genes = 0;
	int max_id = -1;
	int num_max_events = -1;
	int gene_id = -1;
	int num_cur_exons = 0;
	int num_bp_cur_exons = 0;
	struct I *cur_exons;
	bool is_missing = false;
	bool is_maf_missing = false;

	debug_mode = FALSE;
	if( argc == 5 ) {
		debug_mode = TRUE;
		strcpy(gene_name, argv[4]);
	}
	else if( argc != 4 ) {
		fatal("args: non-redundant.gc annot_directory self-alignment_directory\n");
	}

	num_algns = (int *) ckalloc(sizeof(int));
	size1 = (int *) ckalloc(sizeof(int));
	size2 = (int *) ckalloc(sizeof(int));

	strcpy(buf, "");

	if( (cv_f = fopen(argv[1], "r") ) == NULL ) {
		fatal("all.gc file open error!\n");
	}
	else {
		num_init_cv = 0;
		while( fgets(buf, 1000, cv_f) ) num_init_cv++;

		init_cv = (struct cv_list *) ckalloc(num_init_cv * sizeof(struct cv_list));
		cv = (struct cv_list *) ckalloc(num_init_cv * sizeof(struct cv_list));
		fseek(cv_f, 0, SEEK_SET);
		i = 0;
		while( fgets(buf, 1000, cv_f) ) {
			if( (buf[0] == '(' ) || (buf[0] == '#') )  {
				if( buf[0] == '(') {
   	    	strcpy(tree_line, buf);
        	leave_only_taxons(tree_line);
        	num_sp = 0;
        	j = 0;
        	while( tree_line[j] != '\0') {
          	if( tree_line[j] == ',' ) num_sp++;
          	j++;
        	}

        	num_sp++;
        	if( num_sp > 0 ) {
          	sp_id = (struct sp_list *) ckalloc(num_sp * sizeof(struct sp_list));
          	for( j = 0; j < num_sp; j++ ) {
            	strcpy(sp_id[j].name, "");
            	sp_id[j].id = -1;
          	}
          	temp_num = assign_sp_code(tree_line, sp_id, num_sp);
					}
				}
			}
			else {
				sscanf(buf, "%*s %s %d %d %*s %d %d %c %*s %*s %*s %*s %d %d %d %d %d %s %*s %*s %*s %s %*s %*s %*s %d %*s %*s %*s %d", init_cv[i].pr_name, &init_cv[i].s1, &init_cv[i].s2, &init_cv[i].t1, &init_cv[i].t2, &init_cv[i].ori, &init_cv[i].a1, &init_cv[i].a2, &init_cv[i].b1, &init_cv[i].b2, &init_cv[i].dir, init_cv[i].name1, init_cv[i].name2, &init_cv[i].fid, &init_cv[i].status);
				i++;
			}
		}
	}
	fclose(cv_f);
	num_init_cv = i;
	
	printf("species_name\tnum_genes\tnum_events_crit1\tnum_events_crit2\tnum_events_total\tnum_paralog_pairs\tnum_pairs_crit1\tnum_pairs_crit2\tnum_pairs_conv\tnum_dup_bases\tnum_bases_crit1\tnum_bases_crit2\tnum_bases_conv\tnum_coding_bases\tnum_coding_bases_crit1\tnum_coding_bases_crit2\tnum_coding_bases_conv\n");
	for( l = 0; l < num_sp; l++ ) {
		strcpy(species, sp_id[l].name); 
		num_cv = 0;
		j = 0;
		for( i = 0; i < num_init_cv; i++ ) {
			if( strcmp(species, init_cv[i].pr_name) == 0 )  {
				cv[j].fid = init_cv[i].fid;
				cv[j].s1 = init_cv[i].s1;
				cv[j].s2 = init_cv[i].s2;
				cv[j].t1 = init_cv[i].t1;
				cv[j].t2 = init_cv[i].t2;
				cv[j].ori = init_cv[i].ori;
				cv[j].a1 = init_cv[i].a1;
				cv[j].a2 = init_cv[i].a2;
				cv[j].b1 = init_cv[i].b1;
				cv[j].b2 = init_cv[i].b2;
				cv[j].dir = init_cv[i].dir;
				strcpy(cv[j].name1, init_cv[i].name1);
				cv[i].status = init_cv[i].status;
				if( init_cv[i].status == 4 ) {
					num_cv_cr2++;
				}
				else {
					num_cv_cr1++;
				}
				j++;	
			}
		}
		num_cv = j;

		sprintf(filename, "%s/%s.codex", argv[2], species);
		num_exons = 0;
		num_genes = 0;
  	if( (f = fopen(filename, "r")) == NULL ) {
			num_exons = -1;
			num_genes = -1;
		}
		else {
  		while(fgets(buf, 1000, f))
  		{
   			if( buf[0] == '#' ) {}
				else if( (buf[0] == '<') || (buf[0] == '>') ) num_genes++;
    		else num_exons++;
  		}
		}

		if( num_exons > 0 ) {
			exons = (struct I *) ckalloc(sizeof(struct I) * num_exons);
			org_exons = (struct I *) ckalloc(sizeof(struct I) * num_exons);
		}

		if( num_genes > 0 ) {
			genes = (struct g_list *) ckalloc(sizeof(struct g_list) * num_genes);
			init_g_list(genes, 0, num_genes-1);
		}

		if( num_genes > 0 ) {
			fseek(f, 0, SEEK_SET);
 			j = 0;
			i = 0;
 		 	while(fgets(buf, 1000, f))
 		 	{
 		   	if( buf[0] == '#' ) {}
 		   	else if( (buf[0] == '>') || (buf[0] == '<') ) {
					sscanf(buf, "%*s %d %d %s", &genes[i].cdsStart, &genes[i].cdsEnd, genes[i].gname);
					if( strcmp(genes[i].gname, gene_name) == 0 ) {
						gene_id = i;
					}
						
					is_ps = false;
					if( is_pseudo_gene(genes[i].gname) == false ) {
						genes[i].txStart = j;
						genes[i].sid = 0; // the number of events
						if( i > 0 ) {
							genes[i-1].txEnd = j-1;
						genes[i-1].exonCount = genes[i-1].txEnd - genes[i-1].txStart + 1;
	
						}
						i++;
					}
					else {
						is_ps = true;	
					}
			}
 		   	else if( is_ps == false ) {
 		     		sscanf(buf, "%d %d", &b, &e);
 		     		exons[j] = assign_I(b, e);
					org_exons[j] = assign_I(b, e);
 		     		j++;
 		   	}
 		 	}

			if( (i > 0) && (i <= num_genes ) ) {
				genes[i-1].txEnd = j-1;
			}
			num_genes = i;
			num_exons = j;

			j = 0;
			if( ((gene_id != -1) && (genes[gene_id].txEnd - genes[gene_id].txStart + 1) > 0) ) {
				num_cur_exons = genes[gene_id].txEnd - genes[gene_id].txStart + 1;
				cur_exons = (struct I *) ckalloc(sizeof(struct I) * (num_cur_exons) );
				for( i = genes[gene_id].txStart; i <= genes[gene_id].txEnd; i++ ) {
					cur_exons[j] = assign_I(org_exons[i].lower, org_exons[i].upper);
					j++;
				}
			} 
			quick_sort_inc_int(exons, 0, num_exons-1);
		}	
		if( f != NULL ) fclose(f);

		sprintf(filename, "%s/%s.remove_repeats.maf", argv[3], species);
		count = 0;
		is_maf_missing = false;
		if( (f = fopen(filename, "r")) == NULL ) {
			is_maf_missing = true;
//			fatalf("%s not exist\n", filename);	
		}
		else {
			while(fgets(S, BIG, f)) {
			if( S[0] == '#' ) {
				while( S[0] == '#' ) {
					fgets(S, BIG, f);
					if( strncmp(S, "##maf", 5) == 0 ) {
						count = 0;
					}
				}
				count = 0;	
			}

 		 	if( S[0] == 'a' ) {
				count++;
				if ((fgets(S, BIG, f) == NULL) || (fgets(T, BIG, f) == NULL))      
					fatalf("cannot find alignment in %s", argv[1]);    
				if( (sscanf(S, "%*s %s %*s", species) != 1) || (sscanf(T, "%*s %s %*s", species2) != 1)) {}
			}
		}
		fclose(f);

		num_dup_reg = 0;
		num_bp_dup = 0;	
		if( count > 0 ) {
			algns = (struct DotList *) ckalloc(count * (sizeof(struct DotList)) );
 			read_maf(filename, D_MODE, algns, num_algns, size1, size2);

			num_paralogs = *num_algns;
			dup_reg = (struct I *) ckalloc(2 * num_paralogs * (sizeof(struct I)) );
			new_dup_reg = (struct I *) ckalloc(2 * num_paralogs * (sizeof(struct I)) );
			j = 0;
			for( i = 0; i < num_paralogs; i++ ) {
				dup_reg[j] = assign_I(algns[i].x.lower, algns[i].x.upper-1);	
				dup_reg[j+1] = assign_I(algns[i].y.lower, algns[i].y.upper-1);	
				j = j+2;
			}

			num_dup_reg = sort_merge_intervals(dup_reg, 2*num_paralogs, new_dup_reg);
			free(dup_reg);
	
			num_bp_dup = count_bp(new_dup_reg, num_dup_reg);	
			free(new_dup_reg);

			for( i = 0; i < num_cv; i++ ) {
				src = assign_I(cv[i].a1, cv[i].a2);
				dst = assign_I(cv[i].b1, cv[i].b2);
				for( j = 0; j < num_genes; j++ ) {
					is_in = false;
					for( k = genes[j].txStart; k <= genes[j].txEnd; k++ ) {
						if( (proper_overlap(org_exons[k], src) == true) || (proper_overlap(org_exons[k], dst) == true) ) {
							is_in = true;
						}
					}
					
					if( is_in == true ) {
						genes[j].sid = genes[j].sid + 1;
						if( (debug_mode == TRUE) && (j == gene_id) ) {
							printf("%d %d-%d %d-%d %c %s %d %d ", cv[i].fid, src.lower, src.upper, dst.lower, dst.upper, cv[i].ori, genes[j].gname, genes[j].cdsStart, genes[j].cdsEnd);
							for( k = genes[j].txStart; k <= genes[j].txEnd; k++ ) {
								printf("%d-%d ", org_exons[k].lower, org_exons[k].upper);
							}
							printf("\n");
						}
					}
				}
			}
		}
	
		max_id = -1;
		for( i = 0; i < num_genes; i++ ) {
			if( max_id == -1 ) {
				max_id = i;
				num_max_events = genes[i].sid;
			}
			else {
				if( num_max_events < genes[i].sid ) {
					max_id = i;
					num_max_events = genes[i].sid;
				}
			}
		}

//		if( max_id != -1 ) {
//			printf("%s %d %d %d\n", genes[max_id].gname, genes[max_id].cdsStart, genes[max_id].cdsEnd, genes[max_id].sid);
//		}
//		else {
//			printf("none\n");
//		}

		num_pairs_cr1 = 0; // # of pairs via criterion #1
		num_pairs_cr2 = 0; // # of pairs via criterion #2
		num_pairs_cr = 0;
		num1 = 0;
		num2 = 0;

		if( num_cv > 0 ) {
			cv_reg = (struct I *) ckalloc(sizeof(struct I) * (2*num_cv));
			cv_reg_cr1 = (struct I *) ckalloc(sizeof(struct I) * (2*num_cv));
			cv_reg_cr2 = (struct I *) ckalloc(sizeof(struct I) * (2*num_cv));
		}

		for( i = 0; i < num_cv; i++ ) {
			is_identical = false;
			is_identical_other = false;
			src = assign_I(cv[i].s1, cv[i].s2);
			dst = assign_I(cv[i].t1, cv[i].t2);
			cv_reg[2*i] = assign_I(cv[i].a1, cv[i].a2);
			cv_reg[(2*i)+1] = assign_I(cv[i].b1, cv[i].b2);

			if( cv[i].status == 4 ) {
				cv_reg_cr2[2*num2] = assign_I(cv[i].a1, cv[i].a2);
				cv_reg_cr2[(2*num2)+1] = assign_I(cv[i].b1, cv[i].b2);
				num2++;
			}
			else {
				cv_reg_cr1[2*num1] = assign_I(cv[i].a1, cv[i].a2);
				cv_reg_cr1[(2*num1)+1] = assign_I(cv[i].b1, cv[i].b2);
				num1++;
			}

			for( j = (i+1); j < num_cv; j++ ) {
				cmp1 = assign_I(cv[j].s1, cv[j].s2);
				cmp2 = assign_I(cv[j].t1, cv[j].t2);
				if( (((equal(src, cmp1) == true)	&& (equal(dst, cmp2) == true)) || ((equal(src, cmp2) == true) && (equal(dst, cmp1) == true))) && (cv[i].status == cv[j].status)) {
					is_identical = true;
				}

				if( (((equal(src, cmp1) == true)	&& (equal(dst, cmp2) == true)) || ((equal(src, cmp2) == true) && (equal(dst, cmp1) == true))) ) 
				{
					is_identical_other = true;
				}
			}

			if( is_identical == false ) {
				if( cv[i].status == 4 ) {
					num_pairs_cr2++;
				}
				else {
					num_pairs_cr1++;
				}
			}

			if( is_identical_other == false ) {
				num_pairs_cr++;
			}
		}
	
		num_cv_reg = 0;
		num_cv_reg_cr1 = 0;
		num_cv_reg_cr2 = 0;

		num_bp_cv = 0;
		num_bp_cr1 = 0;
		num_bp_cr2 = 0;

		if( num_cv > 0 ) {
			new_cv_reg = (struct I *) ckalloc(sizeof(struct I) * (2*num_cv));
			num_cv_reg = sort_merge_intervals(cv_reg, 2*num_cv, new_cv_reg);
			num_bp_cv = count_bp(new_cv_reg, num_cv_reg);
		}

		if( num1 > 0 ) {
			new_cv_reg_cr1 = (struct I *) ckalloc(sizeof(struct I) * (2*num1));
			num_cv_reg_cr1 = sort_merge_intervals(cv_reg_cr1, 2*num1, new_cv_reg_cr1);
			num_bp_cr1 = count_bp(new_cv_reg_cr1, num_cv_reg_cr1);
		}

		if( num2 > 0 ) {
			new_cv_reg_cr2 = (struct I *) ckalloc(sizeof(struct I) * (2*num2));
			num_cv_reg_cr2 = sort_merge_intervals(cv_reg_cr2, 2*num2, new_cv_reg_cr2);
			num_bp_cr2 = count_bp(new_cv_reg_cr2, num_cv_reg_cr2);
		}
		
		if( num_cv > 0 ) {
			free(cv_reg);
			free(cv_reg_cr1);
			free(cv_reg_cr2);
		}

		if( num_exons > 0 ) {
			num_bp_exons = 0;
			num_bp_cur_exons = 0;
			num_bp_exons_cr1 = 0;
			num_bp_exons_cr2 = 0;
			num_bp_exons_cv = 0;
			num_bp_exons = count_bp(exons, num_exons);
			if( num_cur_exons > 0 ) {
				num_bp_cur_exons = count_bp(cur_exons, num_cur_exons);
			}
			
			if( num_cv_reg_cr1 > 0 ) {
				num_bp_exons_cr1 = count_overlap_bp(exons, num_exons, new_cv_reg_cr1, num_cv_reg_cr1);
			}
			
			if( num_cv_reg_cr2 > 0 ) {
				num_bp_exons_cr2 = count_overlap_bp(exons, num_exons, new_cv_reg_cr2, num_cv_reg_cr2);
			}

			if( num_cv_reg > 0 ) {
				num_bp_exons_cv = count_overlap_bp(exons, num_exons, new_cv_reg, num_cv_reg);
			}
		}

		if( (num_genes == -1) || (num_exons == -1) ) {
			is_missing = true;
			printf("%s\t*\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t*\t*\t*\t*\n", species, num_cv_cr1, num_cv_cr2, num_cv, num_paralogs, num_pairs_cr1, num_pairs_cr2, num_pairs_cr, num_bp_dup, num_bp_cr1, num_bp_cr2, num_bp_cv);
		}
		else {
			printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", species, num_genes, num_cv_cr1, num_cv_cr2, num_cv, num_paralogs, num_pairs_cr1, num_pairs_cr2, num_pairs_cr, num_bp_dup, num_bp_cr1, num_bp_cr2, num_bp_cv, num_bp_exons, num_bp_exons_cr1, num_bp_exons_cr2, num_bp_exons_cv);
		}
	}
		if( num_genes > 0 ) {
			free(genes);
		}

		if( num_exons > 0 ) {
			free(exons);
			free(org_exons);
		}

		if( num_cur_exons > 0 ) {
			free(cur_exons);
		}

		if( is_maf_missing == false ) {
		if( num_cv > 0 ) free(new_cv_reg);
		if( num1 > 0 ) free(new_cv_reg_cr1);
		if( num2 > 0 ) free(new_cv_reg_cr2);
		if( count > 0 ) free(algns);
		}
	}

	if( is_missing == true ) {
		printf("\n\'*\' means missing gene annotations\n");
	}
	
	free(num_algns);
	free(size1);
	free(size2);
	free(cv);
	free(init_cv);
	free(sp_id);
	return EXIT_SUCCESS;
}

int count_bp(struct I *regs, int num_regs)
{
	int i = 0;
	int sum = 0;
	
	for( i = 0; i < num_regs; i++ ) {
		sum = sum + (regs[i].upper - regs[i].lower + 1);
	}
	
	return(sum);
}

int count_overlap_bp(struct I *reg1, int num_reg1, struct I *reg2, int num_reg2)
{
	int i = 0, j = 0;
	int sum = 0;

	for( i = 0; i < num_reg1; i++ ) {
		for( j = 0; j < num_reg2; j++ ) {
			if( proper_overlap(reg1[i], reg2[j]) == true ) {
//				if( debug_mode == TRUE ) {
//					printf("%d-%d vs %d-%d\n", reg1[i].lower, reg1[i].upper, reg2[j].lower, reg2[j].upper);
//				}
				sum = sum + width(intersect(reg1[i], reg2[j]));
			}
		}
	}
	return(sum);
}

int sort_merge_intervals(struct I *regs, int num_regs, struct I *new_regs) // closed interval
{
	int j = 0;
	int i = 0;
	int lo = 0, hi = 0;
	
  if( num_regs > 0 ) {
    quick_sort_inc_int(regs, 0, num_regs-1);
    j =  0;
    new_regs[j] = assign_I(regs[0].lower, regs[0].upper);
    for( i = 1; i < num_regs; i++ ) {
      while( (i < num_regs) && (overlap(new_regs[j], regs[i]) == true ) ) {
        if( new_regs[j].lower < regs[i].lower ) lo = new_regs[j].lower;
        else lo = regs[i].lower;

        if( new_regs[j].upper < regs[i].upper ) hi = regs[i].upper;
        else hi = new_regs[j].upper;

        new_regs[j] = assign_I(lo, hi);
        i++;
      }

      if( i >= num_regs ) {}
      else {
        j++;
        if( j > num_regs ) {
          fatalf("overflow in new_regs[]: %d\n", j);
        }
        new_regs[j] = assign_I(regs[i].lower, regs[i].upper);
      }
    }

    j++;
  }
  return(j);
}	

int find_status(struct cv_list cur_cv, struct cv_list *cv, int num_cv)
{
	int i = 0;
	struct I src1, dst1, src2, dst2;
	int res = -1;
	char name[50];

	src1 = assign_I(cur_cv.a1, cur_cv.a2);
	dst1 = assign_I(cur_cv.b1, cur_cv.b2);
	for( i = 0; i < num_cv; i++ ) {
		src2 = assign_I(cv[i].a1, cv[i].a2);
		dst2 = assign_I(cv[i].b1, cv[i].b2);

    if( (strcmp( cv[i].name1, "NAN" ) != 0) && (strcmp( cv[i].name2, "NAN") != 0) ) {
      strcpy( name, cv[i].name1 );
    }
    else if( (strcmp( cv[i].name1, "NAN" ) == 0) && (strcmp( cv[i].name2, "NAN") != 0) ) {
      strcpy( name, cv[i].name2 );
    }    else if( (strcmp( cv[i].name2, "NAN" ) == 0) && (strcmp( cv[i].name1, "NAN") != 0) ) {
      strcpy( name, cv[i].name1 );
    }
    else {
      fatalf("both out-species not found %s %s\n", cv[i].name1, cv[i].name2);
    }

		if( (cur_cv.fid == cv[i].fid) && ( ((strict_almost_equal(src1, src2) == true) && (strict_almost_equal(dst1, dst2) == true)) || ( (strict_almost_equal(src1, dst2) == true) && (strict_almost_equal(dst1, src2) == true) )) && (strcmp(cur_cv.name1, name) == 0) ) {
			res = cv[i].status;
		}
	}

	if( res == -1 ) {
		fatalf("status not found\n");
	}
	return(res);
}

bool is_pseudo_gene(char *name)
{
	bool res = false;
	int len = 0;
	char end_str[4];

	strcpy(end_str, "");
	len = strlen(name);
	if( len > 3 ) {
		end_str[2] = name[len-1];
		end_str[1] = name[len-2];
		end_str[0] = name[len-3];
		end_str[3] = '\0';
	}

	if( strcmp(end_str, "_ps") == 0 ) {
		res = true;
	}

	return(res);
}

void init_g_list(struct g_list *genes, int b, int e)
{
	int i = 0;

	for( i = b; i <=e; i++ ) 
	{
		genes[i].gid = 0;
		genes[i].sid = 0;
		genes[i].strand = '+';
		genes[i].txStart = 0;
		genes[i].txEnd = 0;
		genes[i].cdsStart = 0;
		genes[i].cdsEnd = 0;
		genes[i].exonCount = 0;
		genes[i].ortho_id = 0;
		strcpy(genes[i].gname, "");
	}
}
