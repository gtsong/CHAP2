#include "main.h"
#include "regions.h"
#include "util.h"
#include "util_i.h"
#include "tree_op.h"

int debug_mode;

int find_status(struct cv_list cur_cv, struct cv_list *cv, int num_cv, char *name1, char *name2);

int main(int argc, char **argv)
{
	FILE *cv_f;
	int *bid;
 	struct cv_list *cv; // the conversion events detected by Chih-Hao's program
	int num_cv;
 	struct cv_list *init_cv; // the conversion events detected by Chih-Hao's program
	int num_init_cv;
	char buf[1000];
	char tree_line[1000];
	char name_line[1000];
	int i = 0, j = 0, k = 0;
//	char species1[100], species2[100];
//	int res = 0;
	char name1[100], name2[100];
//	struct I src1, dst1, src2, dst2;
	int event_number = 0;
	struct p_tree *sp_tree = NULL;
	struct p_tree *lca = NULL;
	struct p_tree *t = NULL;
	int nid = 0;
	int root_id = 0;
	struct sp_list *sp_id; 	
	int num_sp = 0, temp_num = 0;
	int *primary_sp;
	int num_primary = 0;
	bool is_in = false;
	int temp_bid[50];
	int num_bid = 0;

	debug_mode = FALSE;
	strcpy(name_line, "");
	strcpy(tree_line, "");
	if( argc == 5 ) debug_mode = TRUE;
	else if( argc == 4 ) {
		event_number = atoi(argv[3]);
	}
	else if( argc != 3 ) {
		fatal("args: all.gc all.remove_redundancy.gc [event_id]\n");
	}
	else {
		event_number = -1;
	}

	bid = (int *) ckalloc(sizeof(int));
	cv_f = fopen(argv[1], "r");
	if( cv_f == NULL ) {
		fatalf("conv-file open error %s!\n", argv[1]);
	}
	else {
		num_init_cv = 0;
		while( fgets(buf, 1000, cv_f) ) num_init_cv++;

		init_cv = (struct cv_list *) ckalloc(num_init_cv * sizeof(struct cv_list));
		fseek(cv_f, 0, SEEK_SET);
		i = 0;
		while( fgets(buf, 1000, cv_f) ) {
			if( (buf[0] == '(' ) || (buf[0] == '#') )  {
				if( buf[0] == '(' ) {
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
						if( temp_num != num_sp ) {
							fatalf("Number of species mismatched in %d and %d", num_sp, temp_num);
						}
						sp_tree = read_one_line(tree_line);
						root_id = sp_tree->nid;
						assign_sp_id(sp_tree, sp_id, 0, num_sp);
						if( debug_mode == TRUE ) {
							print_tree(sp_tree, TREE_PRINT);
							printf("\n");
						}
					}
				}
				else if( buf[0] == '#' ) {
					strcpy(name_line, buf);
				}
			}
			else {
				if( event_number == -1 ) {
					sscanf(buf, "%d %s %d %d %*s %d %d %c %d %f %d %s %d %d %d %d %d %s %d %d %c %s %d %d %c %d %d %s %s %d", &init_cv[i].oid, init_cv[i].name1, &init_cv[i].s1, &init_cv[i].s2, &init_cv[i].t1, &init_cv[i].t2, &init_cv[i].ori, &init_cv[i].len1, &init_cv[i].pid, &init_cv[i].len2, init_cv[i].pval, &init_cv[i].a1, &init_cv[i].a2, &init_cv[i].b1, &init_cv[i].b2, &init_cv[i].dir, init_cv[i].name2, &init_cv[i].c1, &init_cv[i].c2, &init_cv[i].ori1, init_cv[i].name3, &init_cv[i].d1, &init_cv[i].d2, &init_cv[i].ori2, &init_cv[i].fid, &init_cv[i].bid1, init_cv[i].ortho1, init_cv[i].ortho2, &init_cv[i].status);
					init_cv[i].num_bid2 = 0; // marking unvisited
					i++;
				}
				else {
					sscanf(buf, "%d %s %d %d %*s %d %d %c %d %f %d %s %d %d %d %d %d %s %d %d %c %s %d %d %c %d %*s %s %s %d", &init_cv[i].oid, init_cv[i].name1, &init_cv[i].s1, &init_cv[i].s2, &init_cv[i].t1, &init_cv[i].t2, &init_cv[i].ori, &init_cv[i].len1, &init_cv[i].pid, &init_cv[i].len2, init_cv[i].pval, &init_cv[i].a1, &init_cv[i].a2, &init_cv[i].b1, &init_cv[i].b2, &init_cv[i].dir, init_cv[i].name2, &init_cv[i].c1, &init_cv[i].c2, &init_cv[i].ori1, init_cv[i].name3, &init_cv[i].d1, &init_cv[i].d2, &init_cv[i].ori2, &init_cv[i].fid, init_cv[i].ortho1, init_cv[i].ortho2, &init_cv[i].status);
					if( init_cv[i].fid == event_number ) {
						printf("%s", buf);
					}
					i++;
				}
			}
		}
	}
	fclose(cv_f);
	num_init_cv = i;

	primary_sp = ckalloc(sizeof(int) * num_sp);

	if( event_number == -1 ) {
		for( i = 0; i < num_init_cv; i++ ) {
			if( init_cv[i].num_bid2 == 0 ) { // init_cv[i] unvisited
				strcpy( name1, init_cv[i].name1 );
				num_primary = 0;
				primary_sp[num_primary] = get_id_in_list(name1, sp_id, num_sp);
				num_primary++;
			}

			for( j = (i+1); j < num_init_cv; j++ ) {
				if( (init_cv[j].num_bid2 == 0) && (init_cv[j].fid == init_cv[i].fid) ) {
					strcpy(name2, init_cv[j].name1);
					temp_num = get_id_in_list(name2, sp_id, num_sp);
					is_in = false;
					for(k = 0;  k < num_primary; k++) {
						if( temp_num == primary_sp[k] ) is_in = true;
					}
					if( is_in == false ) {
						if( num_primary >= num_sp ) {
							fatalf("overflow in primary_sp[], %s\n", name2);
						}
						else {
							primary_sp[num_primary] = temp_num;
							num_primary++;
						}
					}
				}
			}
			lca = find_lca(primary_sp, num_primary, sp_tree, &nid);

			if( (lca->nid) > init_cv[i].bid1 ) {
				num_bid = 1;
				temp_bid[k] = -1;
			}
			else if( lca != NULL ) 
			{
				t = lca;
				k = 0;
				temp_num = t->nid;
				temp_bid[k] = t->nid;
				k++;
				while( (k < 50) && (t != sp_tree) && ((t->nid) < root_id) && ((t->nid) <init_cv[i].bid1) ) 
				{
					t = t->parent;
					temp_bid[k] = t->nid;
					k++;
				}
				num_bid = k;
			}
			else {
				fatal("lca not found\n");	
			}

			if( k >= 50 ) {
				printf("Warning: exceed limit 50 branch numbers\n");
			}

			for( j = i; j < num_init_cv; j++ ) {
				if( (init_cv[j].num_bid2 == 0) && (init_cv[j].fid == init_cv[i].fid) ) 
				{
					init_cv[j].num_bid2 = num_bid;
					for( k = 0; k < num_bid; k++ ) {
						init_cv[j].bid2[k] = temp_bid[k];
					}
				}
			}
		}
	}

	cv_f = fopen(argv[2], "r");
	if( cv_f == NULL ) {
		fatalf("conv-file open error %s!\n", argv[2]);
	}
	else {
		num_cv = 0;
		while( fgets(buf, 1000, cv_f) ) num_cv++;

		cv = (struct cv_list *) ckalloc(num_cv * sizeof(struct cv_list));
		fseek(cv_f, 0, SEEK_SET);
		i = 0;
/*
		while( fgets(buf, 1000, cv_f) ) {
			if( (buf[0] == '(' ) || (buf[0] == '#') )  {}
			else {
				sscanf(buf, "%d %s %d %d %*s %d %d %c %*s %*s %*s %*s %d %d %d %d %d %s", &cv[i].fid, species1, &cv[i].s1, &cv[i].s2, &cv[i].t1, &cv[i].t2, &cv[i].ori, &cv[i].a1, &cv[i].a2, &cv[i].b1, &cv[i].b2, &cv[i].dir, species2);
				res = find_status(cv[i], init_cv, num_init_cv, species1, species2);

				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%.2f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%d\t%d\t%s\t%s\t%d\n", init_cv[res].oid, init_cv[res].name1, init_cv[res].s1, init_cv[res].s2, init_cv[res].name1, init_cv[res].t1, init_cv[res].t2, init_cv[res].ori, init_cv[res].len1, init_cv[res].pid, init_cv[res].len2, init_cv[res].pval, init_cv[res].a1, init_cv[res].a2, init_cv[res].b1, init_cv[res].b2, init_cv[res].dir, init_cv[res].name2, init_cv[res].c1, init_cv[res].c2, init_cv[res].ori1, init_cv[res].name3, init_cv[res].d1, init_cv[res].d2, init_cv[res].ori2, init_cv[res].fid, init_cv[res].bid1, init_cv[res].ortho1, init_cv[res].ortho2, init_cv[res].status);
				i++;
			}
		}
*/
	}
	fclose(cv_f);

	if( event_number == -1 ) {
		printf("%s\n", tree_line);
		printf("%s", name_line);
		for( i = 0; i < num_init_cv; i++ ) {
			if( (init_cv[i].num_bid2 == 0) || ((init_cv[i].num_bid2 == 1) && (init_cv[i].bid2[0] == -1)) )  
			{
				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%.2f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%d\t%d,?\t%s\t%s\t%d\n", init_cv[i].oid, init_cv[i].name1, init_cv[i].s1, init_cv[i].s2, init_cv[i].name1, init_cv[i].t1, init_cv[i].t2, init_cv[i].ori, init_cv[i].len1, init_cv[i].pid, init_cv[i].len2, init_cv[i].pval, init_cv[i].a1, init_cv[i].a2, init_cv[i].b1, init_cv[i].b2, init_cv[i].dir, init_cv[i].name2, init_cv[i].c1, init_cv[i].c2, init_cv[i].ori1, init_cv[i].name3, init_cv[i].d1, init_cv[i].d2, init_cv[i].ori2, init_cv[i].fid, init_cv[i].bid1, init_cv[i].ortho1, init_cv[i].ortho2, init_cv[i].status);
			}
			else {
				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%.2f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%d\t", init_cv[i].oid, init_cv[i].name1, init_cv[i].s1, init_cv[i].s2, init_cv[i].name1, init_cv[i].t1, init_cv[i].t2, init_cv[i].ori, init_cv[i].len1, init_cv[i].pid, init_cv[i].len2, init_cv[i].pval, init_cv[i].a1, init_cv[i].a2, init_cv[i].b1, init_cv[i].b2, init_cv[i].dir, init_cv[i].name2, init_cv[i].c1, init_cv[i].c2, init_cv[i].ori1, init_cv[i].name3, init_cv[i].d1, init_cv[i].d2, init_cv[i].ori2, init_cv[i].fid);
				printf("%d", init_cv[i].bid2[0]);
				for( k = 1; k < init_cv[i].num_bid2; k++ ) {
					printf(",%d", init_cv[i].bid2[k]);
				}
				printf("\t%s\t%s\t%d\n", init_cv[i].ortho1, init_cv[i].ortho2, init_cv[i].status);
			}
		}
	}

	free_p_tree(sp_tree);
	free(primary_sp);
	free(sp_id);
	free(bid);
	free(init_cv);
	return EXIT_SUCCESS;
}

int find_status(struct cv_list cur_cv, struct cv_list *cv, int num_cv, char *name1, char *name2)
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
		if( (strcmp( cv[i].name2, "NAN" ) != 0) && (strcmp( cv[i].name3, "NAN") != 0) ) {
			strcpy( name, cv[i].name2 );
		}
		else if( (strcmp( cv[i].name2, "NAN" ) == 0) && (strcmp( cv[i].name3, "NAN") != 0) ) {
			strcpy( name, cv[i].name3 );
		}
		else if( (strcmp( cv[i].name3, "NAN" ) == 0) && (strcmp( cv[i].name2, "NAN") != 0) ) {
			strcpy( name, cv[i].name2 );
		}
		else {
			fatalf("both out-species not found %s %s\n", cv[i].name2, cv[i].name3);
		}

		if( (cur_cv.fid == cv[i].fid) && ( ((strict_almost_equal(src1, src2) == true) && (strict_almost_equal(dst1, dst2) == true)) || ( (strict_almost_equal(src1, dst2) == true) && (strict_almost_equal(dst1, src2) == true) )) && (strcmp(name1, cv[i].name1) == 0) && (strcmp(name2, name) == 0)) {
			res = i;
		}
	}

	if( res == -1 ) {
		fatalf("status not found %d\n", cur_cv.fid);
	}
	return(res);
}
