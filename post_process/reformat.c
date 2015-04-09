#include "main.h"
#include "regions.h"
#include "util.h"
#include "util_i.h"

int debug_mode;

void cal_real_pos(struct n_pair *contigs, int num_contigs, int *loc1, int *loc2,  char *name1, char *name2);
int find_species_id(char *name, struct sp_list *species, int num_sp);
int find_status(struct cv_list cur_cv, struct cv_list *cv, int num_cv, char *name1, char *name2);

int main(int argc, char **argv)
{
	FILE *cv_f;
	int *bid;
 	struct cv_list *cv; // the conversion events detected by Chih-Hao's program
	int num_cv;
 	struct cv_list *init_cv; // the conversion events detected by Chih-Hao's program
	int num_init_cv;
	char buf[2000];
	int i = 0, j = 0;
	char species1[100], species2[100];
	int res = 0;
	DIR *dir;
	struct dirent *dent;
	int num_sp = 0;
	struct n_pair **contigs;
	struct sp_list *species;
	FILE *fp;
	char name1[100], name2[100];
	char name3[100], name4[100];
	char name5[100], name6[100];
	char temp_name[100];
	int len = 0, len_sum = 0;
	char file_name[100];

	debug_mode = FALSE;
	if( argc == 5 ) debug_mode = TRUE;
	else if( argc != 4 ) {
		fatal("args: all.gc all.remove_redundancy.gc contigs_list_directory\n");
	}

	if((dir = opendir(argv[3])) != NULL) {
		while((dent = readdir(dir))) num_sp++;

		contigs = (struct n_pair **) ckalloc(sizeof(struct n_pair *) * num_sp);
		species = (struct sp_list *) ckalloc(sizeof(struct sp_list) * num_sp);

		seekdir(dir, 0);
		j = 0;
		while( (dent = readdir(dir)) ) {
			if( dent->d_name[0] != '.' ) {
				strcpy(file_name, argv[3]);
				strcat(file_name, "/");
				strcat(file_name, dent->d_name);
				fp = ckopen(file_name, "r");
				i = 0;
				while(fgets(buf, 2000, fp)) i++;	
				contigs[j] = (struct n_pair *) ckalloc(sizeof(struct n_pair) * i);
				species[j].id = i;
				fseek(fp, 0, SEEK_SET);
				i = 0;
				while(fgets(buf, 2000, fp)) {
					sscanf(buf, "%s %s %d %d", name1, name2, &len, &len_sum);
					strcpy(contigs[j][i].name1, name1);
					strcpy(contigs[j][i].name2, name2);
					contigs[j][i].id = len;
					contigs[j][i].len = len_sum;
					i++;
				}
				strcpy(species[j].name, name1);
				j++;
			}
		}
	}	

	fclose(fp);
	closedir(dir);

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
//				if( buf[0] == '(' ) {    // keep column headers too  [-CR, 2/2011]
					printf("%s", buf);
//				}
			}
			else {
				sscanf(buf, "%d %s %d %d %*s %d %d %c %d %f %d %s %d %d %d %d %d %s %d %d %c %s %d %d %c %d %s %s %s %d", &init_cv[i].oid, init_cv[i].name1, &init_cv[i].s1, &init_cv[i].s2, &init_cv[i].t1, &init_cv[i].t2, &init_cv[i].ori, &init_cv[i].len1, &init_cv[i].pid, &init_cv[i].len2, init_cv[i].pval, &init_cv[i].a1, &init_cv[i].a2, &init_cv[i].b1, &init_cv[i].b2, &init_cv[i].dir, init_cv[i].name2, &init_cv[i].c1, &init_cv[i].c2, &init_cv[i].ori1, init_cv[i].name3, &init_cv[i].d1, &init_cv[i].d2, &init_cv[i].ori2, &init_cv[i].fid, init_cv[i].bid, init_cv[i].ortho1, init_cv[i].ortho2, &init_cv[i].status);
				i++;
			}
		}
	}
	fclose(cv_f);
	num_init_cv = i;
	
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
		while( fgets(buf, 1000, cv_f) ) {
			if( (buf[0] == '(' ) || (buf[0] == '#') )  {}
			else {
				sscanf(buf, "%d %s %d %d %*s %d %d %c %*s %*s %*s %*s %d %d %d %d %d %s", &cv[i].fid, species1, &cv[i].s1, &cv[i].s2, &cv[i].t1, &cv[i].t2, &cv[i].ori, &cv[i].a1, &cv[i].a2, &cv[i].b1, &cv[i].b2, &cv[i].dir, species2);
				res = find_status(cv[i], init_cv, num_init_cv, species1, species2);

				j = find_species_id(init_cv[res].name1, species, num_sp);
				if( j == -1 ) {
					fatalf("%s not found in the contigs list\n", init_cv[res].name1);
				}

				if( strcmp( init_cv[res].name1, contigs[j][0].name1 ) != 0 ) {
					fatalf("species names not match: %s vs %s\n", init_cv[res].name1, contigs[j][0].name1);
				}

				cal_real_pos(contigs[j], species[j].id, &init_cv[res].s1, &init_cv[res].s2, name1, name1);
				cal_real_pos(contigs[j], species[j].id, &init_cv[res].a1, &init_cv[res].a2, name1, name1);
				strcpy(temp_name, init_cv[res].name1);
				if( strcmp( temp_name, name1) != 0 ) {
					strcat(temp_name, ".");
					strcat(temp_name, name1);
				}
				strcpy(name1, temp_name);

				cal_real_pos(contigs[j], species[j].id, &init_cv[res].t1, &init_cv[res].t2, name2, name2);
				cal_real_pos(contigs[j], species[j].id, &init_cv[res].b1, &init_cv[res].b2, name2, name2);
				strcpy(temp_name, init_cv[res].name1);
				if( strcmp( temp_name, name2) != 0 ) {
					strcat(temp_name, ".");
					strcat(temp_name, name2);
				}
				strcpy(name2, temp_name);

				if( strcmp( init_cv[res].name2, "NAN" ) != 0 ) {
					j = find_species_id(init_cv[res].name2, species, num_sp);
					if( j == -1 ) {
						fatalf("%s not found in the contigs list\n", init_cv[res].name2);
					}
					cal_real_pos(contigs[j], species[j].id, &init_cv[res].c1, &init_cv[res].c2, name3, name4);
					strcat(name3, ":");
					strcat(name4, ":");
				}
				else {
					strcpy(name3, "");
					strcpy(name4, "");
				}

				if( strcmp( init_cv[res].name3, "NAN" ) != 0 ) {
					j = find_species_id(init_cv[res].name3, species, num_sp);
					if( j == -1 ) {
						fatalf("%s not found in the contigs list\n", init_cv[res].name3);
					}
					cal_real_pos(contigs[j], species[j].id, &init_cv[res].d1, &init_cv[res].d2, name5, name6);
					strcat(name5, ":");
					strcat(name6, ":");
				}
				else {
					strcpy(name5, "");
					strcpy(name6, "");
				}

//				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%.2f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%d\t%s\t%s\t%s\t%d\n", init_cv[res].oid, init_cv[res].name1, init_cv[res].s1, init_cv[res].s2, init_cv[res].name1, init_cv[res].t1, init_cv[res].t2, init_cv[res].ori, init_cv[res].len1, init_cv[res].pid, init_cv[res].len2, init_cv[res].pval, init_cv[res].a1, init_cv[res].a2, init_cv[res].b1, init_cv[res].b2, init_cv[res].dir, init_cv[res].name2, init_cv[res].c1, init_cv[res].c2, init_cv[res].ori1, init_cv[res].name3, init_cv[res].d1, init_cv[res].d2, init_cv[res].ori2, init_cv[res].fid, init_cv[res].bid, init_cv[res].ortho1, init_cv[res].ortho2, init_cv[res].status);    // added tabs and float precision to match input  [-CR, 2/2011]
				printf("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%.2f\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s%d\t%s%d\t%c\t%s\t%s%d\t%s%d\t%c\t%d\t%s\t%s\t%s\t%d\n", init_cv[res].oid, name1, init_cv[res].s1, init_cv[res].s2, name2, init_cv[res].t1, init_cv[res].t2, init_cv[res].ori, init_cv[res].len1, init_cv[res].pid, init_cv[res].len2, init_cv[res].pval, init_cv[res].a1, init_cv[res].a2, init_cv[res].b1, init_cv[res].b2, init_cv[res].dir, init_cv[res].name2, name3, init_cv[res].c1, name4, init_cv[res].c2, init_cv[res].ori1, init_cv[res].name3, name5, init_cv[res].d1, name6, init_cv[res].d2, init_cv[res].ori2, init_cv[res].fid, init_cv[res].bid, init_cv[res].ortho1, init_cv[res].ortho2, init_cv[res].status);    // added tabs and float precision to match input  [-CR, 2/2011]
				i++;
			}
		}
	}
	fclose(cv_f);


	for( i = 0; i < num_sp; i++ ) {
		free(contigs[i]);
	}
	
	if( num_sp > 0 ) {
		free(species);
		free(contigs);
	}

	free(bid);
	free(init_cv);
	return EXIT_SUCCESS;
}

void cal_real_pos(struct n_pair *contigs, int num_contigs, int *loc1, int *loc2,  char *name1, char *name2)
{
	int cur_loc1 = 0, cur_loc2 = 0;	
	int i = 0;
	int res1 = -1, res2 = -1;
	int len_sum;

	cur_loc1 = *loc1;
	cur_loc2 = *loc2;

	i = 0;
	while( (res1 == -1) && (i < num_contigs) ) {
		if( cur_loc1 <= contigs[i].len ) {
			res1 = i - 1;
		}
		i++;
	}	
	if( i == num_contigs ) res1 = i - 1;
	
	i = 0;
	while( (res2 == -1) && (i < num_contigs) ) {
		if( cur_loc2 <= contigs[i].len ) {
			res2 = i - 1;
		}
		i++;
	}	
	if( i == num_contigs) res2 = i - 1;

	if( (res1 < 0) || (res1 >= num_contigs) ) {
		fatalf("invalid index %d in reformatting\n", res1);
	}

	if( (res2 < 0) || (res2 >= num_contigs) ) {
		fatalf("invalid index %d in reformatting\n", res2);
	}

	len_sum = contigs[res1].len;
	strcpy(name1, contigs[res1].name2);
	*loc1 = cur_loc1 - len_sum;

	len_sum = contigs[res2].len;
	strcpy(name2, contigs[res2].name2);
	*loc2 = cur_loc2 - len_sum;
}

int find_species_id(char *name, struct sp_list *species, int num_sp)
{
	int res = -1;
	int i = 0;

	while( (res == -1) && (i < num_sp) ) {
		if( strcmp(name, species[i].name) == 0 ) {
			res = i;
		}
		i++;
	} 

	return(res);
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
