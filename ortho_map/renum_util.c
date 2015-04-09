#include "main.h"
#include "renumber.h"
#include "renum_util.h"
#include "util.h"
#include "util_gen.h"
#include "util_input.h"
#include "tree_op.h"

void read_view_file(char *fname, struct mapping_list *ortho_maps)
{
	FILE *f;
	char buf[LEN_NAME];
	char token[LEN_NAME];
	char piece[LEN_NAME];
	int i = 0, j = 0, k = 0, n = 0, loc = 0;
	int num_genes = 0;
	int num_sp = 0, len = 0, temp_len = 0;
	int num_ids = 0;
	int num_parts = 0;
	
	if( (f = ckopen(fname, "r")) != NULL ) {
		fgets(buf, LEN_NAME, f); // the first line is a list of human gene names
		while(fgets(buf, LEN_NAME, f)) {
			i = 0;
			j = 0;
			num_genes = count_tokens(buf);
			num_genes--;
			ortho_maps[num_sp].num_genes = num_genes;
			if( num_genes > 0 ) {
				ortho_maps[num_sp].o_id = (struct o_list *) ckalloc(num_genes * sizeof(struct o_list));
				init_o_list(ortho_maps[num_sp].o_id, num_genes);
				len = strlen(buf);
				k = 0;
				while( (j < (num_genes+1)) && (i < len) && (buf[i] != '\0') && (buf[i] != '\n')) {
					i = concat_tokens(buf, i, token);	

					if( j == 0 ) {
						strcpy(ortho_maps[num_sp].sp_name, token);
					}
					else {
						num_parts = count_pieces(token);
						if( num_parts > 0 ) {
							ortho_maps[num_sp].o_id[k].num_parts = num_parts;
							ortho_maps[num_sp].o_id[k].pieces = (struct part_list *) ckalloc(num_parts * sizeof(struct part_list));
							init_part_list(ortho_maps[num_sp].o_id[k].pieces, num_parts);
							n = 0;
							loc = 0;
							temp_len = strlen(token);
							while( (n < (num_parts)) && (loc < temp_len) && (token[loc] != '\0') && (token[loc] != '\n') ) {
								loc = concat_piece(token, loc, piece);
								num_ids = count_id(piece);
								if( num_ids > 0 ) {
									ortho_maps[num_sp].o_id[k].pieces[n].num_ids = num_ids;
									ortho_maps[num_sp].o_id[k].pieces[n].ids = (int *) ckalloc(num_ids * sizeof(int));

									initialize_int_list(ortho_maps[num_sp].o_id[k].pieces[n].ids, 0, num_ids-1);
									if( piece[0] == '\\' ) {
										ortho_maps[num_sp].o_id[k].is_pseudo = true;
										assign_ids(piece+1, ortho_maps[num_sp].o_id[k].pieces[n].ids, num_ids);
									}
									else assign_ids(piece, ortho_maps[num_sp].o_id[k].pieces[n].ids, num_ids);
									quick_sort_num_list_inc(ortho_maps[num_sp].o_id[k].pieces[n].ids, 0, num_ids-1);
									n++;
								}
								else {
									fatalf("incorrect piece: %s\n", piece);
								}
							}
							k++;
						}
						else {
							fatalf("incorrect token: %s\n", token);
						}
					}
					j++;
				}
			}
			else if ( num_genes == 0 ) {
				i = concat_tokens(buf, i, token);
				if( j == 0 ) {
					strcpy(ortho_maps[num_sp].sp_name, token);
				}
			}
			num_sp++;
		}
		fclose(f);
	}
}

void assign_ids(char *token, int *ids, int num_ids)
{
	int len = 0, i = 0, j = 0, k = 0;
	char temp[LEN_NAME];
	int temp_id = 0;

	len = strlen(token);

	while( (i < len) && (token[i] != '\0') && (token[i] != '\n') ) {
		if( token[i] == ',' ) {
			temp[j] = '\0';
			temp_id = atoi(temp);
			if( k >= num_ids ) {
				fatalf("%d is over %d\n", k, num_ids);
			}
			ids[k] = temp_id;
			k++;
			j = 0;
		}
		else {
			temp[j] = token[i];	
			j++;
		}
		i++;
	}
	temp[j] = '\0';
	temp_id = atoi(temp);
	if( k >= num_ids ) {
		fatalf("%d is over %d\n", k, num_ids);
	}
	ids[k] = temp_id;
}

int count_id(char *token) // token is separated by ','
{
	int len = 0, i = 0;
	int count = 0;

	len = strlen(token);

	while( (i < len) && (token[i] != '\0') && (token[i] != '\n') ) {
		if( token[i] == ',' ) count++;
		i++;
	}
	count++;
	return(count);
}

int count_pieces(char *token) // token is separated by ','
{
	int len = 0, i = 0;
	int count = 0;

	len = strlen(token);

	while( (i < len) && (token[i] != '\0') && (token[i] != '\n') ) {
		if( token[i] == ';' ) count++;
		i++;
	}
	return(count);
}

void init_ortho_maps(struct mapping_list * maps, int num_list)
{
  int i = 0;

  for( i = 0; i < num_list; i++ ) {
    maps[i].id = -1;
		strcpy(maps[i].sp_name, "");
		maps[i].num_genes = 0;
		maps[i].o_id = NULL;
  }
}

void init_o_list(struct o_list * oids, int num_oids)
{
  int i = 0;

  for( i = 0; i < num_oids; i++ ) {
    oids[i].id = -1;
		oids[i].num_parts = 0;
		oids[i].pieces = NULL;
		oids[i].is_pseudo = false;
  }
}

void init_part_list(struct part_list * pieces, int num_parts)
{
	int i = 0;

	for( i = 0; i < num_parts; i++ ) {
		pieces[i].num_ids = 0;
		pieces[i].ids = NULL;
	}
}

void print_view_file(struct p_tree *t, struct sp_list *genes, int num_genes, struct mapping_list *ortho_maps, int num_list)
{
	int i = 0;

	printf("gene_name");
	for( i = 0; i < num_genes; i++ ) {
		printf("\t%s", genes[i].name);
	}
	printf("\n");

	postorder_print(t, ortho_maps, num_list);	
}

void postorder_print(struct p_tree *t, struct mapping_list *ortho_maps, int num_list)
{
	int i = 0, k = 0, j = 0, n = 0;
	if( t-> left != NULL ) postorder_print(t->left, ortho_maps, num_list);
	if( t-> right != NULL ) postorder_print(t->right, ortho_maps, num_list);
	if( (is_leaf_node(t) == true) && (t->sp_code != -1) ) {
		i = t->sp_code;
		printf("%s", ortho_maps[i].sp_name);
		for( j = 0; j < ortho_maps[i].num_genes; j++ ) {
			if( ortho_maps[i].o_id[j].is_pseudo == true ) {
				printf("\t\\");
			}
			else printf("\t");
	
			for( k = 0; k < ortho_maps[i].o_id[j].num_parts; k++ ) {
				for( n = 0; n < ortho_maps[i].o_id[j].pieces[k].num_ids; n++ ) {
					if( n == 0 ) {
						printf("%d", ortho_maps[i].o_id[j].pieces[k].ids[n]);
					}
					else {
						if( ortho_maps[i].o_id[j].pieces[k].ids[n] != ortho_maps[i].o_id[j].pieces[k].ids[n-1] ) {
							printf(",%d", ortho_maps[i].o_id[j].pieces[k].ids[n]);
						}
					}
				}
				printf(";");
			}
		}
		printf("\n");
	}
}

void assign_sp_id_mapping_list(struct p_tree *t, struct mapping_list *list, int depth, int num_sp)
{
  struct p_tree *temp = NULL;
  int count;
  int i;

  if( t->left != NULL ) assign_sp_id_mapping_list(t->left, list, depth+1, num_sp);
  if( t->right != NULL ) assign_sp_id_mapping_list(t->right, list, depth+1, num_sp);

  if( is_leaf_node(t) || is_done_node(t) ) {
      t->depth = depth;
      t->sp_code = get_id_in_mapping_list(t->name, list, num_sp);
      t->num_csp = 1;
      t->ch_sp = (int *) ckalloc(sizeof(int));
      t->ch_sp[0] = t->sp_code; // the initialization of a species list
//      if( debug_mode == true) printf("%d: %s %d\n", t->nid, t->name, t->sp_code);
  }
  else
  {
    t->depth = depth;
    temp = t->left;
    count = 0;
    if(  ( (t->left->num_csp) + (t->right->num_csp) ) < num_sp ) {
      t->ch_sp = (int *) ckalloc(((t->left->num_csp) + (t->right->num_csp)) * sizeof(int));
    }
    else {
      t->ch_sp = (int *) ckalloc(num_sp * sizeof(int));
    }

    for( i = 0; i < temp->num_csp; i++ ) {
      t->ch_sp[count] = temp->ch_sp[i];
      count++;
    }
    temp = t->right;
    for( i = 0; i < temp->num_csp; i++ ) {
      if( is_in_list(temp->ch_sp[i], t->ch_sp, count) == false ) {
        t->ch_sp[count] = temp->ch_sp[i];
        count++;
      }
    }
    t->num_csp = count;
  }
}

void flip_tree(struct p_tree *t, int ref_id) 
{
	struct p_tree *temp = NULL;

	if( is_leaf_node(t) ) {}
	else { 
		if(is_in_list(ref_id, t->left->ch_sp, t->left->num_csp) == false)
		{
			temp = t->left;
			t->left = t->right;
			t->right = temp;	
			flip_tree(t->left, ref_id);
		}
		else {
			flip_tree(t->left, ref_id);
		}
	}
}

int get_id_in_mapping_list(char *name, struct mapping_list *list, int num_sp)
{
  int i;

  i = 0;
  while( ( i < num_sp) && (strcmp(name, list[i].sp_name) != 0) ) i++;

  if( i >= num_sp ) {
//    fatalf("error: no species name exist in %s\n", name);
//    return(num_sp);
		return(-1);
  }
  else return(i);
}

int concat_piece(char *line, int loc, char *temp_name)
{
  int len = 0, i = 0, j = 0;

  len = strlen(line);
  i = loc;
  if( line[i] == ';' ) {
    while((line[i] != '\0') && (line[i] != '\n') && (line[i] == ';')) i++;
  }

  while((line[i] != '\0') && (line[i] != '\n') && (line[i] != ';'))
  {
    temp_name[j] = line[i];
    j++;
    if( j >= LEN_NAME ) {
      fatalf("Over the max length for a gene name ( j characters ) in %s\n", line);
    }
    i++;
  }
  temp_name[j] = '\0';

	if( line[i] == ';' ) i++;

  return(i);
}

