#include "main.h"
#include "renumber.h"
#include "renum_util.h"
#include "regions.h"
#include "util.h"
#include "util_i.h"
#include "util_input.h"
#include "tree_op.h"

#define DPI 72
#define WIDTH 6
#define GFS 12
#define SFS 14
#define MIN_GFS 9
#define MIN_SFS 9
#define GSIZE 30

int debug_mode;
int count_node; 
char S[BIG], T[BIG];

int main(int argc, char **argv)
{
	char buf[LEN_NAME];
	char tree_line[LEN_NAME];
	char name_line[LEN_NAME];
	struct p_tree *sp_tree = NULL;
	int root_id = 0;
	int i = 0, j = 0, k = 0;
	int num_sp = 0, num_lines = 0;
	struct mapping_list *ortho_maps;
	struct sp_list *genes;
	int num_genes = 0;
	FILE *f;
	int len = 0;
	int num_chars = 0;
	int max_gene_num = 0;
	int max_chars = 0;
	int gene_size = 0, line_size = 0, font_size = 0;
	int s_font_size = 0, g_font_size = 0, angle = 0;
	int ref_id = 0;
	
	count_node = 0;
	debug_mode = FALSE;
	strcpy(name_line, "");
	strcpy(tree_line, "");
	if( argc == 6 ) {
		if( strcmp(argv[4], "debug_mode") == 0 ) debug_mode = TRUE;
		else {
			fatal("args: ortho.view sp_tree parameters.txt\n");
		}
	}
	else if( argc != 5 ) {
		fatal("args: ortho.view sp_tree(in non-redundant.gc) param.txt ref_sp\n");
	}

	if((f = ckopen(argv[2], "r")) != NULL ) {
		if( fgets(buf, LEN_NAME, f) == NULL ) {
			fatalf("nothing in file %s\n", argv[2]);
		}
		else if( buf[0] == '(' ) {
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
				sp_tree = read_one_line(tree_line);
				root_id = sp_tree->nid;
				if( debug_mode == TRUE ) {
					print_tree(sp_tree, TREE_PRINT);
					printf("\n");
				}
			}
		}
		else {
			fatalf("not newick format: %s", buf);
		}
		fclose(f);
	}
	else {
		fatalf("file open error %s\n", argv[1]);
	}

	f = ckopen(argv[1], "r");
	num_lines = count_lines(f);
	num_lines--;
	ortho_maps = (struct mapping_list *) ckalloc(num_lines * sizeof(struct mapping_list));
	init_ortho_maps(ortho_maps, num_lines);
	fseek(f, 0, SEEK_SET);
	fgets(buf, LEN_NAME, f);
	num_chars = strlen(buf);
	num_genes = count_tokens(buf); // the first column is a row name
	num_genes--;
	max_gene_num = num_genes;
	genes = (struct sp_list *) ckalloc(num_genes * sizeof(struct sp_list));	
	len = strlen(buf);
	i = 0;
	j = 0;
	while( (j < num_genes) && (i < len) && (buf[i] != '\0') && (buf[i] != '\n') ) 
	{
		genes[j].id = j;
		if( j == 0 ) {
			i = concat_tokens(buf, i, genes[j].name); // gene_name
		}
		i = concat_tokens(buf, i, genes[j].name);
		j++;
	}
	fclose(f);
	read_view_file(argv[1], ortho_maps);

	max_chars = 0;
	for( i = 0; i < num_lines; i++ ) {
		if( ortho_maps[i].num_genes > max_gene_num ) {
			max_gene_num = ortho_maps[i].num_genes;
		}

		if( (int)strlen(ortho_maps[i].sp_name) > max_chars ) {
			max_chars = (int)strlen(ortho_maps[i].sp_name);
			if( strcmp(ortho_maps[i].sp_name, argv[4]) == 0 ) {
				ref_id = i;
			}
		}
	}
	num_chars = num_chars + max_chars;
	line_size = DPI * WIDTH;
	font_size = (int)((float)(1.5*(float)line_size)/(float)(num_chars+1));

	if( font_size >= SFS ) {
		s_font_size = SFS;
		g_font_size = GFS;
	}
	else {
		s_font_size = font_size+1;
		g_font_size = font_size;
	}

	if( (s_font_size <= MIN_SFS) || (g_font_size <= MIN_GFS) ) {
		s_font_size = MIN_SFS;
		g_font_size = MIN_GFS;
		angle = 45;
	}
	else {
		angle = 0;
	}

	gene_size = (int)((float)line_size/(float)(max_gene_num));
	if( (int)((float)gene_size * (float)2 / (float)3) >= GSIZE ) {
		gene_size = GSIZE;
	}
	else {
		gene_size = (int)((float)gene_size * (float)2 / (float)3);
	}

	f = ckopen(argv[3], "w");
	fprintf(f, "%d %d %d %d\n", gene_size, s_font_size, g_font_size, angle);
	fclose(f);

	assign_sp_id_mapping_list(sp_tree, ortho_maps, 0, num_lines);

	if( is_in_list(ref_id, sp_tree->ch_sp, sp_tree->num_csp) == true ) {
		flip_tree(sp_tree, ref_id);
		print_view_file(sp_tree, genes, num_genes, ortho_maps, num_lines);	
	}

	for( i = 0; i < num_lines; i++ ) {
		for( j = 0; j < ortho_maps[i].num_genes; j++ ) {
			for( k = 0; k < ortho_maps[i].o_id[j].num_parts; k++ ) {
				if( ortho_maps[i].o_id[j].pieces[k].num_ids > 0 ) {
					free(ortho_maps[i].o_id[j].pieces[k].ids);
				}
			}
			if( ortho_maps[i].o_id[j].num_parts > 0 ) {
				free(ortho_maps[i].o_id[j].pieces);
			}
		}
		if( ortho_maps[i].num_genes > 0 ) {
			free(ortho_maps[i].o_id);
		}
	}
	free(ortho_maps);
	free(genes);
	free_p_tree(sp_tree);
	return(EXIT_SUCCESS);
}	
