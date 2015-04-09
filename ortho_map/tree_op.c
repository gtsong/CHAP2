#include "main.h"
#include "tree_op.h"
#include "util.h"
#include "util_i.h"

extern int count_node;

void init_tree(struct p_tree *t)
{
	t->left = NULL;
	t->right = NULL;
	t->parent = NULL;
	t->reg = assign_I(0,1);
	t->name = NULL;
	t->b_len = (double) 0;
  t->d_mode = SP; 
  t->od = 0; // orientation for printing orthologous alignments
  t->sp_code = -1; // a code number of self-alignment: species id, seq id for orthologous mappings
  t->gid = -1; // a gene identifier 
  t->nid = -1; // a node identifier
  t->val = 0;
  t->ch_sp = NULL; // the list of children nodes(species for a species tree)
  t->num_csp = 0; // the number of child species
  t->depth = 0; // the depth of each node
  t->visited = false;
}

bool is_sep_char(char ch) 
{
	if (ch == ',' || ch == '(' || ch == ')' || ch == ';' || ch == ':' || ch == '[' || ch == ']')
		return true;
	return false;
}

bool is_done_node(struct p_tree *node)
{
	if( (node->d_mode == DONE) || (node->d_mode == LEAF) )
		return true;
	return false;
}

bool is_leaf_node(struct p_tree *node)
{
	if( (node->left == NULL) && (node->right == NULL) )
		return true;
	return false;
}

struct p_tree *read_tree(char *filename) 
{
	FILE *tree_file = fopen(filename, "r");
	struct p_tree *root = NULL;
	char *str;
	char *one_line;

	str = (char *) ckalloc(MAX_LEN*sizeof(char));
	one_line = (char *) ckalloc(MAX_NEWICK*sizeof(char));
	str[0] = '\0';
	one_line[0] = '\0';
	while( fgets(str, MAX_LEN, tree_file) ) {
		strcat(one_line, str);
	}
	
	leave_only_taxons(one_line);
	count_node = 0;
	root = read_one_line(one_line); // in here, it is assumed that the last taxon is an out-group and the branch of the last taxon and others is a root. So re-rooting might be necessary.

	free(str);
	free(one_line);
	return root;
}

struct p_tree *read_one_line(char *one_line) {
	struct p_tree *p = NULL, *q = NULL, *t = NULL;
	struct p_tree *stack[STACK_SIZE];
	char buf[100];
	char *pt;
	int i, top = 0;
	double b_len = 1;
//	int count_gene = 0;
	int state;
	int temp_len;
	int temp_id = -1;

	pt = one_line;
	while(*pt != '\0' && *pt != ';') {
		if( !is_sep_char(*pt) ) {
			i = 0;
			while(!is_sep_char(*pt)) {
				buf[i++] = *pt;
				pt++;
			}
			buf[i] = '\0';
//			count_gene++;
			state = NAME;

			p = (struct p_tree *) ckalloc(sizeof(struct p_tree));	
			init_tree(p);
			p->nid = count_node;
			count_node++;
			temp_len = strlen(buf) + 1;
			p->name = (char *) ckalloc(temp_len * sizeof(char));
			strcpy(p->name, buf);
			buf[0] = '\0';
//			p->gid = count_gene;
			p->left = NULL;
			p->right = NULL;
		}

		switch (*pt) 
		{
			case ':': 
				pt++;
				if (sscanf(pt, "%lf", &b_len) != 1) {
					printf("branch length error %s\n", pt);
					exit(EXIT_FAILURE);
				}
				while (!is_sep_char(*pt))
					pt++;
				break;
				
			case '(': 
				p = (struct p_tree *) ckalloc(sizeof(struct p_tree));
				init_tree(p);
				p->nid = count_node;
				count_node++;
				p->left = NULL;
				p->right = NULL;
				stack[top++] = p;
				pt++;
				state = LPAR;
				break;
			
			case '[':
				pt++;
				i = 0;
				while(!is_sep_char(*pt)) {
					buf[i++] = *pt;
					pt++;
				}
				buf[i] = '\0';
				temp_id = atoi(buf);
				buf[0] = '\0';
				state = LBR;
				break;

			case ']':
				pt++;
				state = RBR;
				break;

			case ',': 
				p->b_len = b_len;
				q = stack[top-1];
				if( q->left == NULL ) {
					if( temp_id != -1 ) {
						p->sp_code = temp_id;
						temp_id = -1;
					}
					q->left = p;
				}
				else if( q->right == NULL ) {
					p->sp_code = temp_id;
					temp_id = -1;
					q->right = p;
				}
				else {
					t = (struct p_tree *) ckalloc(sizeof(struct p_tree));
					init_tree(t);
					t->nid = count_node;
					count_node++;
					if( temp_id != -1 ) {
						q->sp_code = temp_id;
						temp_id = -1;
					}
					t->left = q;
					t->right = NULL;
					stack[top-1] = t;
				}

				q = stack[top-1];
				p->parent = q;
				pt++;
				state = COMMA;
				break;
				
			case ')': 
				if( top == 0 ) {
					printf("tree description unbalanced\n");
					exit(EXIT_FAILURE);
				}
				else
				{
					q = stack[--top];
					p->parent = q;
					p->b_len = b_len;
					if( q->right == NULL ) 
					{
						if( temp_id != -1 ) {
							p->sp_code = temp_id;
							temp_id = -1;
						}
						q->right = p;
						p = q;
					}
					else {
						t = (struct p_tree *) ckalloc(sizeof(struct p_tree));
						init_tree(t);
						t->nid = count_node;
						count_node++;
						if( temp_id != -1 ) {
							p->sp_code = temp_id;
							temp_id = -1;
						}
						t->left = q;
						t->right = p;
						stack[top-1] = t;
						p = t;
					}
					pt++;
				}
				state = RPAR;
				break;
				
			case ';':
				break;

			default:
				printf("undefined symbol\n");
				exit(EXIT_FAILURE);
				break;
		}
	}

	p->b_len = 1;
	if( temp_id != -1 ) {
		p->sp_code = temp_id;
	}
	return p;
}

void leave_only_taxons(char *one_line)
{
	int i;
	int len;
	int j = 0;

	len = strlen(one_line);

	for( i = 0; i < len; i++ )
	{
		if( !isspace(one_line[i]) ) {
			one_line[j] = one_line[i];
			j++;
		}
	}

	if( one_line[j-1] != ';' ) {
		one_line[j] = ';';
		one_line[j+1] = '\0';
	}
	else one_line[j]='\0';
}

void free_p_tree(struct p_tree *t)
{ 
  if(t == NULL) return;
  else
  { 
    if( t->left != NULL) free_p_tree(t->left);
    if( t->right != NULL) free_p_tree(t->right);
		if( t->ch_sp != NULL ) free(t->ch_sp);
		if( t->name != NULL ) free(t->name);
    free(t); 
  }
} 

void postorder(struct p_tree *t)
{
	if( t->left != NULL ) postorder(t->left);
	if( t->right != NULL ) postorder(t->right);
	if(t->d_mode == DP) printf("%d - %f : duplication %d | ", t->num_csp, t->b_len, t->depth);
	else if(t->d_mode == SP) printf("%d - %f : speciation %d | ", t->num_csp, t->b_len, t->depth);
	else if(t->d_mode == DP_W_LOSS) printf("%d - %f : dup_loss %d | ", t->num_csp, t->b_len, t->depth);
	else if(t->d_mode == SP_W_LOSS) printf("%d - %f : sp_loss %d | ", t->num_csp, t->b_len, t->depth);
	else if(t->d_mode == CV) printf("%d - %f : conversion %d | ", t->num_csp, t->b_len, t->depth);
	else if(t->d_mode == DONE) printf("%d - %f : done %s %d | ", t->num_csp, t->b_len, t->name, t->depth);
	else printf("%d - %f : %s %d | ", t->ch_sp[0], t->b_len, t->name, t->depth);
}

void mapping_tree(struct p_tree *gt, struct p_tree *sp)
{
	if( gt->left != NULL && gt->d_mode != DONE ) mapping_tree(gt->left, sp);
	if( gt->right != NULL && gt->d_mode != DONE ) mapping_tree(gt->right, sp);

	if( gt->d_mode == DONE ) {}
	else if( is_leaf_node(gt) == true ) {
		gt->d_mode = LEAF;
		gt->visited = false;
	}
	else { 
		gt->d_mode = assign_dmode(gt, sp);
		gt->visited = false;
	}
}

int assign_dmode(struct p_tree *gt, struct p_tree *sp)
{
	int l_nid, r_nid; // a node id of a left child and a right child
	int s_list[NUM_SP];	
	int num_sp;
	int i;
	struct p_tree *t;
	struct p_tree *l, *r;

	t = gt->left;
	num_sp = t->num_csp; 
	for( i = 0; i < num_sp; i++ ) {
		s_list[i] = t->ch_sp[i];
	}
	l = r = NULL;
	l_nid = r_nid = -1;
	l = find_lca(s_list, num_sp, sp, &l_nid);

	t = gt->right;
	num_sp = t->num_csp; 
	for( i = 0; i < num_sp; i++ ) {
		s_list[i] = t->ch_sp[i];
	}
	r = find_lca(s_list, num_sp, sp, &r_nid);

	if( l_nid == -1 || r_nid == -1 ) {
		printf("a gene includes uncovered species\n");
		exit(EXIT_FAILURE);
	}
	else if( l == r ) {
		return(DP);
	}
	else if( is_descendant(l_nid, r_nid, sp) == true || is_descendant(r_nid, l_nid, sp) == true ) {
		return(DP_W_LOSS);
	}
	else if( l->parent == r->parent ) {
		return(SP);
	}
	else {
		return(SP_W_LOSS);
	}
}

struct p_tree *find_lca(int *list, int num_sp, struct p_tree *sp, int *nid)
{
	int *s_list = NULL;
	int cur_num_sp = 0;
	int i = 0;
	struct p_tree *lca = NULL;

	cur_num_sp = sp->num_csp;
	s_list = (int *) ckalloc(sizeof(int) * cur_num_sp);
	for(i = 0; i < cur_num_sp; i++ ) {
		s_list[i] = sp->ch_sp[i];
	}

	if( is_covering(list, num_sp, s_list, cur_num_sp) == true ) {
		if( sp->left != NULL) lca = find_lca(list, num_sp, sp->left, nid);
		else {
			lca = NULL;
		}

		if( lca == NULL ) {
			if( sp->right != NULL ) lca = find_lca(list, num_sp, sp->right, nid);
			else lca = NULL;
		}

		free(s_list);
		if( lca == NULL ) {
			*nid = sp->nid;
			return(sp);
		}
		else return(lca);
	}
	else {
		free(s_list);
		return(NULL);
	}
}

bool is_covering(int *list, int num_sp, int *cmp, int num_cmp)
{
	int i = 0;

	if( num_sp > num_cmp ) return(false);
	else
	{
		while((i < num_sp) && (is_in_list(list[i], cmp, num_cmp) == true)) i++; 
		if( i >= num_sp ) return(true);
		else return(false);
	}
}

bool is_descendant(int nid1, int nid2, struct p_tree *t) // is nid1 a descendant of nid2 in tree t?
{
	bool res = false;

	if( t->nid == nid2 ) {
		if( t->left != NULL ) res = is_node_in_tree(nid1, t->left);

		if( res == false ) {
			if( t->right != NULL ) res = is_node_in_tree(nid1, t->right);
		}
	}	
	else { 
		if( t->left != NULL ) res = is_descendant(nid1, nid2, t->left);
		
		if( res == false ) {
			if( t->right != NULL ) res = is_descendant(nid1, nid2, t->right);
		}
	}
	return(res);
}

bool is_node_in_tree(int nid, struct p_tree *t)
{
	bool res = false;

	if( t->nid == nid ) return(true);
	else {
		if( t->left != NULL ) res = is_node_in_tree(nid, t->left);
		else res = false;

		if( res == false ) {
			if( t->right != NULL ) res = is_node_in_tree(nid, t->right);
			else res = false;
		}
		return(res);
	}
}

void assign_sp_id(struct p_tree *t, struct sp_list *list, int depth, int num_sp)
{
	struct p_tree *temp = NULL;
	int count;
	int i;

	if( t->left != NULL ) assign_sp_id(t->left, list, depth+1, num_sp);
	if( t->right != NULL ) assign_sp_id(t->right, list, depth+1, num_sp);

	if( is_leaf_node(t) || is_done_node(t) ) {
			t->depth = depth;
			t->sp_code = get_id_in_list(t->name, list, num_sp);				
			t->num_csp = 1;
			t->ch_sp = (int *) ckalloc(sizeof(int));
			t->ch_sp[0] = t->sp_code; // the initialization of a species list
//			if( debug_mode == true)	printf("%d: %s %d\n", t->nid, t->name, t->sp_code);
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

void assign_gid(struct p_tree *t, struct g_list *genes, int num_genes)
{
	int i;

	if( t->left != NULL ) assign_gid(t->left, genes, num_genes);
	if( t->right != NULL ) assign_gid(t->right, genes, num_genes);

	if( is_leaf_node(t) == true ) {
		i = 0;
		while( (i < num_genes) && (strcmp(t->name, genes[i].gname) != 0) ) i++;

		if( i >= num_genes ) {
			printf("%s is not in the gene list\n", t->name);
			exit(EXIT_FAILURE);
		}
		else t->gid = genes[i].gid;	
	}
	else {
		t->gid = -1;
	}
}

int get_id_in_list(char *name, struct sp_list *list, int num_sp)
{
	int i;

	i = 0;
	while( ( i < num_sp) && (strcmp(name, list[i].name) != 0) ) i++;

	if( i >= num_sp ) {
		fatalf("error: no species name exist in %s\n", name);	
		return(num_sp);
	}
	else return(list[i].id);
}

bool is_in_list(int id, int *list, int num_id)
{
	int i = 0;

	while( (i < num_id) && (id != list[i]) ) i++;

	if( i >= num_id ) return(false);
	else return(true);
}

void get_deepest_sp(struct p_tree *sp, int *sp1, int *sp2, int *out_sp, int *cur_max) // sp1 is a left node and sp2 is a right node
{
	struct p_tree *temp;
	struct p_tree *par;

	if( sp->left != NULL ) get_deepest_sp(sp->left, sp1, sp2, out_sp, cur_max);
	if( sp->right != NULL ) get_deepest_sp(sp->right, sp1, sp2, out_sp, cur_max);

	if( is_leaf_node(sp) == false ) {
		if( is_leaf_node(sp->left) && is_leaf_node(sp->right) ) {
			if( (*cur_max) < sp->depth ) {
				*cur_max = sp->depth;
				temp = sp->left;
				*sp1 = temp->sp_code;
				temp = sp->right;
				*sp2 = temp->sp_code;	
				par = sp->parent;
				temp = par->left;
				if( temp->nid == sp->nid ) *out_sp = get_leaf_node(par->right, LEFT);
				else *out_sp = get_leaf_node(par->left, RIGHT);
			}
		}
	}
}

int get_leaf_node(struct p_tree *sp, int dir)
{
	struct p_tree *temp;

	temp = sp;
	while( (!is_leaf_node(temp)) && (temp != NULL) ) {
		if( dir == RIGHT ) temp = temp->right;
		else if( dir == LEFT ) temp = temp->left; 
	}
	
	if( temp == NULL ) return(-1);		
	else return(temp->sp_code);
}

struct p_tree *get_lca_node(struct p_tree *cp, int sp_id)
{
	struct p_tree *t;
	t = cp;
	while( !(is_in_list(sp_id, t->ch_sp, t->num_csp)) && (t->d_mode != ROOT) ) t = t->parent;
	return(t);
}

void swap_nodes(struct p_tree *a, struct p_tree *b)
{
	struct p_tree *t = NULL;

	t = (struct p_tree *) ckalloc(sizeof(struct p_tree));
	init_tree(t);
	assign_node(t, a, NO_DEALLOC);
	t->parent = a->parent;
	t->depth = a->depth;
	assign_node(a, b, NO_DEALLOC);
	assign_node(b, t, NO_DEALLOC);

	free(t);
}

void make_new_node(struct p_tree *ch_l, struct p_tree *ch_r)
{
	struct p_tree *node = NULL;
	int i;

	node = (struct p_tree *) ckalloc(sizeof(struct p_tree));
	init_tree(node);
	ch_r->depth = ch_r->depth + 1;
	assign_node(node, ch_l, NO_DEALLOC);
	node->depth = ch_l->depth + 1;
	node->nid = count_node;
	count_node++;
	node->parent = ch_l;
	ch_r->parent = ch_l;

	ch_l->left = node;
	ch_l->right = ch_r;
	ch_l->d_mode = SP;

	for( i = 0; i < (ch_r->num_csp); i++ ) {
		if( !is_in_list(ch_r->ch_sp[i], ch_l->ch_sp, ch_l->num_csp) ) {
			ch_l->ch_sp[ch_l->num_csp] = ch_r->ch_sp[i];
			ch_l->num_csp = (ch_l->num_csp) + 1;
		}
	}
	strcpy(ch_l->name, "\0");
}

void mark_visited(struct p_tree *c_node)
{
	c_node->visited = true;
	if( c_node->left != NULL ) mark_visited(c_node->left);
	if( c_node->right != NULL ) mark_visited(c_node->right);
}

void print_tree(struct p_tree *t, int mode)
{
	if( !is_leaf_node(t) )
	{
		if( mode == TREE_PRINT ) printf("(");
		if( !is_leaf_node(t) ) print_tree(t->left, mode);
		if( mode == TREE_PRINT ) printf(",");
		if( !is_leaf_node(t) ) print_tree(t->right, mode);
		if( mode == TREE_PRINT ) printf(")");
	}
//	else printf("%s", t->name);
	else {
		if( mode == TREE_PRINT ) {
			if( t->sp_code != -1 ) printf("%s%d", t->name, t->sp_code);
			else printf("%s", t->name);
		}
		else if( mode == INT_PRINT ) {
			if( t->sp_code != -1 ) {
				if( t->od == 0 ) printf("s %s%d + %d %d\n", t->name, t->sp_code, t->reg.lower, (t->reg.upper)-1); 
				else if( t->od == 1 ) printf("s %s%d - %d %d\n", t->name, t->sp_code, t->reg.lower, (t->reg.upper)-1); 
			}
			else {
				if( t->od == 0 ) printf("s %s + %d %d\n", t->name, t->reg.lower, (t->reg.upper)-1); 
				else if( t->od == 1 ) printf("s %s - %d %d\n", t->name, t->reg.lower, (t->reg.upper)-1); 
			}
		}
	}
}

void delete_one_sp(struct p_tree *t, int sp)
{
	struct p_tree *left, *right;
	struct p_tree *parent;	

	left = t->left;
	right = t->right;
	parent = t->parent;

	if( (left != NULL) && (is_leaf_node(left)) && left->sp_code == sp ) {
		if( (right != NULL) && (is_leaf_node(right)) && right->sp_code == sp ) {
			if( t->left->ch_sp != NULL ) free(t->left->ch_sp);
			free(t->left);
			if( t->right->ch_sp != NULL ) free(t->right->ch_sp);
			free(t->right);
			left = NULL;
			right = NULL;
		}
		else {
			if( t->left->ch_sp != NULL ) free(t->left->ch_sp);
			free(t->left);
			if( t->right != NULL ) {
				assign_node(t, t->right, DEALLOC);
				left = t->left;
				right = t->right;
			}
		}
	}
	else if( (right != NULL) && (is_leaf_node(right)) && right->sp_code == sp ) {
		if( t->right->ch_sp != NULL ) free(t->right->ch_sp);
		free(t->right);
		if( t->left != NULL ) {
			assign_node(t, t->left, DEALLOC);
			left = t->left;
			right = t->right;
		}
	}

	if( (left != NULL) && !is_leaf_node(left) ) delete_one_sp(left, sp);
	if( (right != NULL) && !is_leaf_node(right) ) delete_one_sp(right, sp);

	if( t->d_mode != LEAF && left == NULL && right == NULL ) {
		if( parent->left == t ) assign_node(parent, parent->right, DEALLOC);
		else if( parent->right == t) assign_node(parent, parent->left, DEALLOC);
		if( t->num_csp > 0 ) free(t->ch_sp);
		free(t);
	}
}

void assign_node(struct p_tree *a, struct p_tree *b, int mode)
{
	int i;

	if( b->name != NULL ) {
		a->name = (char *) ckalloc(sizeof(char)*(strlen(b->name)+1));
		strcpy(a->name, b->name);
		if( mode == DEALLOC ) free(b->name);
	}		
	a->d_mode = b->d_mode;
	a->od = b->od;
	a->sp_code = b->sp_code;
	a->gid = b->gid;
	a->val = b->val;
	a->num_csp = b->num_csp;
	for( i = 0; i < a->num_csp; i++ ) a->ch_sp[i] = b->ch_sp[i];
	a->visited = b->visited;
	a->left = b->left;
	a->right = b->right;
	if( a->left != NULL ) a->left->parent = a;
	if( a->right != NULL ) a->right->parent = a;
	if( mode == DEALLOC ) free(b);
}

void write_bk_gname(struct p_tree *gt, struct g_list *genes)
{
	if( gt->left != NULL ) write_bk_gname(gt->left, genes);
	if( gt->right != NULL ) write_bk_gname(gt->right, genes);

	if( is_leaf_node(gt) ) {
		strcpy(gt->name, genes[gt->gid].gname);
	}
}

void extract_gtree_3sp(struct p_tree *t, int sp1, int sp2, int out_sp, FILE *f)
{
	if( (!is_leaf_node(t)) && ( is_3sp_in_list(sp1, sp2, out_sp, t->left) && is_3sp_in_list(sp1, sp2, out_sp, t->right) ) )
	{
		fprintf(f, "(");
		if( !is_leaf_node(t) ) extract_gtree_3sp(t->left, sp1, sp2, out_sp, f);
		fprintf(f, ",");
		if( !is_leaf_node(t) ) extract_gtree_3sp(t->right, sp1, sp2, out_sp, f);
		fprintf(f, ")");
	}
	else if( (!is_leaf_node(t)) && ( is_3sp_in_list(sp1, sp2, out_sp, t->left) ) )
	{
		extract_gtree_3sp(t->left, sp1, sp2, out_sp, f);
	}
	else if( (!is_leaf_node(t)) && ( is_3sp_in_list(sp1, sp2, out_sp, t->right) ) )
	{
		extract_gtree_3sp(t->right, sp1, sp2, out_sp, f);
	}
	else if( is_leaf_node(t) ) fprintf(f, "%s", t->name);
	else {}
}

bool is_3sp_in_list(int sp1, int sp2, int out_sp, struct p_tree *t)
{
	if( is_in_list(sp1, t->ch_sp, t->num_csp) || is_in_list(sp2, t->ch_sp, t->num_csp) || is_in_list(out_sp, t->ch_sp, t->num_csp) ) return(true);
	else return(false);
}

struct p_tree *add_node(struct p_tree *cur, int id, int num_list, char *name, int sp_code, int sign) {
	struct p_tree *t;

	t = (struct p_tree *) ckalloc(sizeof(struct p_tree));
	init_tree(t);
	t->left = NULL;
	t->right = NULL;
	if( cur != NULL ) t->parent = cur;
	else t->parent = NULL;
	t->ch_sp = (int *) ckalloc(sizeof(int) * num_list);
	t->name = (char *) ckalloc(sizeof(char) * LEN_NAME);
	strcpy(t->name, name);
	t->ch_sp[0] = id;
	t->sp_code = sp_code;
	t->od = sign;
	t->num_csp = 1;
	t->gid = id;
	return(t);
}

bool is_outgroup(struct p_tree *p, int sp1_code, int sp2_code)
{
	struct p_tree *node;	
	int list[2];
	int nid;
	int cur_num_sp;
	int *s_list;
	bool res = false;
	int i = 0;

	list[0] = sp1_code;
	list[1] = sp2_code;
	node = find_lca(list, 2, p, &nid);
	cur_num_sp = node->num_csp;
	s_list = (int *) ckalloc(sizeof(int) * cur_num_sp);
	for(i = 0; i < cur_num_sp; i++ ) {
		s_list[i] = node->ch_sp[i];
	}	
		
	res = is_covering(list, 2, s_list, cur_num_sp); 
	
	free(s_list);
	return(res);
}

bool is_within_branch(int out_code, int sp1_code, int sp2_code, struct p_tree *p)
{
	struct p_tree *node;
	int list[2];
	int cur_num_sp = 0;
	int i = 0;
	int cur_id = 0;
	bool res = false;
	int nid = 0;

	list[0] = sp1_code;
	list[1] = sp2_code;
	node = find_lca(list, 2, p, &nid);
	cur_num_sp = node->num_csp;
	for(i = 0; i < cur_num_sp; i++ ) {
		cur_id = node->ch_sp[i];
		if( (cur_id != sp1_code) && (cur_id != sp2_code) && (cur_id  == out_code )) 
		{
			res = true;	
		}
	}	

	return(res);
}

int assign_sp_code(char *line, struct sp_list *sp_code, int num_sp_code)
{
	char *pt;
	char buf[100];
	int i = 0;
	int num_sp = 0;

	pt = line;
	while(((*pt) != '\0') && ((*pt) != ';')) {
		buf[0] = '\0';
		if( !is_sep_char(*pt) ) {
			i = 0;
			while(!is_sep_char(*pt)) {
				buf[i++] = *pt;
				pt++;
			}
			buf[i] = '\0';
		}

		if( (buf[0] != '\0') && ((*pt) == '[') ) {
			strcpy(sp_code[num_sp].name, buf);
			pt++;
			i = 0;
			while(!is_sep_char(*pt)) {
				buf[i++] = *pt;
				pt++;
			}
			buf[i] = '\0';
			sp_code[num_sp].id = atoi(buf);
			num_sp++;
			if( num_sp > num_sp_code ) {
				fatalf("overflow: number of species %d vs. %d\n", num_sp, num_sp_code);
			}
			buf[0] = '\0';
		}
		pt++;
	}
	return(num_sp);
}
