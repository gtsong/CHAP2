/* check_species_tree - check species tree vs. file of sequence names */
#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "maf.h"
#include "multi_util.h"

struct TreeNode {
	struct TreeNode *left, *right, *parent;
	char speciesname[100];
};

struct TreeNode *stree;
struct TreeNode *species_nodes[100];
int species_num = 0;
char species_name[100][100];

char* gettoken(char *s, int *n) {
	char *t;

    t = s + *n;
    if(isalnum(s[*n])) {
		while(isalnum(s[*n]) || s[*n] == '_' || s[*n] == '-') 
			(*n)++;
//		if(s[*n] == ':')    // might not be 8 characters; better to
//			*n += 8;       //  prohibit branch lengths for CHAP
		return t;
    }
    if((s[*n] == '(')  || (s[*n] == ')' ) || (s[*n] == ',')) {
		(*n)++;
//		if(s[*n] == ':')
//			(*n) += 8;
		return t;
    }
	else {
		printf("Invalid character '%c' in species tree\n", s[*n]);
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
		strcpy(species_name[species_num], c->speciesname);
		species_num++;
			
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
		printf("Can't open species tree\n");
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

int find_species(char species[100]) {
	int i;
	
	for(i=0; i<species_num; i++) {
		if(strcmp(species_name[i], species) == 0) {
			strcpy(species_name[i], "");    // mark this one as used
			return i;
		}
	}

	printf("Can't find species \"%s\" in the species tree\n", species);
	exit(-1);
}

void check_species(char filename[100]) {
	FILE *s_file;
	char species[100];
	int i;

	s_file = fopen(filename, "r");
	while(fscanf(s_file, "%s", species) == 1)
		find_species(species);
	
	fclose(s_file);

	for(i=0; i<species_num; i++) {
		if(strcmp(species_name[i], "") != 0) {
			printf("Missing sequence for species \"%s\"\n", species_name[i]);
			exit(-1);
		}
	}

	return;
}

int main(int argc, char *argv[]) { 
	if ( argc != 3) {
		printf("Usage: check_species_tree species_tree species_file\n");
		return 1;
	}

	extract_species_tree(argv[1]);
	check_species(argv[2]);
	
	return 0;
}
