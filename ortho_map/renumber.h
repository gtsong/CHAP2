#ifndef RENUMBER_H
#define RENUMBER_H

struct mapping_list {
	int id;
	char sp_name[100];
	int num_genes;
	struct o_list *o_id;
};

struct o_list {
	int id;
	int num_parts;
	bool is_pseudo;
	struct part_list *pieces; 
};

struct part_list {
	int num_ids;
	int *ids; // a negative value indicates a pseudogene and its absolute value is a human orthologous gene id
};

#endif /* RENUMBER_H */
