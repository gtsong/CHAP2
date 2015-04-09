#include "main.h"
#include "adjust_for_indels.h"
#include "util.h"
#include "util_i.h"

void adjust_for_indels(char op, struct I src, struct I dst, int len, int flag)
{
	char *ops;
	struct I *origin;
	struct I *dup;
	int num_ops = 0;
	char buf[500];
	char a;
	int b, c, d, e;
	int i, j;
	int *adj_len;
	struct I *temp_origin;
	struct I *temp_dup;
	FILE *outfile;
	FILE *temp_output; 
	 
	outfile = fopen("predicted_ops", "r");
	while(fgets(buf, 500, outfile) != NULL) num_ops++;

	ops = (char *) ckalloc(sizeof(char) * num_ops);
	origin = (struct I *) ckalloc(sizeof(struct I) * num_ops);
	dup = (struct I *) ckalloc(sizeof(struct I) * num_ops);
	temp_origin = (struct I *) ckalloc(sizeof(struct I) * num_ops);
	temp_dup = (struct I *) ckalloc(sizeof(struct I) * num_ops);
	adj_len = (int *) ckalloc(sizeof(int) * num_ops);

	fseek(outfile, 0, SEEK_SET);
	num_ops = 0;
	while(fgets(buf, 500, outfile) != NULL)
	{
		sscanf(buf, "%c %d %d %d %d", &a, &b, &c, &d, &e);
		ops[num_ops] = a;
		origin[num_ops] = assign_I(b, c);
		dup[num_ops] = assign_I(d, e);
		temp_origin[num_ops] = assign_I(b, c);
		temp_dup[num_ops] = assign_I(d, e);
	}
	fclose(outfile);

	ops[num_ops] = op;
	origin[num_ops] = assign_I(src.lower, src.upper);
	dup[num_ops] = assign_I(dst.lower, dst.upper);
	temp_origin[num_ops] = assign_I(src.lower, src.upper);
	temp_dup[num_ops] = assign_I(dst.lower, dst.upper);
	num_ops++;

	if(flag == DEL)
	{
		adj_len[num_ops-1] = len;
	}
	else if(flag == INS)
	{
		adj_len[num_ops-1] = (-1) * len;
	}

	for( i = (num_ops - 1); i >= 0; i-- )
	{
		if( temp_origin[i].lower >= temp_dup[i].lower )
		{
			temp_origin[i] = assign_I(temp_origin[i].lower + width(temp_dup[i]), temp_origin[i].upper + width(temp_dup[i]));
		}

		for( j = (num_ops - 1); j > i; j-- )
		{
			if( temp_dup[j].lower >= temp_dup[i].lower )  
			{
				temp_dup[j] = assign_I(temp_dup[j].lower + width(temp_dup[i]), temp_dup[j].upper + width(temp_dup[i]));	
			}
			
			if( temp_origin[j].lower >= temp_dup[i].lower )
			{
				temp_origin[j] = assign_I(temp_origin[j].lower + width(temp_dup[i]), temp_origin[j].upper + width(temp_dup[i]));
			}
		}
	}

	reflect_diff( origin, dup, temp_origin, temp_dup, num_ops );

	temp_output = fopen("predicted_ops", "w");
	for( i = 0; i < (num_ops-1); i++ )
	{
		fprintf(temp_output, "%c %d %d %d %d\n", ops[i], origin[i].lower, origin[i].upper, dup[i].lower, dup[i].upper);
	}
	fclose(temp_output);

	temp_output = fopen("temp.output", "w");
	for( i = 0; i < num_ops; i++ )
	{
		fprintf(temp_output, "%c %d %d %d %d\n", ops[i], temp_origin[i].lower, temp_origin[i].upper, temp_dup[i].lower, temp_dup[i].upper);
	}
	fclose(temp_output);

	free(ops);
	free(origin);
	free(dup);
	free(temp_origin);
	free(temp_dup);
	free(adj_len);
}

void reflect_diff(struct I *origin, struct I *dup, struct I *temp_origin, struct I *temp_dup, int num_ops)
{
	int i;
	int len;
	struct I dup_reg;

	dup_reg = assign_I((temp_dup[num_ops-1].lower), (temp_dup[num_ops-1].lower) + width(temp_origin[num_ops-1]));

	if( width(temp_origin[num_ops-1]) == width(temp_dup[num_ops-1]) ) {}
	else
	{
		len = width(temp_origin[num_ops-1]) - width(temp_dup[num_ops-1]);

		for( i = (num_ops - 2); i >= 0; i-- )
		{
			if( (temp_origin[i].lower) >= (temp_dup[num_ops-1].lower) )
			{
				temp_origin[i] = assign_I((temp_origin[i].lower) + len, (temp_origin[i].upper) + len);
				origin[i] = assign_I((origin[i].lower) + len, (origin[i].upper) + len);
			}
			else if( subset( temp_origin[i], dup_reg ) == true )
			{

			}
			else if( in( dup_reg.upper, temp_origin[i] ) == true )
			{
				temp_origin[i] = assign_I((temp_origin[i].lower) + len, (temp_origin[i].upper) + len);
				origin[i] = assign_I((origin[i].lower) + len, (origin[i].upper) + len);
			}

			if( (temp_dup[i].lower) >= (temp_dup[num_ops-1].lower) )
			{
				temp_dup[i] = assign_I((temp_dup[i].lower) + len, (temp_dup[i].upper) + len);
				dup[i] = assign_I((dup[i].lower) + len, (dup[i].upper) + len);
			}
		}
	}
}
