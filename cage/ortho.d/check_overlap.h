#ifndef CHECK_OVERLAP_H
#define CHECK_OVERLAP_H

bool find_overlapped_region(struct I *temp_ptr, struct I *pair_ptr, struct I cur_reg, int id, struct DotList *dots, int sign, bool is_x, int distance);
struct I change_value(struct I temp, int distance, int sign, bool is_x);

#endif /* CHECK_OVERLAP_H */
