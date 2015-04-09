#ifndef CAL_PID_CONV_H
#define CAL_PID_CONV_H

void cal_pid_conv(struct DotList *algns, int num_algns, struct cv_list *cv, int num_conv, FILE *f);
int make_old_dup_list(struct DotList *algns, int num_algns, struct cv_list *init_dup_cv, int num_dup_cv, int *old_dups);

#endif /* CAL_PID_CONV_H */
