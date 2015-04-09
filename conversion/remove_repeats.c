/* remove_repeats - remove paralogous pairs which contain tandem repeats. */

#include <stdio.h>
#include <stdlib.h>
#include "maf.h"

static int REPEAT_MAX_DISTANCE=300;
static int REPEAT_MIN_COVERAGE=80; // 80%

struct rAli {
  struct mafAli* ali;
  int b1, e1, b2, e2;
  char repeat;
};

int getAliCount(struct mafAli* head) {
  int AliCount=0;

  for (; head!=NULL; head=head->next)
    AliCount++;
  return AliCount;
}

int compar_ali_start(void* a, void* b) {
	if( ( (*((struct mafAli**)a))->components->start - (*((struct mafAli**)b))->components->start ) != 0 ) 
	{
  	return ( (*((struct mafAli**)a))->components->start - (*((struct mafAli**)b))->components->start );
	}
	else if( ( (*((struct mafAli**)a))->components->next->start - (*((struct mafAli**)b))->components->next->start ) != 0 ) 
	{
		return( (*((struct mafAli**)a))->components->next->start - (*((struct mafAli**)b))->components->next->start );
	}	
	else if( ((*((struct mafAli**)a))->components->size - (*((struct mafAli**)b))->components->size ) != 0 )
	{
  	return ( (*((struct mafAli**)a))->components->size - (*((struct mafAli**)b))->components->size );
	}
	else {
  	return ( (*((struct mafAli**)a))->components->next->size - (*((struct mafAli**)b))->components->next->size );
	}
}

void remove_same_alignment(struct mafAli** phead) {
  struct mafAli** aliArr, *ali, *head=*phead;
  int count, i, j;
  
  count = getAliCount(head);
  aliArr = (struct mafAli**)malloc(count*sizeof(struct mafAli*));
  for (i=0, ali=head; i<count; i++, ali=ali->next)
    aliArr[i] = ali;
  for (i=0; i<count; i++)
    aliArr[i]->next = NULL;
  
  qsort((void*)aliArr, count, sizeof(struct mafAli*), (void *)compar_ali_start);
  
  for (i=1; i<count; i++)
    for (j=i-1; j>=0; j--) {
      if ( aliArr[j] == NULL )
	continue;
      if ( aliArr[j]->components->start != aliArr[i]->components->start )
	break;
      if ( aliArr[i]->components->next->start == aliArr[j]->components->next->start 
	   && aliArr[i]->components->size == aliArr[j]->components->size
	   && aliArr[i]->components->next->size == aliArr[j]->components->next->size ) {
	mafAliFree(&aliArr[i]);
	aliArr[i] = NULL;
	break;
      }
    }
  head = NULL;
  for (i=count-1; i>=0; i--)
    if ( aliArr[i] != NULL ) {
      aliArr[i]->next = head;
      head = aliArr[i];
    }
  free(aliArr);
  *phead = head;
}


struct rAli* aliList2rAliArr(struct mafAli* head, int count) {
  struct rAli* rAliArr;
  struct mafComp* comp;
  int i;

  if ( count == 0 )
    return NULL;
  rAliArr = (struct rAli*)malloc(count*sizeof(struct rAli));
  for (i=0; head!=NULL; head=head->next, i++) {
    rAliArr[i].ali = head;
    comp = head->components;
    rAliArr[i].b1  = comp->start;
    rAliArr[i].e1  = comp->start + comp->size - 1;
    comp = comp->next;
    rAliArr[i].b2  = comp->start; 
    rAliArr[i].e2  = comp->start + comp->size - 1;
    rAliArr[i].repeat = 'n';
  }
  for (i=0; i<count; i++)
    rAliArr[i].ali->next = NULL;
  return rAliArr;
}

void free_rAliArr(struct rAli* rAliArr) {
  free(rAliArr);
}

// 4 poinst: 1, 2, 3, 4
int repeat(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
  int overbeg, overend, overlen, midx, midy1, midy2, left, right;

  if ( x2 < x3 || x4 < x1 ) // no overlap on x
    return 0;
  if ( y2 < y3 || y4 < y1 ) // no overlap on y
    return 0;
  
  overbeg = ( x1 > x3 ? x1 : x3 );
  overend = ( x2 < x4 ? x2 : x4 );
  overlen = overend - overbeg + 1;
  overlen = overlen*100;

  if ( x2 == x1 || x4 == x3 )
    return 1;
  midx = (overbeg + overend)/2;
  midy1 = (y2 - y1)*(midx - x1)/(x2 - x1) + y1;
  midy2 = (y4 - y3)*(midx - x3)/(x4 - x3) + y3;

  if ( midy2 - midy1 > REPEAT_MAX_DISTANCE || midy1 - midy2 > REPEAT_MAX_DISTANCE )
    return 0;

  left  = overlen > REPEAT_MIN_COVERAGE*(x2 - x1 + 1);
  right = overlen > REPEAT_MIN_COVERAGE*(x4 - x3 + 1);
  
  if ( left && right ) 
    return 3;
  if ( left ) // the first line is contained in second
    return 1;
  if ( right ) // the second line is contained in first
    return 2;
  return 0;
}

int compar_rAli(void* a, void* b) {
	if( ( (*((struct rAli**)a))->ali->components->start - (*((struct rAli**)b))->ali->components->start ) != 0 ) 
	{
  	return ( (*((struct rAli**)a))->ali->components->start - (*((struct rAli**)b))->ali->components->start );
	}
	else if( ( (*((struct rAli**)a))->ali->components->next->start - (*((struct rAli**)b))->ali->components->next->start ) != 0 ) 
	{
		return( (*((struct rAli**)a))->ali->components->next->start - (*((struct rAli**)b))->ali->components->next->start );
	}	
	else if( ((*((struct rAli**)a))->ali->components->size - (*((struct rAli**)b))->ali->components->size ) != 0 )
	{
  	return ( (*((struct rAli**)a))->ali->components->size - (*((struct rAli**)b))->ali->components->size );
	}
	else {
  	return ( (*((struct rAli**)a))->ali->components->next->size - (*((struct rAli**)b))->ali->components->next->size );
	}
}

void remove_repeat(struct rAli* rAliArr, int AliCount) {
  struct rAli **pArr;
  int *rightMost, i, j, ret, done;

  pArr = (struct rAli**)malloc(AliCount*sizeof(struct rAli*));
  rightMost = (int*)malloc(AliCount*sizeof(int));
  for (i=0; i<AliCount; i++)
    pArr[i] = &(rAliArr[i]);
  
  qsort((void*)pArr, AliCount, sizeof(struct rAli*), (void *)compar_rAli);

  for (i=1, rightMost[0]=pArr[0]->e1; i<AliCount; i++)
    rightMost[i] =  pArr[i]->e1 > rightMost[i-1] ? pArr[i]->e1 : rightMost[i-1];
  
  for (i=0; i<AliCount; i++) {
    if ( pArr[i]->repeat == 'y' )
      continue;
    done = 0;
    for (j=i+1; j<AliCount && pArr[j]->b1 <= pArr[i]->e1; j++) {
      ret = repeat(pArr[i]->b1, pArr[i]->b2, pArr[i]->e1, pArr[i]->e2, 
		   pArr[j]->b1, pArr[j]->b2, pArr[j]->e1, pArr[j]->e2);
      if ( ret == 3 || ret == 1 )
	pArr[i]->repeat = 'y';
      if ( ret == 3 || ret == 2 )
	pArr[j]->repeat = 'y';
      if ( ret == 1 || ret == 3 ) {
	done = 1;
	break;
      }
      ret = repeat(pArr[i]->b2, pArr[i]->b1, pArr[i]->e2, pArr[i]->e1,
		   pArr[j]->b2, pArr[j]->b1, pArr[j]->e2, pArr[j]->e1);
      if ( ret == 3 || ret == 1 )
	pArr[i]->repeat = 'y';
      if ( ret == 3 || ret == 2 )
	pArr[j]->repeat = 'y';
      if ( ret == 1 || ret == 3 ) {
	done = 1;
	break;
      }
    }
    if ( done == 1 )
      continue;
    for (j=i-1; j>=0 && rightMost[j] >= pArr[i]->b1; j--) {
      ret = repeat(pArr[i]->b1, pArr[i]->b2, pArr[i]->e1, pArr[i]->e2, 
		   pArr[j]->b1, pArr[j]->b2, pArr[j]->e1, pArr[j]->e2);
      if ( ret == 3 || ret == 1 )
	pArr[i]->repeat = 'y';
      if ( ret == 3 || ret == 2 )
	pArr[j]->repeat = 'y';
      if ( ret == 1 || ret == 3 )
	break;
      ret = repeat(pArr[i]->b2, pArr[i]->b1, pArr[i]->e2, pArr[i]->e1,
		   pArr[j]->b2, pArr[j]->b1, pArr[j]->e2, pArr[j]->e1);
      if ( ret == 3 || ret == 1 )
	pArr[i]->repeat = 'y';
      if ( ret == 3 || ret == 2 )
	pArr[j]->repeat = 'y';
      if ( ret == 1 || ret == 3 ) 
	break;
    }
  }
  free(rightMost);
  free(pArr);
}

struct mafAli* rAliArr2aliList(struct rAli* rAliArr, int count) {
  struct mafAli* head=NULL;
  int i;

  for (i=count-1; i>=0; i--) {
    if ( rAliArr[i].repeat == 'n') {
      rAliArr[i].ali->next = head;
      head = rAliArr[i].ali;
    }
    else
      mafAliFree(&(rAliArr[i].ali));
    rAliArr[i].ali = NULL;
  }
  free(rAliArr);
  return head;
}

void remove_repeats(struct mafAli** pHead) {
  struct rAli* rAliArr;
  int count;

  count = getAliCount(*pHead);
  if ( count <= 0)
    return;
  rAliArr = aliList2rAliArr(*pHead, count);
  remove_repeat(rAliArr, count);
  *pHead = rAliArr2aliList(rAliArr, count);
}

#define REPEAT_TEST
#ifdef REPEAT_TEST
int main(int argc, char** argv) {
  struct mafFile* maf;
  struct mafAli* head;
  
  if ( argc != 2) {
	printf("remove_repeats MAF-file\n");
	return 1;
  }

  mafWriteStart(stdout, "norepeat");
  maf = mafReadAll(argv[1], 1);
  head = maf->alignments;
  maf->alignments = NULL;
  mafFileFree(&maf);

  remove_same_alignment(&head);
  remove_repeats(&head);

  for (; head!=NULL; head=head->next)
    mafWrite(stdout, head);

  return 0;
}
#endif
