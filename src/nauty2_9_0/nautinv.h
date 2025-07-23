/*****************************************************************************
* This is the header file for versions 2.7 of nautinv.c.                     *
*                                                                            *
*   Copyright (1984-) Brendan McKay.                                         *
*   Subject to the conditions and disclaimers in the file COPYRIGHT.         *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       20-Apr-01 : initial creation out of naututil.h                       *
*       10-Nov-10 : remove types shortish and permutation                    *
*                                                                            *
*****************************************************************************/

#include "nauty.h"              /* which includes stdio.h */

#ifdef __cplusplus
extern "C" {
#endif

void adjacencies(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void adjtriang(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void cellcliq(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void cellfano(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void cellfano2(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void cellind(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void cellquads(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void cellquins(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void celltrips(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void cellstarts(int*,int,set*,int,int);
void cliques(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void distances(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void getbigcells(int*,int,int,int*,int*,int*,int);
void indsets(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void nautinv_check(int,int,int,int);
void nautinv_freedyn(void);
void quadruples(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void refinvar(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void setnbhd(graph*,int,int,set*,set*);
void triples(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
void twopaths(graph*,int*,int*,int,int,int,int*,int,boolean,int,int);
#ifdef __cplusplus
}
#endif
