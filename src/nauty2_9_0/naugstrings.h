/* naugstrings.h : Write graph6 or sparse6 strings into array. */
/* Version 1.1, Jun 2015. */

#include "gtools.h"

#ifdef __cplusplus
extern "C" {
#endif

void gtog6string(graph*,char**,int,int);
void gtos6string(graph*,char**,int,int);
void gtod6string(graph*,char**,int,int);
void sgtos6string(sparsegraph*,char**);
void sgtog6string(sparsegraph*,char**);
void sgtod6string(sparsegraph*,char**);
void gtois6string(graph*,graph*,char**,int,int);

#ifdef __cplusplus
}
#endif
