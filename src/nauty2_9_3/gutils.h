/* gutils.h - procedure declarations for gutil1.c and gutil2.c */

#ifndef  _GUTILS_H_    /* only process this file once */
#define  _GUTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

void degstats(graph*,int,int,
                     unsigned long*,int*,int*,int*,int*,boolean*);
void degstats2(graph*,boolean,int,int,unsigned long*,int*l,
     int*,int*,int*,int*, int*, int*,int*,int*, boolean*);
void degstats3(graph*,int,int,
                     unsigned long*,int*,int*,int*,int*,int*);
void diamstats(graph*,int,int,int*,int*);
void find_dist(graph*,int,int,int,int*);
void find_dist2(graph*,int,int,int,int,int*);
int numcomponents(graph*,int,int);
int numcomponents1(graph*,int);
int girth(graph*,int,int);
boolean isbiconnected1(graph*,int);
boolean isbiconnected(graph*,int,int);
boolean isbipartite(graph*,int,int);
int bipartiteside(graph*,int,int);
boolean twocolouring(graph*,int*,int,int);
boolean isconnected1(graph*,int);
boolean isconnected(graph*,int,int);
boolean issubconnected(graph*,set*,int,int); 
long maxcliques(graph*,int,int);
int maxcliquesize(graph*,int,int);
int maxindsetsize(graph*,int,int);
void sources_sinks(graph*,int,int,int*,int*);

long digoncount(graph*,int,int);
int loopcount(graph*,int,int);
long pathcount1(graph*,int,setword,setword);
long cyclecount1(graph*,int);
long cyclecount1lim(graph*,long,int);
long cyclecount(graph*,int,int);
long cyclecountlim(graph*,long,int,int);
long indpathcount1(graph*,int,setword,setword);
long indcyclecount1(graph*,int);
long indcyclecount(graph*,int,int);
void commonnbrs(graph*,int*,int*,int*,int*,int,int);
void contract1(graph*,graph*,int,int,int);
int cstarcontent(graph*,int,int);
long numtriangles1(graph*,int);
long numtriangles(graph*,int,int);
long numtriangles1(graph*,int);
long numind3sets(graph*,int,int);
long numind3sets1(graph*,int);
long numdirtriangles(graph*,int,int);
long numdirtriangles1(graph*,int);
long numsquares(graph*,int,int);
long numdiamonds(graph*,int,int);
long numpentagons(graph*,int,int);
long numhexagons1(graph*,int);
long numhexagons(graph*,int,int);
void delete1(graph*,graph*,int,int);
int conncontent(graph*,int,int);
boolean stronglyconnected(graph*,int,int);
int ktreeness1(graph*,int);
int ktreeness(graph*,int,int);

/* extern int diameter_sg(sparsegraph*,int*,int*); */

#ifdef __cplusplus
}
#endif

#endif
