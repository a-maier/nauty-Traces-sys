/* gtnauty.c :  nauty-related routines for gtools.

   Jan 15, 2001 : Increased graph order limit from 2^16-1 to 2^22-1.
   Aug  9, 2001 : Added fgroup_inv() and fcanonise_inv()
   Sep 15, 2004 : Completed prototypes
   Oct 16, 2004 : DEAFULTOPTIONS_GRAPH
   Nov 17, 2005 : Added fcanonise_inv_sg()
   May 11, 2010 : use sorttemplates.c
   Sep  5, 2013 : Unify format processing and remove 2^22 limit
   Oct 14, 2017 : Include code for n=0
   Sep 28, 2019 : Define breakcellwt
   Mar 26, 2024 : Increase workspace to match densenauty()
   Apr 13, 2024 : Make setlabptnfmt global and generalize.

**************************************************************************/

#include "gtools.h"   /* which includes naututil.h, nausparse.h, stdio.h */
#include "nautinv.h"

static boolean issymm;
static set *g0;
static int gm;
static int fuzz2[] = {006532,070236,035523,062437};
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3])

int gt_numorbits;

#ifdef REFINE
void REFINE(graph*,int*,int*,int,int*,int*,set*,int*,int,int);
#endif

#define MIN_SCHREIER 33  /* If n is this large, schreier will be
                            turned on. */

/**************************************************************************/

#define SORT_OF_SORT 3
#define SORT_NAME sortwt
#define SORT_TYPE1 int
#define SORT_TYPE2 int
#include "sorttemplates.c"
/* Creates static void sortwt(int *lab, int *wt, int n) */

void
setlabptn(int *weight, int *lab, int *ptn, int n)
/* Define (lab,ptn) with cells in increasing order of weight. */
{
    int i;

    if (n == 0) return;

    for (i = 0; i < n; ++i) lab[i] = i;

    if (weight)
    {
        sortwt(lab,weight,n);
        for (i = 0; i < n-1; ++i)
        {
            if (weight[lab[i]] != weight[lab[i+1]])
                ptn[i] = 0;
            else
                ptn[i] = 1;
        }
        ptn[n-1] = 0;
    }
    else
    {
        for (i = 0; i < n-1; ++i) ptn[i] = 1;
        ptn[n-1] = 0;
    }
}

int
breakcellwt(int *weight, int *lab, int *ptn, int n1, int n2)
/* Break (lab[n1..n2-1],ptn[n1..n2-1]) into cells in increasing
   order of weight.  If is assumed that lab[n1..n2-1] are defined
   but ptn[n1..n2-1] are ignored.
   The weight of lab[i] is weight[lab[i]]. 
   The number of cells is returned. */
{
    int i,nc;

    if (n2 <= n1) return 0;

    nc = 1;
    if (weight)
    {
        sortwt(lab+n1,weight,n2-n1);
        for (i = n1; i < n2-1; ++i)
        {
            if (weight[lab[i]] != weight[lab[i+1]])
            {
                ptn[i] = 0;
                ++nc;
            }
            else
                ptn[i] = 1;
        }
        ptn[n2-1] = 0;
    }
    else
    {
        for (i = n1; i < n2-1; ++i) ptn[i] = 1;
        ptn[n2-1] = 0;
    }

    return nc;
}

int
setlabptnfmt(char *fmt, int *lab, int *ptn, set *active, int m, int n)
/* Define (lab,ptn,active) according to format string.
   The format string is a string of characters assumed extended forever
   with 'z'. The shortcut c^N, where N is an integer, means c repeated
   N times.

   Grammar:  ( '-' | ) (<char> | <char> '^' <number>)^*

   The cells of the partition are defined by the character for each
   vertex in ascii order. If fmt starts with '-', vertices are counted
   backwards starting from vertices n-1 and reverse ascii order is used.

   active is not set if it is NULL.
   Return the number of cells.  */
{
    int i,j,mult,a,b,nc;
    unsigned char *s,*t;
    boolean minus;
#if MAXN
    int wt[MAXN];
#else
    DYNALLSTAT(int,wt,wt_sz);

    DYNALLOC1(int,wt,wt_sz,n,"setlabptnfmt");
#endif

    if (n == 0) return 0;

    if (active != NULL)
    {
        EMPTYSET(active,m);
        ADDELEMENT(active,0);
    }

    if (fmt != NULL && fmt[0] != '\0')
    {
#if !MAXN
        DYNALLOC1(int,wt,wt_sz,n,"setlabptnfmt");
#endif
        if (fmt[0] == '-')
        {
            minus = TRUE;
            s = (unsigned char*)(fmt+1);
        }
        else
        {
            minus = FALSE;
            s = (unsigned char*)fmt;
        }

        i = 0;
        while (i < n && *s != '\0')
        {
            if (*(s+1) == '^' && *(s+2) >= '0' && *(s+2) <= '9')
            {
                t = s+2;
                mult = 0;
                while (*t >= '0' && *t <= '9')
                {
                    mult = mult*10 + (*t - '0');
                    ++t;
                }
            }
            else
            {
                mult = 1;
                t = s+1;
            }
            for (j = 0; j < mult && i < n; ++i, ++j)
                wt[i] = (int)(*s);
            s = t;
        }

        for ( ; i < n; ++i) wt[i] = (int)('z');

        for (i = 0; i < n; ++i) lab[i] = i;
        if (minus)
        {
            for (i = 0, j = n-1; i <= j; ++i, --j)
            {
                a = wt[i];
                b = wt[j];
                wt[i] = -b;
                wt[j] = -a;
            }
        }

        sortwt(lab,wt,n);

        nc = 1;
        for (i = 0; i < n-1; ++i)
        {
            if (wt[lab[i]] != wt[lab[i+1]])
            {
                ptn[i] = 0;
                ++nc;
            }
            else
                ptn[i] = 1;
        }
        ptn[n-1] = 0;

        if (active != NULL)
        {
            for (i = 0; i < n-1; ++i)
                if (ptn[i] == 0) ADDELEMENT(active,i+1);
        }
    }
    else
    {
        for (i = 0; i < n; ++i)
        {
            lab[i] = i;
            ptn[i] = 1;
        }
        ptn[n-1] = 0;
        nc = 1;
    }

    return nc;
}

/**************************************************************************/

static boolean
hasloops(graph *g, int m, int n)
/* Test for loops */
{
    int i;
    set *gi;

    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i)) return TRUE;

    return FALSE;
}

static boolean
hasloops_sg(sparsegraph *sg)
{
    size_t *v,vi,j;
    int *d,*e,n,i;

    n = sg->nv;
    SG_VDE(sg,v,d,e);
    for (i = 0; i < n; ++i)
    {
        vi = v[i];
        for (j = vi; j < vi + d[i]; ++j)
            if (e[vi] == i) return TRUE;
    }

    return FALSE;
}

void
fcanonise(graph *g, int m, int n, graph *h, char *fmt, boolean digraph)
/*  canonise g under format fmt; result in h.
   fmt is either NULL (for no vertex classification) or is a string
   with char-valued colours for the vertices.  If it ends early, it
   is assumed to continue with the colour 'z' indefinitely. */
{
#if MAXN
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int count[MAXN];
    set active[MAXM];
    setword workspace[2*500*MAXM];
#else
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,count,count_sz);
    DYNALLSTAT(set,active,active_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
#endif
    int i;
    int numcells,code;
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);

    if (n == 0) return;

#if MAXN
    if (n > MAXN || m > MAXM)
    {
        fprintf(stderr,">E fcanonise: m or n too large\n");
        NAUTY_ABORT(">E fcanonise");
    }
#else
    DYNALLOC1(int,lab,lab_sz,n,"fcanonise");
    DYNALLOC1(int,ptn,ptn_sz,n,"fcanonise");
    DYNALLOC1(int,orbits,orbits_sz,n,"fcanonise");
    DYNALLOC1(int,count,count_sz,n,"fcanonise");
    DYNALLOC1(set,active,active_sz,m,"fcanonise");
    DYNALLOC1(setword,workspace,workspace_sz,2*500*m,"fcanonise");
#endif

    digraph = digraph || hasloops(g,m,n);

    numcells = setlabptnfmt(fmt,lab,ptn,active,m,n);

    if (m == 1)
        refine1(g,lab,ptn,0,&numcells,count,active,&code,1,n);
    else
        refine(g,lab,ptn,0,&numcells,count,active,&code,m,n);

    if (numcells == n || (numcells == n-1 && !digraph))
    {
        for (i = 0; i < n; ++i) count[i] = lab[i];
        updatecan(g,h,count,0,m,n);
        gt_numorbits = numcells;
    }
    else
    {
        options.getcanon = TRUE;
        options.defaultptn = FALSE;
        options.digraph = digraph;
#ifdef REFINE
        options.userrefproc = REFINE;
#endif
        if (n >= MIN_SCHREIER) options.schreier = TRUE;

        EMPTYSET(active,m);
        nauty(g,lab,ptn,active,orbits,&options,&stats,
                                              workspace,2*500*m,m,n,h);
        gt_numorbits = stats.numorbits;
    }
}

/**************************************************************************/

void
fcanonise_inv(graph *g, int m, int n, graph *h, char *fmt,
   void (*invarproc)(graph*,int*,int*,int,int,int,int*,int,
    boolean,int,int), int mininvarlevel, int maxinvarlevel,
    int invararg, boolean digraph)
/* Canonise g under format fmt; result in h.
   fmt is either NULL (for no vertex classification) or is a string
   with char-valued colours for the vertices.  If it ends early, it
   is assumed to continue with the colour 'z' indefinitely.
   This is like fcanonise() except that a invariant and its arguments
   can be specified. */
{
#if MAXN
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int count[MAXN];
    set active[MAXM];
    setword workspace[2*500*MAXM];
#else
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,count,count_sz);
    DYNALLSTAT(set,active,active_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
#endif
    int i;
    int numcells,code;
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);

    if (n == 0) return;

#if MAXN
    if (n > MAXN || m > MAXM)
    {
        fprintf(stderr,">E fcanonise: m or n too large\n");
        NAUTY_ABORT(">E fcanonise");
    }
#else
    DYNALLOC1(int,lab,lab_sz,n,"fcanonise");
    DYNALLOC1(int,ptn,ptn_sz,n,"fcanonise");
    DYNALLOC1(int,orbits,orbits_sz,n,"fcanonise");
    DYNALLOC1(int,count,count_sz,n,"fcanonise");
    DYNALLOC1(set,active,active_sz,m,"fcanonise");
    DYNALLOC1(setword,workspace,workspace_sz,2*500*m,"fcanonise");
#endif

    numcells = setlabptnfmt(fmt,lab,ptn,active,m,n);
    digraph = digraph || hasloops(g,m,n);

    if (m == 1)
        refine1(g,lab,ptn,0,&numcells,count,active,&code,1,n);
    else
        refine(g,lab,ptn,0,&numcells,count,active,&code,m,n);

    if (numcells == n || (!digraph && numcells >= n-1))
    {
        for (i = 0; i < n; ++i) count[i] = lab[i];
        updatecan(g,h,count,0,m,n);
        gt_numorbits = numcells;
    }
    else
    {
        options.getcanon = TRUE;
        options.digraph = digraph;
        options.defaultptn = FALSE;
        if (invarproc)
        {
            options.invarproc = invarproc;
            options.mininvarlevel = mininvarlevel;
            options.maxinvarlevel = maxinvarlevel;
            options.invararg = invararg;
        }
#ifdef REFINE
        options.userrefproc = REFINE;
#endif
        if (n >= MIN_SCHREIER) options.schreier = TRUE;

        EMPTYSET(active,m);
        nauty(g,lab,ptn,active,orbits,&options,&stats,workspace,2*500*m,m,n,h);
        gt_numorbits = stats.numorbits;
    }
}

/**************************************************************************/

void
fcanonise_inv_sg(sparsegraph *g, int m, int n, sparsegraph *h, char *fmt,
   void (*invarproc)(graph*,int*,int*,int,int,int,int*,int,
    boolean,int,int), int mininvarlevel, int maxinvarlevel,
    int invararg, boolean digraph)
/*  canonise g under format fmt; result in h.
   fmt is either NULL (for no vertex classification) or is a string
   with char-valued colours for the vertices.  If it ends early, it
   is assumed to continue with the colour 'z' indefinitely.
   This is like fcanonise() except that a invariant and its arguments
   can be specified.  Version for sparse graphs. */
{
#if MAXN
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    int count[MAXN];
    set active[MAXM];
    setword workspace[2*500*MAXM];
#else
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,count,count_sz);
    DYNALLSTAT(set,active,active_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
#endif
    int i;
    int numcells,code;
    statsblk stats;
    static DEFAULTOPTIONS_SPARSEGRAPH(options);

    if (n == 0)
    {
        h->nv = 0;
        h->nde = 0;
        return;
    }

#if MAXN
    if (n > MAXN || m > MAXM)
    {
        fprintf(stderr,">E fcanonise: m or n too large\n");
        NAUTY_ABORT(">E fcanonise");
    }
#else
    DYNALLOC1(int,lab,lab_sz,n,"fcanonise");
    DYNALLOC1(int,ptn,ptn_sz,n,"fcanonise");
    DYNALLOC1(int,orbits,orbits_sz,n,"fcanonise");
    DYNALLOC1(int,count,count_sz,n,"fcanonise");
    DYNALLOC1(set,active,active_sz,m,"fcanonise");
    DYNALLOC1(setword,workspace,workspace_sz,2*500*m,"fcanonise");
#endif

    numcells = setlabptnfmt(fmt,lab,ptn,active,m,n);
    digraph = digraph || hasloops_sg(g);

    refine_sg((graph*)g,lab,ptn,0,&numcells,count,active,&code,1,n);

    if (numcells == n || (!digraph && numcells == n-1))
    {
        for (i = 0; i < n; ++i) count[i] = lab[i];
        updatecan_sg((graph*)g,(graph*)h,count,0,m,n);
        gt_numorbits = numcells;
    }
    else
    {
        options.getcanon = TRUE;
        options.digraph = digraph;
        options.defaultptn = FALSE;
        if (invarproc)
        {
            options.invarproc = invarproc;
            options.mininvarlevel = mininvarlevel;
            options.maxinvarlevel = maxinvarlevel;
            options.invararg = invararg;
        }
#ifdef REFINE
        options.userrefproc = REFINE;
#endif
        if (n >= MIN_SCHREIER) options.schreier = TRUE;

        EMPTYSET(active,m);
        nauty((graph*)g,lab,ptn,active,orbits,&options,&stats,
                                         workspace,2*500*m,m,n,(graph*)h);
        gt_numorbits = stats.numorbits;
    }
}

/**************************************************************************/

void
fgroup(graph *g, int m, int n, char *fmt, int *orbits, int *numorbits)  
/* Find the orbits of undirected graph g stabilised by format fmt.
   The orbits are put into orbits[] and the number of them into *numorbits
   fmt is either NULL (for no vertex classification) or is a string
   with char-valued colours for the vertices.  If it ends early, it
   is assumed to continue with the colour 'z' indefinitely. */
{
#if MAXN
    int lab[MAXN],ptn[MAXN];
    int count[MAXN];
    set active[MAXM];
    setword workspace[2*500*MAXM];
#else
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,count,count_sz);
    DYNALLSTAT(set,active,active_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
#endif
    int i,j;
    int orbrep;
    int numcells,code;
    boolean digraph;
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);

    if (n == 0)
    {
        *numorbits = 0;
        return;
    }

#if MAXN
    if (n > MAXN || m > MAXM)
    {
        fprintf(stderr,">E fcanonise: m or n too large\n");
        NAUTY_ABORT(">E fcanonise");
    }
#else
    DYNALLOC1(int,lab,lab_sz,n,"fcanonise");
    DYNALLOC1(int,ptn,ptn_sz,n,"fcanonise");
    DYNALLOC1(int,count,count_sz,n,"fcanonise");
    DYNALLOC1(set,active,active_sz,m,"fcanonise");
    DYNALLOC1(setword,workspace,workspace_sz,2*500*m,"fcanonise");
#endif

    numcells = setlabptnfmt(fmt,lab,ptn,active,m,n);
    digraph = hasloops(g,m,n);

    if (m == 1)
        refine1(g,lab,ptn,0,&numcells,count,active,&code,1,n);
    else
        refine(g,lab,ptn,0,&numcells,count,active,&code,m,n);

    if (cheapautom(ptn,0,digraph,n))
    {
        for (i = 0; i < n; )
        {
            if (ptn[i] == 0)
            {
                orbits[lab[i]] = lab[i];
                ++i;
            }
            else
            {
                orbrep = n;
                j = i;
                do
                {
                    if (lab[j] < orbrep) orbrep = lab[j];
                } while (ptn[j++] != 0);

                for (; i < j; ++i) orbits[lab[i]] = orbrep;
            }
        }
        *numorbits = gt_numorbits = numcells;
    }
    else
    {
        options.getcanon = FALSE;
        options.defaultptn = FALSE;
        options.digraph = digraph;
#ifdef REFINE
        options.userrefproc = REFINE;
#endif
        if (n >= MIN_SCHREIER) options.schreier = TRUE;

        EMPTYSET(active,m);
        nauty(g,lab,ptn,active,orbits,&options,&stats,workspace,2*500*m,m,n,NULL);
        *numorbits = gt_numorbits = stats.numorbits;
    }
}

/**************************************************************************/

void
fgroup_inv(graph *g, int m, int n, char *fmt, int *orbits, int *numorbits,
      void (*invarproc)(graph*,int*,int*,int,int,int,int*,int,
       boolean,int,int), int mininvarlevel, int maxinvarlevel, int invararg)  
/* Find the orbits of undirected graph g stabilised by format fmt.
   The orbits are put into orbits[] and the number of them into *numorbits
   fmt is either NULL (for no vertex classification) or is a string
   with char-valued colours for the vertices.  If it ends early, it
   is assumed to continue with the colour 'z' indefinitely.
   This is like fgroup() except that a invariant and its arguments
   can be specified. */
{
#if MAXN
    int lab[MAXN],ptn[MAXN];
    int count[MAXN];
    set active[MAXM];
    setword workspace[2*500*MAXM];
#else
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,count,count_sz);
    DYNALLSTAT(set,active,active_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
#endif
    int i,j;
    int orbrep;
    boolean digraph;
    int numcells,code;
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);

    if (n == 0)
    {
        *numorbits = 0;
        return;
    }

#if MAXN
    if (n > MAXN || m > MAXM)
    {
        fprintf(stderr,">E fcanonise: m or n too large\n");
        NAUTY_ABORT(">E fcanonise");
    }
#else
    DYNALLOC1(int,lab,lab_sz,n,"fcanonise");
    DYNALLOC1(int,ptn,ptn_sz,n,"fcanonise");
    DYNALLOC1(int,count,count_sz,n,"fcanonise");
    DYNALLOC1(set,active,active_sz,m,"fcanonise");
    DYNALLOC1(setword,workspace,workspace_sz,2*500*m,"fcanonise");
#endif

    numcells = setlabptnfmt(fmt,lab,ptn,active,m,n);
    digraph = hasloops(g,m,n);

    if (m == 1)
        refine1(g,lab,ptn,0,&numcells,count,active,&code,1,n);
    else
        refine(g,lab,ptn,0,&numcells,count,active,&code,m,n);

    if (cheapautom(ptn,0,digraph,n))
    {
        for (i = 0; i < n; )
        {
            if (ptn[i] == 0)
            {
                orbits[lab[i]] = lab[i];
                ++i;
            }
            else
            {
                orbrep = n;
                j = i;
                do
                {
                    if (lab[j] < orbrep) orbrep = lab[j];
                } while (ptn[j++] != 0);

                for (; i < j; ++i) orbits[lab[i]] = orbrep;
            }
        }
        *numorbits = gt_numorbits = numcells;
    }
    else
    {
        options.getcanon = FALSE;
        options.defaultptn = FALSE;
        options.digraph = digraph;
        if (invarproc)
        {
            options.invarproc = invarproc;
            options.mininvarlevel = mininvarlevel;
            options.maxinvarlevel = maxinvarlevel;
            options.invararg = invararg;
        }
#ifdef REFINE
        options.userrefproc = REFINE;
#endif
        if (n >= MIN_SCHREIER) options.schreier = TRUE;

        EMPTYSET(active,m);
        nauty(g,lab,ptn,active,orbits,&options,&stats,workspace,2*500*m,m,n,NULL);
        *numorbits = gt_numorbits = stats.numorbits;
    }
}

/**************************************************************************/

static void
userlevel(int *lab, int *ptn, int level, int *orbits, statsblk *stats,
      int tv, int index, int tcellsize, int numcells, int cc, int n)
{
    int i0,i;

    if (level != 2) return;

    issymm = TRUE;

    i0 = nextelement(g0,gm,-1);
    if (i0 >= 0)
        for (i = i0; (i = nextelement(g0,gm,i)) >= 0;)
            if (orbits[i] != i0)
            {
                issymm = FALSE;
                return;
            }
}
 
/*******************************************************************/

/* istransitive(g,m,n,h)

   g   is an input undirected graph without loops
   m,n of standard meaning.  
   h   is a place to put an output graph.

   If g is transitive, return 1 or 2 and put a canonically labelled
       version of g into h.  The value is 2 for symmetric graphs, 
       and 1 for other transitive graphs.
   If g is not transitive, return 0.  In that case h may or 
       may not have something in it.
*/
int
istransitive(graph *g, int m, int n, graph *h)
{
    int i,inv;
    set *gw;
    short wt;
    int d,inv0,v,w;
    statsblk stats; 
    static DEFAULTOPTIONS_GRAPH(options);
#if MAXN
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    setword workspace[2*500*MAXM];
    set workset[MAXM];
    set sofar[MAXM],frontier[MAXM];
#else
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
    DYNALLSTAT(set,workset,workset_sz);
    DYNALLSTAT(set,sofar,sofar_sz);
    DYNALLSTAT(set,frontier,frontier_sz);
#endif

    if (n == 0) return 2;

#if MAXN
    if (m > MAXM || n > MAXN)
    {
        fprintf(stderr,
            ">E istransitive: bad input parameters (n=%d m=%d)\n",n,m);
        exit(1);
    }
#else
    DYNALLOC1(int,lab,lab_sz,n,"istransitive");
    DYNALLOC1(int,ptn,ptn_sz,n,"istransitive");
    DYNALLOC1(int,orbits,orbits_sz,n,"istransitive");
    DYNALLOC1(setword,workspace,workspace_sz,2*500*m,"istransitive");
    DYNALLOC1(set,workset,workset_sz,m,"istransitive");
    DYNALLOC1(set,sofar,sofar_sz,m,"istransitive");
    DYNALLOC1(set,frontier,frontier_sz,m,"istransitive");
#endif

    inv0 = 0;
    for (v = 0; v < n; ++v)
    {
        inv = 0;
        EMPTYSET(sofar,m);
        ADDELEMENT(sofar,v);
        EMPTYSET(frontier,m);
        ADDELEMENT(frontier,v);
        for (d = 1; d < n; ++d)
        {
            EMPTYSET(workset,m);
            wt = 0;
            for (w = -1; (w = nextelement(frontier,m,w)) >= 0;)
            {
                ++wt;
                gw = GRAPHROW(g,w,m);
                for (i = m; --i >= 0;) workset[i] |= gw[i];
            }
            if (wt == 0) break;
            wt += (short)(0x73 ^ d);
            wt = (short)FUZZ2(wt);
            inv += wt;
            for (i = m; --i >= 0;)
            {
                frontier[i] = workset[i] & ~sofar[i];
                sofar[i] |= frontier[i];
            }
        }
        if (v == 0) inv0 = inv;
        else if (inv != inv0) return 0;
    }

    options.getcanon = TRUE;
    options.userlevelproc = userlevel;
    if (hasloops(g,m,n)) options.digraph = TRUE;
#ifdef REFINE
    options.userrefproc = REFINE;
#endif
    if (n >= MIN_SCHREIER) options.schreier = TRUE;
 
    issymm = TRUE;
    g0 = (set*) g;
    gm = m;

    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,2*500*m,m,n,h);

    if (stats.numorbits != 1) return 0;
    else if (!issymm)         return 1;
    else                      return 2;
}

/**************************************************************************/

void 
tg_canonise(graph *g, graph *h, int m, int n)
/* Canonise vertex-transitive graph */
{
    int i;
#if MAXN
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    set active[MAXM];
    setword workspace[2*500*MAXM];
#else
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(set,active,active_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
#endif
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(options);

#if MAXN
    if (n > MAXN || m > MAXM)
    {
        fprintf(stderr,">E tg_canonise: m or n too large\n");
        NAUTY_ABORT(">E tg_canonise");
    }
#else
    DYNALLOC1(int,lab,lab_sz,n,"tg_canonise");
    DYNALLOC1(int,ptn,ptn_sz,n,"tg_canonise");
    DYNALLOC1(int,orbits,orbits_sz,n,"tg_canonise");
    DYNALLOC1(set,active,active_sz,m,"tg_canonise");
    DYNALLOC1(setword,workspace,workspace_sz,2*500*m,"tg_canonise");
#endif

    if (n == 0) return;

    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    if (hasloops(g,m,n)) options.digraph = TRUE;
#ifdef REFINE
    options.userrefproc = REFINE;
#endif

    for (i = 0; i < n; ++i)
    {
        lab[i] = i;
        ptn[i] = 1;
    }
    ptn[0] = ptn[n-1] = 0;

    active[0] = 0;
    EMPTYSET(active,m);
    ADDELEMENT(active,0);

    if (n >= MIN_SCHREIER) options.schreier = TRUE;
    nauty(g,lab,ptn,active,orbits,&options,&stats,workspace,2*500*m,m,n,h);
}

/***********************************************************************/

typedef struct arc_st { int from,to; } arc;

DYNALLSTAT(arc,arclist,arclist_sz);
DYNALLSTAT(size_t,arcorbits,arcorbits_sz);
static TLS_ATTR size_t numarcs,numarcorbits;
static TLS_ATTR graph *g_save;
static TLS_ATTR int m_save;

static size_t
findarc(arc *a, int na, int from, int to)
/* Find position of from->to.  (Must be present.) */
{
    size_t lo,mid,hi;

    lo = 0;
    hi = na - 1;

    while (lo <= hi)
    {
        mid = lo + (hi - lo) / 2;
        if (a[mid].from == from)
        {
            if (a[mid].to == to)
                return mid;
            else if (a[mid].to > to)
                hi = mid - 1;
            else
                lo = mid + 1;
        }
        else if (a[mid].from > from)
            hi = mid - 1;
        else
            lo = mid + 1;
    }
    gt_abort(">E findarc error\n");
}

void
arcorbitjoin(int ngens, int *p, int *orbs, int numorbs, int stab, int n)
/* p is a vertex permutation. Apply it to arcorbits. */
{
    size_t i,pi,j1,j2,t;
    int ii,jj,m;
    graph *gi;

    if (ngens == 1)
    {
        m = m_save;
        DYNALLOC1(arc,arclist,arclist_sz,numarcs,"countorbits");
        DYNALLOC1(size_t,arcorbits,arcorbits_sz,numarcs,"countorbits");

        t = 0;
        for (ii = 0, gi = g_save; ii < n; ++ii, gi += m)
        {
            for (jj = -1; (jj = nextelement(gi,m,jj)) >= 0; )
            {
                arclist[t].from = ii;
                arclist[t].to = jj;
                ++t;
            }
        }
        for (i = 0; i < numarcs; ++i) arcorbits[i] = i;

        numarcorbits = 0;
        for (i = 0; i < numarcs; ++i)
        if (arcorbits[i] == i)
        {
            ++numarcorbits;
            pi = i;
            do
            {
                pi = findarc(arclist,numarcs,
                        p[arclist[pi].from],p[arclist[pi].to]);
                arcorbits[pi] = i;
            } while (pi != i);
        }

        return;
    }

    for (i = 0; i < numarcs; ++i)
    {
        pi = findarc(arclist,numarcs,p[arclist[i].from],p[arclist[i].to]);

        if (pi != i)
        {
            j1 = arcorbits[i];
            while (arcorbits[j1] != j1) j1 = arcorbits[j1];
            j2 = arcorbits[pi];
            while (arcorbits[j2] != j2) j2 = arcorbits[j2];

            if (j1 < j2)      arcorbits[j2] = j1;
            else if (j1 > j2) arcorbits[j1] = j2;
        }
    }

    numarcorbits = 0;
    for (i = 0; i < numarcs; ++i)
        if ((arcorbits[i] = arcorbits[arcorbits[i]]) == i) ++numarcorbits;
}

static size_t
arcorbtoedgeorb(void)
/* Convert arc orbits to edge orbits; strictly undirected only */
{
    size_t i,pi,j1,j2,numedgeorbits;

    for (i = 0; i < numarcs; ++i)
    {
        if (arclist[i].from < arclist[i].to)
        {
            pi = findarc(arclist,numarcs,arclist[i].to,arclist[i].from);

            j1 = arcorbits[i];
            while (arcorbits[j1] != j1) j1 = arcorbits[j1];
            j2 = arcorbits[pi];
            while (arcorbits[j2] != j2) j2 = arcorbits[j2];

            if (j1 < j2)      arcorbits[j2] = j1;
            else if (j1 > j2) arcorbits[j1] = j2;
        }
    }

    numedgeorbits = 0;
    for (i = 0; i < numarcs; ++i)
        if ((arcorbits[i] = arcorbits[arcorbits[i]]) == i) ++numedgeorbits;

    return numedgeorbits;
}

void
countorbits_sg(sparsegraph *sg, boolean digraph,
        double *grpsize1, int *grpsize2,
        int *vorbits, int *fixedpts, size_t *eorbits, size_t *aorbits)
/* Find:
   (grpsize1,grpsize2) = group size as in the nauty options structure.
   vorbits = number of orbits on vertices
   fixedpts = number of fixed points on vertices
   eorbits = number of orbits on undirected edges
   aorbits = number of orbits on directed edges
   In the case of digraphs, eorbits is set equal to aorbits.
*/ 
{
    gt_abort(">E countorbits_sg is not implemented\n");
}

void
countorbits(graph *g, int m, int n, boolean digraph,
        double *grpsize1, int *grpsize2, int *vorbits,
        int *fixedpts, size_t *eorbits, size_t *aorbits)
/* Find:
   (grpsize1,grpsize2) = group size as in the nauty options structure.
   vorbits = number of orbits on vertices
   fixedpts = number of fixed points on vertices
   eorbits = number of orbits on undirected edges
   aorbits = number of orbits on directed edges
   In the case of digraphs, eorbits is set equal to aorbits.
*/ 
{
    graph *gi;
    size_t t;
    int loops,i,fixed;
    statsblk stats;
    static DEFAULTOPTIONS_GRAPH(goptions);
    static DEFAULTOPTIONS_DIGRAPH(doptions);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(setword,work,work_sz);

    numarcs = 0;
    for (t = 0; t < m*(size_t)n; ++t) numarcs += POPCOUNT(g[t]);

    if (numarcs == 0)
    {
        *grpsize1 = 1.0;
        *grpsize2 = 0;
        for (i = 2; i <= n; ++i) MULTIPLY(*grpsize1,*grpsize2,i);
        *vorbits = 1;
        *fixedpts = (n == 1 ? 1 : 0);
        *eorbits = 1;
        *aorbits = 1;
        return;
    }
    loops = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i)) ++loops;

    g_save = g;
    m_save = m;

    DYNALLOC1(int,lab,lab_sz,n,"countorbits");
    DYNALLOC1(int,ptn,ptn_sz,n,"countorbits");
    DYNALLOC1(int,orbits,orbits_sz,n,"countorbits");
    DYNALLOC1(setword,work,work_sz,2*500*m,"countorbits");

    if (digraph)
    {
        doptions.userautomproc = arcorbitjoin;
        nauty(g,lab,ptn,NULL,orbits,&doptions,&stats,work,2*500*m,m,n,NULL);
    }
    else
    {
        goptions.userautomproc = arcorbitjoin;
        if (loops > 0) goptions.digraph = TRUE; 
        nauty(g,lab,ptn,NULL,orbits,&goptions,&stats,work,2*500*m,m,n,NULL);
    }
 
    *grpsize1 = stats.grpsize1;
    *grpsize2 = stats.grpsize2;
    *vorbits = stats.numorbits;

    if (*vorbits == n)
    {
        *aorbits = numarcs;
        *eorbits = (digraph ? numarcs : (numarcs + loops) / 2);
    }
    else
    {
        *aorbits = numarcorbits;
        *eorbits = (digraph ? numarcorbits : arcorbtoedgeorb());
    }

    for (i = 0; i < n; ++i) ptn[i] = 0;
    fixed = stats.numorbits;
    for (i = 0; i < n; ++i)
        if (++ptn[orbits[i]] == 2) --fixed;
    *fixedpts = fixed;

    if (n > 128)
    {
        DYNFREE(lab,lab_sz);
        DYNFREE(ptn,ptn_sz);
        DYNFREE(orbits,orbits_sz);
        DYNFREE(work,work_sz);
        DYNFREE(arclist,arclist_sz);
        DYNFREE(arcorbits,arcorbits_sz);
    }
}
