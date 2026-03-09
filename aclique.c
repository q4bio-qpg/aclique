// NB: b14 looks like it's a better speed vs iterations tradeoff.
//     b14 looks to find max sooner given the same iteration limit.

// TODO: switch to a bit-vector for clique
// TODO: implement filter_nodes, good for sparse data
// TODO: node aging.  Track when we last tried a node, to ensure more even
//       random distribution when selecting.
// TODO: targetted mutation.  Replacing one node with another compatible node
//       as it may have a better improve metric.  Maybe not useful as we just
//       "remove+improve" now, which is equivalent to a swap anyway?

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

static int time_limit = INT_MAX; // in seconds

#ifndef MAX_CLIQUE
#  define MAX_CLIQUE 10000
#endif

//#define BAIL_OUT
#ifndef BO_TGROWTH    // Max no. of iterations without finding bigger clique
#  define BO_TGROWTH  50000
#endif
#ifndef BO_TMEMBER    // Max consecutive rep of prior-found clique && TGROWTH
#  define BO_TMEMBER  1000
#endif
#ifndef BO_TMEMBER2   // Max rep of prior-found clique, regardless of TGROWTH
#  define BO_TMEMBER2 50000
#endif

#ifndef MAX_BEST
#  define MAX_BEST 1024
#endif

// A small benefit
#define USE_NODE_FREQ

//#define POPULATION 256
//#define NITER 10000

#define MAX_POPULATION 1024

#ifndef POPULATION
#    define POPULATION 96
#endif

#ifndef NITER
#    define NITER 10000
#endif

#ifndef AGGRATE
#    define AGGRATE 0.2
#endif

//#define POPULATION 512
//#define NITER 100000

#define roundup(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#ifndef MIN
#    define MIN(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef MAX
#    define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#define HSIZE 22
#define HMASK ((1<<HSIZE)-1)
static uint8_t chash[1<<HSIZE] = {0};
static uint8_t bhash[1<<HSIZE] = {0};

typedef struct {
    int i1, i2;
} ac_pair;

int pair_cmp(const void *vp1, const void *vp2) {
    const ac_pair *p1 = (const ac_pair *)vp1;
    const ac_pair *p2 = (const ac_pair *)vp2;
    return p2->i2 - p1->i2;
}

int pair_cmp_rev(const void *vp1, const void *vp2) {
    const ac_pair *p1 = (const ac_pair *)vp1;
    const ac_pair *p2 = (const ac_pair *)vp2;
    return p1->i2 - p2->i2;
}

static int timestamp = 0;
typedef struct {
    char **conn;
    int *weight;
    int nnodes, nfiltered;
    int max_possible;
    int *filtered;
    int *degree;
    int random_sz;
    int *node_map;  // for picking nodes biased to high degree
    int *node_map2; // for picking nodes biased to low degree
    int degree_order[MAX_CLIQUE];
    int *neighbour_count;
    int **neighbour;
    // Probabilities of the node being useful
    int *node_used;
    int *node_unused;
    int *node_date;
} graph;

typedef struct {
    int size;
    int is_clique;
    char used[MAX_CLIQUE];
    int32_t hash;
} clique;

void random_clique(graph *g, clique *q, int size);
void improve(graph *g, clique *q, int cycles, int aggressive);

void reset_clique(clique *q) {
    memset(q->used, 1, MAX_CLIQUE);
    q->size = 0;
    q->is_clique = 0;
    q->hash = 0;
}

static clique pop[MAX_POPULATION];

uint32_t fnv1a_32(const void *data, size_t len) {
    const uint8_t *bytes = (const uint8_t *)data;
    uint32_t hash = 0x811C9DC5u;           // FNV offset basis
    const uint32_t prime = 0x01000193u;    // FNV prime

    for (size_t i = 0; i < len; i++) {
        hash ^= bytes[i];
        hash *= prime;
    }
    return hash;
}


void remove_chash(graph *g, clique *q) {
    q->hash = fnv1a_32(q->used, g->nnodes);
    //printf("Remove hash %u\n", q->hash);
    chash[q->hash & HMASK] = 0;
}

void add_chash(graph *g, clique *q) {
    q->hash = fnv1a_32(q->used, g->nnodes);
    //printf("Add hash %u\n", q->hash);
    chash[q->hash & HMASK]++;
}

// Returns 0 if the clique is new, otherwise 1.
int check_chash(graph *g, clique *q) {
    // DIMACS_subset_ascii/gen400_p0.9_55.clq, 256, 10000
    // 0:  43 43 46 46 47 47 47 49 49 49 (try 1)
    // 0:  46 46 47 47 47 48 48 48 49 48 (try 2)
    // 1:  46 46 46 47 47 47 48 48 48 48
    // 2:  46 46 46 46 46 47 47 48 48 48
    // 5:  46 46 47 47 47 47 47 47 48 49
    // 99: 45 47 47 47 47 48 48 48 48 49

    // pop=32,iter=20000
    // 0:  45 46 48 48 48 49 49 49 49 51
    // 2:  46 46 46 48 48 48 49 49 49 49 

    // pop=16,iter=20000
    // 0:  45 46 46 46 47 47 47 47 48 48 

    // pop=4,iter=50000 (approx same time)
    // 0:  45 46 46 48 48 48 48 49 49 50
    // => maybe we're not "learning" or evolving at all.  It's just a
    // numbers game with random start points!

    // latest, pop 32, iter 20000
    // 0:  49 49 50 50 50 50 50 50 51 51

    // Allow some dups as it can be helpful to have multiple good copies?
    return chash[q->hash & HMASK] > 0;
}

int check_bhash(graph *g, clique *q) {
    return bhash[q->hash & HMASK] > 0;
}

// Find the number of connections to a node
int degree(graph *g, int n) {
    int d = 0;
    for (int i = 0; i < g->nnodes; i++) {
	if (i != n)
	    d += g->conn[n][i];
    }
    return d;
}

// Failed idea given we start from Cherry's sub-graph.
// The notion is any linked pair of graphs with less than lower "bound"
// shared neighbours cannot be part of a clique together.  Therefore their
// edge can be removed and degree reduced.  Repeat.
void lower_bound_prune(graph *g, int bound) {
    int n = g->nnodes;

    int lowest_shared = n;
    int nd = n/50+.99;
    for (int i = 0; i < n; i++) {
	if (i%nd == nd-1)
	    putchar('.');fflush(stdout);
	for (int j = i+1; j < n; j++) {
	    if (!g->conn[i][j])
		continue;
//#define PRUNE3  // much too slow
#ifdef PRUNE3
	    for (int k = j+1; k < n; k++) {
		if (!g->conn[i][k] || !g->conn[j][k])
		    continue;
		int shared = 2;
#else
		int shared = 1;
#endif
		for (int z = 0; z < n; z++)
#ifdef PRUNE3
		    if (g->conn[i][z] && g->conn[j][z] && g->conn[k][z])
#else
		    if (g->conn[i][z] && g->conn[j][z])
#endif
			shared++;

		if (lowest_shared > shared)
		    lowest_shared = shared;

//		printf("%d(%d) %d(%d) %d(%d) shared degree %d\n",
//		       i, g->degree[i],
//		       j, g->degree[j],
//		       k, g->degree[k],
//		       shared);
#ifdef PRUNE3
	    }
#endif
	}
    }
    printf("\nLowest shared %d\n", lowest_shared);

    exit(0);
}

void init_graph_degree(graph *g) {
    int n = g->nnodes;
    g->degree = (int *)malloc(n * sizeof(*g->degree));
    int highest_degree = 0;
    int lowest_degree = g->nnodes;
    for (int i = 0; i < n; i++) {
	g->degree[i] = degree(g, i);
	if (highest_degree < g->degree[i])
	    highest_degree = g->degree[i];
	if (lowest_degree  > g->degree[i])
	    lowest_degree  = g->degree[i];
    }

    printf("Highest degree: %d\n", highest_degree);
    printf("Lowest  degree: %d\n", lowest_degree);

//    lower_bound_prune(g, 257);

    // Remove leading diagonal.  It's put in by compat.cpp, but we don't
    // work well with it.
    for (int i = 0; i < g->nnodes; i++)
	g->conn[i][i] = 0;

    // Reduce degree by the degree of things it connects to.
    // A node with degree 100 has edges to a node with degree 90
    // means this node can only be in a clique of maximum 99, not 100.
    // The actual clique size is likely less, but we can use a smaller
    // upper limit when optimising the order in which to add nodes.
    for (int i = 0; i < n; i++) {
	int d = g->degree[i];
	for (int j = 0; j < n; j++) {
	    if (g->conn[i][j] && d > g->degree[j])
		d--;
	}
	g->degree[i] = d;
    }

    // Sort by degree
    ac_pair p[MAX_CLIQUE];
    int total_degree = 0;
// Experimental, but pow of 1 (ie not) seems better.
#define PPP 1.0
    for (int i = 0; i < n; i++) {
	p[i].i1 = i;
	p[i].i2 = g->degree[i];
	//total_degree += p[i].i2;
	total_degree += pow(p[i].i2, PPP);
    }
    printf("Total degree = %d\n", total_degree);

    total_degree++;
    g->random_sz = total_degree;
    roundup(g->random_sz);
    double scale = (double)g->random_sz / (total_degree+g->nnodes);
    g->node_map = (int *)malloc(g->random_sz * sizeof(g->node_map));
    g->node_map2 = (int *)malloc(g->random_sz * sizeof(g->node_map2));

    qsort(p, n, sizeof(*p), pair_cmp);
    for (int i = 0; i < n; i++)
	g->degree_order[i] = p[i].i1;

    // map random value to node number such that high degrees are first.
    int k, i;
    for (k = i = 0; i < g->random_sz && k < n; k++) {
	for (int j = 0;
	     j < pow(p[k].i2,PPP) * scale && i < g->random_sz;
	     j++, i++) {
	    g->node_map[i] = p[k].i1;
	}
    }
    assert(k==n);
    printf("i=%d random_sz=%d total_degreepow=%d\n", i, g->random_sz, total_degree);
    while (i < g->random_sz)
	g->node_map[i++] = p[0].i1;
//    for (int i = 0; i < g->random_sz; i++) {
//	printf("%d\t%d\t%d\n", i, g->node_map[i], g->degree[g->node_map[i]]);
//    }


    // Inverse likelihood, with low degrees more frequently found
    qsort(p, n, sizeof(*p), pair_cmp_rev);
    for (k = i = 0; i < g->random_sz && k < n; k++) {
	for (int j = 0;
	     j < pow(p[k].i2,PPP) * scale && i < g->random_sz;
	     j++, i++) {
	    g->node_map2[i] = p[k].i1;
	}
    }
    assert(k==n);
    while (i < g->random_sz)
	g->node_map2[i++] = p[0].i1;
}

void init_graph_neighbours(graph *g) {
    int n = g->nnodes;
    g->neighbour_count=(int *)calloc(g->nnodes, sizeof(int));
    g->neighbour=(int **)calloc(g->nnodes, sizeof(int *));

    ac_pair p[MAX_CLIQUE];

    for (int i = 0; i < n; i++) {
	int c = 0;
	for (int j = 0; j < n; j++) {
	    if (g->conn[i][j])
		c++;
	}
	g->neighbour_count[i] = c;
	g->neighbour[i] = (int *)calloc(c, sizeof(int));

#if 1
	for (int j = 0, c = 0; j < n; j++) {
	    if (g->conn[i][j])
		g->neighbour[i][c++] = j;
	}
#else
	c = 0;
	for (int j = 0; j < n; j++) {
	    if (g->conn[i][j]) {
		p[c].i1 = j;
		p[c].i2 = g->degree_order[j];
		c++;
	    }
	}
	qsort(p, c, sizeof(*p), pair_cmp);
	for (int j = 0; i < c; i++)
	    g->neighbour[i][j] = p[j].i1;
#endif

	
    }
}

graph *load_graph(char *fn) {
    FILE *fp = fopen(fn, "r");
    int total_degree = 0;
    if (!fp) {
	perror(fn);
	return NULL;
    }

    graph *g = (graph *)malloc(sizeof(*g));
    int n, e;
    for (;;) {
	char line[1024];
	if (!fgets(line, sizeof(line), fp))
	    goto err;
	if (sscanf(line, "p edge %d %d\n", &n, &e) == 2)
	    break;
	if (sscanf(line, "p col %d %d\n", &n, &e) == 2)
	    break;
	if (sscanf(line, "p %d %d\n", &n, &e) == 2)
	    break;
    }
    printf("Nodes:  %d\n", n);
    printf("Edgess: %d\n", e);

    g->nnodes = n;
    g->max_possible = n;
    g->filtered = (int *)calloc(n, sizeof(*g->filtered));
    g->conn = (char **)calloc(n, sizeof(*g->conn));
    g->weight = (int *)malloc(n * sizeof(*g->weight));
    g->node_used = (int *)calloc(n, sizeof(int));
    g->node_unused = (int *)calloc(n, sizeof(int));
    g->node_date = (int *)calloc(n, sizeof(int));
    for (int i = 0; i < n; i++) {
	g->conn[i] = (char *)calloc(n, 1);
	g->weight[i] = 1;
	//g->conn[i][i] = 1; // why does this not work?
    }

    int n1, n2, d, w;
    while ((d = fscanf(fp, "v %d %d\n", &n1, &w) >= 2)) {
	n1--;
	if (n1 < 0 || n1 >= n) {
	    fprintf(stderr, "Node number out of range: %d > %d\n", n1, n);
	    goto err;
	}
	g->weight[n1] = w;
    }
    printf("d=%d\n", d);
    
    while ((d = fscanf(fp, "e %d %d", &n1, &n2)) >= 2) {
	int c;
	do {
	    c = fgetc(fp);
	} while (c != '\n' && c != EOF);

	n1--;
	n2--;
	if (n1 < 0 || n1 >= n ||
	    n2 < 0 || n2 >= n) {
	    fprintf(stderr, "Node number out of range: %d, %d > %d\n", n1, n2, n);
	    goto err;
	}

	g->conn[n1][n2] = 1;
	g->conn[n2][n1] = 1;
    }

    init_graph_degree(g);
    init_graph_neighbours(g);
    
    fclose(fp);
    return g;

 err:
    free(g);
    fclose(fp);
    return NULL;
}

#ifdef USE_NODE_FREQ
void update_node_stats(graph *g, clique *q) {
    // Doesn't help
    for (int i = 0; i < g->nnodes; i++) {
	if (q->used[i])
	    g->node_used[i]++;
	else
	    g->node_unused[i]++;
	
	if (g->node_used[i] + g->node_unused[i] > 1024) {
	    g->node_used[i] >>= 1;
	    g->node_unused[i] >>= 1;
	}
    }
}

int random_neighbour(graph *g, int node) {
    return g->neighbour_count[node]
	? g->neighbour[node][rand() % g->neighbour_count[node]]
	: 0;
}

int random_node(graph *g) {
    int n;
    for (int i = 0; i < 100; i++) {
	const double adj = 100; // min pick chance of 10%.
	n = g->node_map[random() % g->random_sz];
	// used/(used+unused) = rate of being used in a clique.
	// Thus rand < this is selecting for useful nodes.
	if (drand48() < (g->node_used[n]+adj) / (g->node_used[n]+g->node_unused[n]))
	    break;
    }
    return n;
    //return rand() % g->nnodes;
}

int random_node_low(graph *g) {
    return g->node_map2[random() % g->random_sz];
    //return rand() % g->nnodes;
}
#else
int random_neighbour(graph *g, int node) {
    return g->neighbour_count[node]
	? g->neighbour[node][rand() % g->neighbour_count[node]]
	: 0;
}

int random_node(graph *g) {
    return g->node_map[random() % g->random_sz];
    //return rand() % g->nnodes;
}

int random_node_low(graph *g) {
    return g->node_map2[random() % g->random_sz];
    //return rand() % g->nnodes;
}
#endif

void free_graph(graph *g) {
    if (g->conn) {
	for (int i = 0; i < g->nnodes; i++)
	    free(g->conn[i]);
	free(g->conn);
    }
    if (g->neighbour) {
	for (int i = 0; i < g->nnodes; i++)
	    free(g->neighbour[i]);
	free(g->neighbour);
	free(g->neighbour_count);
    }
    free(g->node_map);
    free(g->node_map2);
    free(g->degree);
    free(g->filtered);
    free(g->weight);
    free(g->node_used);
    free(g->node_unused);
    free(g->node_date);
    free(g);
}

// Remove all nodes with edge connectivity lower than "degree".
// We do this to simplify the solution.  If we've got 20 nodes and
// we've found a subset that's a clique of 12 nodes, then the
// remaining 8 must have at least 12 neighbours otherwise they
// can't be part of a larger clique.
void filter_graph(graph *g, int degree) {
    int n = g->nnodes;
    for (int i = 0; i < n; i++) {
	if (g->degree[i] < degree) {
	    fprintf(stderr, "Filter node %d degree %d\n", i, g->degree[i]);
	    g->filtered[i] = 1;
	}
    }
}

// A clique is an array of size nnodes of 0 or 1 indicating presence in a clique
// Test if it forms a clique.
// Returns 1 if clique, 0 if not.
int test_clique(graph *g, clique *q) {
    int n = g->nnodes;

    // Collapse to used nodes for faster loop below
    int cnodes[MAX_CLIQUE];
    int csize = 0;
    for (int i = 0; i < n; i++) {
	if (q->used[i])
	    cnodes[csize++] = i;
    }

    for (unsigned int i = 0; i < csize; i++) {
	int csum2 = 1;
	char *conn_row = g->conn[cnodes[i]];
	for (unsigned int j = i+1; j < csize; j++)
	    csum2 &= conn_row[cnodes[j]];

	if (!csum2) {
	    q->is_clique = 0;
	    return 0;
	}
    }

    q->is_clique = 1;
    return 1;
}

// Returns 0 if clique, or a node number if not.
// The node number is an example failing node.
int test_clique2(graph *g, clique *q) {
    int n = g->nnodes;

    // Collapse to used nodes for faster loop below
    int cnodes[MAX_CLIQUE];
    int csize = 0;
    for (int i = 0; i < n; i++) {
	if (q->used[i])
	    cnodes[csize++] = i;
    }

    // randomise cnodes order
    for (int i = 0; i < csize; i++) {
	int j = (int)(drand48() * csize);
	int t = cnodes[i];
	cnodes[i] = cnodes[j];
	cnodes[j] = t;
    }

    int worst_csum=0, worst_csum_i=-1;

    // Unrolled version, although there's barely much difference
    // with clang18 so maybe it's not worth it.
#if 0
    // 18124 gcc15
    for (unsigned int i = 0; i < csize; i++) {
	int csum2 = 1;
	char *conn_row = g->conn[cnodes[i]];
	for (unsigned int j = i+1; j < csize; j++) {
	    csum2 &= conn_row[cnodes[j]];
	}

	if (!csum2) {
	    q->is_clique = 0;
	    //printf("del %d\n", cnodes[i]);
	    return cnodes[i];
	}
    }
#else
    // 14764 gcc15
    // Unroll in i & j.
    unsigned int i;
    for (i = 0; i+4 < csize; i+=4) {
	int csum0 = 1, csum1 = 1, csum2 = 1, csum3 = 1;
	char *conn_row0 = g->conn[cnodes[i]];
	char *conn_row1 = g->conn[cnodes[i+1]];
	char *conn_row2 = g->conn[cnodes[i+2]];
	char *conn_row3 = g->conn[cnodes[i+3]];

	unsigned int j;
	for (j= i+1; j<i+4 && j < csize; j++) {
	                 csum0 &= conn_row0[cnodes[j]];
	    if (j > i+1) csum1 &= conn_row1[cnodes[j]];
	    if (j > i+2) csum2 &= conn_row2[cnodes[j]];
	    if (j > i+3) csum3 &= conn_row3[cnodes[j]];
	}
	for (j = i+4; j+4 < csize; j+=4) {
	    csum0 &= conn_row0[cnodes[j]]
		  &  conn_row0[cnodes[j+1]]
		  &  conn_row0[cnodes[j+2]]
		  &  conn_row0[cnodes[j+3]];
	    csum1 &= conn_row1[cnodes[j]]
		  &  conn_row1[cnodes[j+1]]
		  &  conn_row1[cnodes[j+2]]
		  &  conn_row1[cnodes[j+3]];
	    csum2 &= conn_row2[cnodes[j]]
		  &  conn_row2[cnodes[j+1]]
		  &  conn_row2[cnodes[j+2]]
		  &  conn_row2[cnodes[j+3]];
	    csum3 &= conn_row3[cnodes[j]]
		  &  conn_row3[cnodes[j+1]]
		  &  conn_row3[cnodes[j+2]]
		  &  conn_row3[cnodes[j+3]];
	}
	for (; j < csize; j++) {
	    csum0 &= conn_row0[cnodes[j]];
	    csum1 &= conn_row1[cnodes[j]];
	    csum2 &= conn_row2[cnodes[j]];
	    csum3 &= conn_row3[cnodes[j]];
	}

	if (!csum0) {
	    q->is_clique = 0;
	    return cnodes[i];
	}
	if (!csum1) {
	    q->is_clique = 0;
	    return cnodes[i+1];
	}
	if (!csum2) {
	    q->is_clique = 0;
	    return cnodes[i+2];
	}
	if (!csum3) {
	    q->is_clique = 0;
	    return cnodes[i+3];
	}
    }
    for (; i < csize; i++) {
	int csum2 = 1;
	char *conn_row = g->conn[cnodes[i]];
	for (unsigned int j = i+1; j < csize; j++) {
	    csum2 &= conn_row[cnodes[j]];
	}

	if (!csum2) {
	    q->is_clique = 0;
	    return cnodes[i];
	}
    }
#endif

    q->is_clique = 1;
    return 0;
}

// Clique union
void merge_clique(graph *g, clique *out, clique *in1, clique *in2) {
    //remove_chash(g, out);

    int n = g->nnodes;
    for (int i = 0; i < n; i++)
	out->used[i] = in1->used[i] | in2->used[i];
    out->size = 0;
    for (int i = 0; i < n; i++)
	out->size += out->used[i] * g->weight[i];

//    if (check_chash(g, out)) {
//	random_clique(g, out, 2);
//	improve(g, out, 9999);
//    } else {
//	add_chash(g, out);
//    }
}


// Pure intersection only
void intersect_clique2(graph *g, clique *out, clique *in1, clique *in2) {
    int n = g->nnodes;
    for (int i = 0; i < n; i++)
	out->used[i] = in1->used[i] & in2->used[i];

    out->size = 0;
    for (int i = 0; i < n; i++)
	out->size += out->used[i] * g->weight[i];

    //printf("%d x %d => %d, ", in1->size, in2->size, out->size);
}


// Clique intersection
void intersect_clique(graph *g, clique *out, clique *in1, clique *in2) {
    //remove_chash(g, out);

    int n = g->nnodes;
    for (int i = 0; i < n; i++)
	out->used[i] = in1->used[i] & in2->used[i];

    for (int i = 0; i < n; i++)
	out->used[i] |= i&1 ? in1->used[i] : in2->used[i];

    out->size = 0;
    int union_sz = 0;
    for (int i = 0; i < n; i++) {
	out->size += out->used[i] * g->weight[i];
	union_sz += (in1->used[i] | in2->used[i]) * g->weight[i];
    }

    if (out->size == 0) {
	out->used[random_node(g)] = 1;
	out->used[random_node(g)] = 1;
	out->used[random_node(g)] = 1;
	out->used[random_node(g)] = 1;
    }

    for (int i = 0; i < union_sz - out->size; i++) {
	int p1, p2;
	do {
	    do {
		p1 = random_node_low(g);
	    } while (!out->used[p1]);
	    
	    p2 = random_node_low(g);
	} while (p1 == p2 || !g->conn[p1][p2]);
	out->used[p1] = 1;
	out->used[p2] = 1;
    }

    out->size = 0;
    for (int i = 0; i < n; i++)
 	out->size += out->used[i] * g->weight[i];

//    void improve_(graph *g, clique *q, int cycles);
//    improve_(g, out, 4);

    //add_chash(g, out);
}

// Successive removal of nodes until we get a clique again
void make_clique(graph *g, clique *q) {
    //remove_chash(g, q);

    int x;
    // FIXME: in a loop this is slow as we repeatedly rebuild cnodes and
    // repeatedly randomise.  Maybe built it once, randomise it once,
    // then update as we go.
    // ANSWER: can't get it working well.  It's slower and poorer
    while ((x = test_clique2(g, q))) {
	int n;
// random, 10k cycles, 320 pop, p_hat700-1 approx cycles/100 to find clique=11
// 36 05 44 41 61 xx xx 84 59 36 21 18 15 79 84 xx 49 64 xx 64 xx 33 20 13 / 19

// test_clique2
// 08 xx 10 21 17 xx 23 15 20 20 xx 16 xx 15 16 08 46 42 20 12 09 16 21 27 / 20

// flat random, pow=1.5
// 10 10 18 33 xx xx 49 33 44 05 xx 56 51 49 33 23 36 59 05 92 26 xx 56 15 / 21

// flat random, pow=2
// xx xx xx xx xx 10 18 59 18 41 38 xx xx xx 31 79 xx xx 41 15 xx 26 36 18 / 13

// random_low, pow=1
// 08 xx 15 xx 23 77 59 xx xx xx xx 20 13 99 08 xx 41 23 44 20 61 xx 49 36 / 16

// random_low, pow=1.5
// 26 18 13 33 xx xx 26 26 28 20 18 15 xx 31 31 31 xx xx 26 08 31 99 26 xx / 19

// random_low, pow=2
// xx 20 13 31 15 15 15 10 xx 33 xx 18 41 26 69 20 xx 28 33 26 18 64 xx 46 / 19

// random_high, pow=1.5
// 20 20 10 54 46 28 05 xx 67 xx 31 54 28 20 77 54 54 15 13 15 xx 18 xx 13 / 21

#if 1
	n = x;
#elif 1
	// Or consider doing both with a "if (rand()&1)" to alternate
	// ways of choosing nodes.
	do {
	    n = rand() % g->nnodes;
	} while (!q->used[n]);
#else
	do {
	    // high: 282 282 283 283 283 283 283 284 284 285
	    // flat: 270 271 281 282 282 283 284 284 285 285
	    // low:  281 283 283 283 283 283 284 284 284 285 
	    // f/h:  277 281 281 283 283 283 283 284 284 284 alternate
	    //n = random_node_low(g);
	    n = random_node(g);
	} while (!q->used[n]);
#endif

	int connected = 1;
	for (int i = 0; i < g->nnodes; i++) {
	    if (i !=n && q->used[i] && !g->conn[i][n]) {
		connected = 0;
		break;
	    }
	}
	//printf("Remove %d\n", n);
	q->size-=g->weight[n];
	q->used[n] = 0;
    }

    q->hash = fnv1a_32(q->used, g->nnodes);
    //add_chash(g, q);
}

// TODO: Use pointers instead of copying?
void copy_clique(clique *out, clique *in) {
    memcpy(out, in, sizeof(*in));
}

static char str[1000000];
char *print_clique(graph *g, clique *q) {
    char *cp = str;
    for (int i = 0; i < g->nnodes; i++)
	if (q->used[i])
	    cp += sprintf(cp, " %d", i);
    return str;
}

void random_clique(graph *g, clique *q, int size) {
    remove_chash(g, q);

    memset(q->used, 0, g->nnodes);
    q->hash = 0;

    //int p0 = rand() % g->nnodes;
    // 44 45 45 46 46 47 48 48 48 49
    // 47 47 47 47 47 47 48 48 48 48
    // 44 46 46 46 47 48 48 48 48 48 

    // random_node doesn't help much here
    // 45 45 45 45 45 47 48 49 49 49  **1
    // 46 46 46 47 47 47 48 48 48 48  **1
    // 45 45 46 46 46 46 49 49 49 49  **1
    // 45 45 46 46 46 47 48 48 49 49  **1.5
    // 45 45 45 45 48 47 47 48 48 49  **2
    int p0 = random_node(g);

    q->used[p0] = 1;
    q->size = g->weight[p0];
    for (int i = 1; i < size; i++) {
	int p = rand() % g->nnodes, pe = (p-1 + g->nnodes) % g->nnodes;
	while ((q->used[p] || !g->conn[p0][p]) && p != pe)
	    p = (p+1) % g->nnodes;
	if (p == pe) {
	    //printf("No connection: decrement size\n");
	    size--;
	    break;
	}

	q->used[p] = 1;
	q->size += g->weight[p];
    }

    add_chash(g, q);
}

static int imp_neigh[MAX_CLIQUE];
static int imp_neigh_cnt;
static int cnodes[MAX_CLIQUE];
static int csize;

// Hill climb
int improve_(graph *g, clique *q, int cycles, double aggressive) {
    timestamp++;

    int n = g->nnodes;
    // FIXME: compute in improve() and update each cycle of improve_() instead
    if (imp_neigh_cnt < 0) {
	csize = 0;
	for (int i = 0; i < n; i++) {
	    if (q->used[i])
		cnodes[csize++] = i;
	}
    }

    int add = -1;

#define TS_DELTA 9

    if (drand48() < aggressive) {
	// Don't do this every loop
	static int last = 0;

	if (++last > 400) {
	    last = 0;
	    ac_pair pp[MAX_CLIQUE];
	    int method = rand()&1;
	    for (int i = 0; i < n; i++) {
		pp[i].i1 = i;
		//pp[i].i2 = 1000*(g->node_used[i]+1.0)/(g->node_used[i]+g->node_unused[i]);
		//pp[i].i2 = g->degree[i] + g->node_used[i];
		switch(method) {
		case 0:
		case 3:
		    pp[i].i2 = g->degree[i];
		    break;
		case 1:
		    pp[i].i2 = 1000*(g->node_used[i]+1.0)/(g->node_used[i]+g->node_unused[i]);
		    break;

		case 2: // unused with method=rand()&1
		    pp[i].i2 = g->degree[i] + g->node_used[i];
		    break;
		    //double uf = (g->node_used[i]+1.0)/(g->node_used[i]+g->node_unused[i]);
		    //pp[i].i2 = g->degree[i] + g->degree[i]*0.5 * uf;
		}
	    }
	    qsort(pp, n, sizeof(*pp), pair_cmp);
	    for (int i = 0; i < n; i++) {
		g->degree_order[i] = pp[i].i1;
		// improves non-aggressive mode initially, but harms in later
		// cycles.
		//g->degree[i] = pp[i].i2;
	    }
	}


	// OLD improver from aclique.10.c.
	int p[MAX_CLIQUE];
	int psize = 0;
	for (int i = 0; i < n; i++) {
	    if (!q->used[g->degree_order[i]]) {
		p[psize] = g->degree_order[i];
		psize++;
	    }
	}

	if (psize == 0) {
	    // one giant clique and we have nothing left
	    return 0;
	}

	for (int i = 0; i < n/10; i++) {
	    int p1 = i;
	    int p2 = rand() % psize;

	    int t = p[p1];
	    p[p1] = p[p2];
	    p[p2] = t;
	}

	// Look for a connected node.
	// TODO: We should be able to AND the entire connection row against
	// the node numbers and see what's left?
	// Or at least as a check that *ANY* node connects as a filter.
	// However good for high csize, but bad for low csize.
	for (unsigned int i = 0; i < psize; i++) {
	    int connected = 1;
	    char *conn_row = g->conn[p[i]];
	    for (unsigned int j = 0; j < csize; j++) {
		if (!conn_row[cnodes[j]]) {
		    connected = 0;
		    break;
		}
	    }

	    if (connected) {// && timestamp - g->node_date[p[i]] > TS_DELTA) {
		//cycles++;
		q->used[p[i]] = 1;
		if (!check_chash(g, q)) {
		    add = p[i];
		    break;
		} else {
		    add = -1;
		    q->used[p[i]] = 0;
		}
	    }
	}

    } else {
	// New improver from aclique.11.c.  This is more hill-climbing,
	// but sometimes hits local maxima

	// Keep a lookup table of all nodes that aren't in q->used but are neighbours
	// of everything in q->used.  Ie all candidates for increasing clique by 1.

	// DEBUG
	//    int imp_neigh_X[MAX_CLIQUE], imp_neigh_X_cnt;
	//    if (imp_neigh_cnt >= 0) {
	//	memcpy(imp_neigh_X, imp_neigh, imp_neigh_cnt*sizeof(int));
	//	imp_neigh_X_cnt = imp_neigh_cnt;
	//	imp_neigh_cnt = -1;
	//    } else {
	//	imp_neigh_X_cnt = -1;
	//    }

	if (imp_neigh_cnt < 0) {
	    memset(imp_neigh, 0, g->nnodes * sizeof(*imp_neigh));
	    imp_neigh_cnt = 0;
	    for (int i = 0; i < g->nnodes; i++) {
		if (q->used[i] || timestamp - g->node_date[i] < TS_DELTA)
		    continue;
		char *conn_row = g->conn[i];
		int connected = 1;
		for (int j = 0; j < csize; j++) {
		    if (!conn_row[cnodes[j]]) {
			connected = 0;
			break;
		    }
		}
		if (connected)
		    imp_neigh[imp_neigh_cnt++] = i;
	    }
	}

	//    if (imp_neigh_X_cnt >= 0) {
	//	if (imp_neigh_cnt != imp_neigh_X_cnt ||
	//	    memcmp(imp_neigh_X, imp_neigh, imp_neigh_cnt * sizeof(int)))
	//	    abort();
	//    }

	if (imp_neigh_cnt == 0)
	    return 0;

	// With degree-orient choice: 200 runs in 13m27s
	// 291 x 74 1.54255 37%
	// 290 x 99 1.27017 50%

	// With 10k iter 256 pop and random reset adjustments, in 13m56
	// 291 x 114 2.20731 57%
	// 290 x 84  1.57241 42%
	// A bit slower overall, but better hit ratio.

	// With 10k iter 128 pop, in 13m13s (no_gain >= 2000) <<<
	// 291 x 129 1.44519  65%
	// 291 x 70  0.837162 35%

	// With 10k iter 128 pop, in 12m5s (no_gain >= 3000) <<<
	// 291 x 129 1.376850 64%
	// 290 x  67 0.953453 34%

	// With 10k iter 128 pop, in 11m31s (no_gain >= 4000) <<<
	// 291 x 126 1.40977  63
	// 290 x  74 0.770625 37

	// With 10k iter 128 pop, in 11m11s (no_gain >= 5000)
	// 291 x 118 1.294660 59%
	// 290 x  79 0.871137 40%

	// With 12k iter 128 pop, in 13m14s
	// 291 x 125 1.42819 62%
	// 290 x 75 0.917589 38%

	// With 14k iter 64 pop, in 12m44s
	// 291 x 102 1.53306 51%
	// 290 x  94 0.80474 47%

	// With 5k iter 256 pop, in 7m47s
	// 291 x 67 1.74315  33%
	// 290 x 131 1.47643 65%

	// ---- aclique.11.c
	// a11: with 10000i 32p in 15m11s
	// 291 x 73 1.73196 36%
	// 290 x 97 1.43204 49%

	// a11: with 10000i 48p in 15m16s <<<<
	// 291 x 181 1.25394 90%
	// 290 x  16 1.01864 8%

	// a11: with 10000i 64p in 17m21s <<<
	// 291 x 176 1.33582 88%
	// 290 x  24 0.91559 12%
	//
	// a11: with 10000i 128p in 29m26s
	// 291 x 130 3.58043 65%
	// 290 x  66 2.25658 33%

	// ---- aclique.12.c.  A bit poorer on 291u, but much better on 209u
	// With 10000i 128p in 21m32s
	// 291 x  98 2.85465 49%
	// 290 x 100 1.56793 50%
	//
	// With 10000i 96p in 19m9s
	// 291 x 134 1.96163  67%
	// 290 x  66 0.922818 33%
	//
	// With 10000i 64p in 16m46s <<< (latest at 15m27s)
	// 291 x 155 1.33804  77%
	// 290 x  46 0.60152  23%
	//
	// With 10000i 48p in 14m25s
	// 291 x 140 1.22861  70%
	// 290 x  60 0.553577 30%

	// ---- aclique.14.c
	// 10000i 64p in 14m49s, but better on other data sets
	// 291 x 130 1.35654 65%
	// 290 x  68 1.10835 34%

	// Pick a random node from imp_neigh and recurse.
	int mode = rand()%2;
	int max_deg = mode == 0 ? MAX_CLIQUE : 0;
	int max_deg_cnt = 0;
	int max_deg_arr[MAX_CLIQUE];
	for (int i = 0; i < imp_neigh_cnt; i++) {
	    int x = g->degree[imp_neigh[i]] + sqrt(g->node_used[imp_neigh[i]]);
	    // NB: switched to min degree as an experiment!  BETTER
	    int keep = mode == 0 ? max_deg > x : max_deg < x;
	    if (keep) {
		max_deg = x;
		//max_deg_cnt = 0;
		// Or, keep a selection of good items
		for (int j = 0; j < max_deg_cnt; j+=2)
		    max_deg_arr[j/2] = max_deg_arr[j];
		max_deg_cnt >>= 1;
	    }
	    if (max_deg == x) {
		max_deg_arr[max_deg_cnt++] = imp_neigh[i];
	    }

	}

	// Mixture of methods: poorer
//	add = drand48() < aggressive           // * 76%/24%; 16m46
//	    ? imp_neigh[rand()%imp_neigh_cnt]  // L 74%/26%; 16m19
//	    : max_deg_arr[rand()%max_deg_cnt]; // R 75%/25%; 14m17

	// By random pick of candidates; best
	add = max_deg_arr[rand()%max_deg_cnt];

//	// By oldest first; poorest
//	int oldest = 0;
//	for (int i = 0; i < max_deg_cnt; i++) {
//	    if (oldest < timestamp - g->node_date[max_deg_arr[i]]) {
//		oldest = timestamp - g->node_date[max_deg_arr[i]];
//		add = max_deg_arr[i];
//	    }
//	}
    }

    if (add >= 0) {
	cnodes[csize++] = add;
	q->used[add] = 1;
	//printf("%d %d\n", g->node_date[add], timestamp);
	g->node_date[add] = timestamp;
	q->size+=g->weight[add];

	if (imp_neigh_cnt >= 0) {
	    // Update imp_neigh lookup table.
	    // We know everything in it matches cnodes[], but check 'add' too.
	    int j = 0;
	    for (int i = 0; i < imp_neigh_cnt; i++) {
		if (imp_neigh[i] == add)
		    continue;
		if (!g->conn[add][imp_neigh[i]])
		    continue;
		imp_neigh[j++] = imp_neigh[i];
	    }
	    imp_neigh_cnt = j;
	}


	//printf("Add node %d, size %d\n", add, q->size);
	if (cycles > 0)
	    return improve_(g, q, --cycles, aggressive);
    }
    return cycles;
}

void improve(graph *g, clique *q, int cycles, int aggressive) {
    remove_chash(g, q);
    //printf("curr=%d ", q->size);
    aggressive |= drand48()<AGGRATE;
    imp_neigh_cnt = -1; // marker for recalculate the cache.
    int c = improve_(g, q, cycles, aggressive);
    //printf("improved cycles used=%d of %d => %d\n", cycles-c, cycles, q->size);
    add_chash(g, q);
}

// Mutates by randomly removing nodes and then adding some back again via
// a few rounds of the improve function.
void mutate(graph *g, clique *q, int cnt, int aggressive) {
    remove_chash(g, q);

    int n = g->nnodes;
    int cnodes[MAX_CLIQUE];
    int csize = 0;
    for (int i = 0; i < n; i++) {
	if (q->used[i])
	    cnodes[csize++] = i;
    }

    // Remove nodes
    chash[q->hash & HMASK] = 0;
    while (cnt && q->size > 1) {
	int node;
	do {
	    node = cnodes[rand()%csize];
	} while (!q->used[node]);
	q->used[node] = 0;
	q->size-=g->weight[node];
	cnt--;
    }

    // Improve, which adds them back again
    improve(g, q, cnt, aggressive);
}

clique *find_cliques(int population, int niter, int time_limit,
		     graph *g, int *ncliques, int *top_score) {
    clock_t start_clock = clock();
    clock_t best_clock = 0;
    FILE *outfp = stderr; // for compat stdout gets swallowed up

    //int bail_iter = BO_TGROWTH;

    if (!niter) {
	//niter = exp(7+g->nnodes/150.0);
	niter = MIN(exp(4.5+g->nnodes/100.0),
		    exp(9.0+g->nnodes/500.0));
	fprintf(outfp, "Choosing %d iterations\n", niter);
    }

    int time_growth = 0; // iterations since new max best size
    int time_member = 0; // iterations since a new member of the current max
    //printf("GRAPH nnode=%d\n", g->nnodes);
    //for (int i = 0; i < g->nnodes; i++) {
    //	printf("GRAPH %d: ", g->weight[i]);
    //	for (int j=0; j < g->nnodes; j++)
    //	    putchar('0'+g->conn[i][j]);
    //	putchar('\n');
    //}

    int best_score = 0;
    static clique best[MAX_BEST];
    int best_hash = 0;
    int best_count = 0;

    if (population > MAX_POPULATION)
	population = MAX_POPULATION;

    // Seed with clique size 2
    fprintf(stderr, "Creating population\n");
    for (int i = 0; i < population; i++) {
	do {
	    random_clique(g, &pop[i], 2);
	    //printf("pop %d: %s\n", i, print_clique(g, &pop[i]));
	    improve(g, &pop[i], 9999, i < population/10);
	} while (!test_clique(g, &pop[i]));

	if (best_score < pop[i].size) {
	    best_score = pop[i].size;
	    best_hash = pop[i].hash;
	    bhash[pop[i].hash & HMASK]++;
	    best_count = 0;
	    copy_clique(&best[best_count++], &pop[i]);
	    fprintf(outfp, "New best clique of size %d time %f iter 0\n",
		    best_score,
		    (clock() - start_clock)/(double)CLOCKS_PER_SEC);
	    best_clock = clock();
	}
    }
    fprintf(stderr, "Done\n");
//    for (int j = 1; j < population; j++) {
//	printf("INIT: %d %u %d\n", j, pop[j].hash, check_chash(g, &pop[j]));
//    }
    //filter_graph(g, best_score);

    // Merge random cliques
    int no_gain = 0, skip_repop = 0, iter_end = niter;
    for (int i = 0; i < niter; i++) {
	if ((clock() - start_clock) /(double)CLOCKS_PER_SEC > time_limit)
	    break;

	time_growth++;
	time_member++;
	if ((i & 0xff) == 0)
	    fprintf(stderr, "Iter %d / %d\n", i, niter);
	int p0;
	for (int i = 0; i < 10; i++) {
	    //p0 = rand() % (population/2 + population/4);
	    p0 = rand() % population;
	    if (pop[p0].size <= best_score * 0.99)
		break;
	}

	int p1, p2;
	p1 = rand() % population;
	p2 = rand() % population;

	clique q;
	merge_clique(g, &q, &pop[p1], &pop[p2]);

	// Also try intersection of p1 and p2.
	// Followed by successive explorations of left over nodes in p1 and p2
	// Ie some characteristics of both.
	// 
	//intersect_clique(g, &q, &pop[p1], &pop[p2]);

	//intersect_clique2(g, &q, &pop[p1], &pop[p2]);

	// It's probably not a clique now, so remove things so it is.
	// 1100.1010: finds 270 in ~2min (try 1)
	//            finds 270 in ~1min (try 2)
	make_clique(g, &q);

	imp_neigh_cnt = -1; // marker for recalculate the cache.
	improve_(g, &q, 9999, (drand48()<AGGRATE?1:0) || (p1<population/2 && p2<population/2));
	//printf("%d\n", q.size);
	//improve_(g, &q, 9999, drand48()<AGGRATE?1:0);

	//fprintf(stderr, "Merged %d %d to %d\n", pop[p1].size, pop[p2].size, q.size);

	//printf("merged %d(%d) and %d(%d) into %d(%d).  Is_clique=%d %s, present %d, hash %u\n",
	//       p1, pop[p1].size,
	//       p2, pop[p2].size,
	//       p0, q.size,
	//       test_clique(g, &q),
	//       "",//print_clique(g, &q),
	//       chash[q.hash & HMASK],
	//       q.hash);

	if (best_score < q.size && q.is_clique)
	    no_gain = 0;
	else
	    no_gain++;

	//if (q.is_clique && q.size == best_score)
	//    printf("Size %d time_growth=%d\n", q.size, time_growth);

	if (q.is_clique && !check_chash(g, &q)) {
	    //improve(g, &q, 9999);
	    //printf("Improved to %d\n", q.size);
	    remove_chash(g, &pop[p0]);
	    copy_clique(&pop[p0], &q);
	    add_chash(g, &pop[p0]);

#ifdef USE_NODE_FREQ
	    if (best_score >= q.size) // a15c >= size.  a15d >= size-5
		update_node_stats(g, &q);
#endif

#if 0
	// Fails.
	}

	if (q.is_clique && !check_bhash(g, &q)) {
#endif

	    if (best_score <= q.size && best_hash != pop[p0].hash) {
		if (best_score < q.size) {
		    best_count = 0;
		    fprintf(outfp, "New best clique of size %d time %f iter %d\n",
			    q.size,
			    (clock() - start_clock)/(double)CLOCKS_PER_SEC,
			    i);
		    best_clock = clock();

		    memset(bhash, 0, (HMASK+1)*sizeof(*bhash));
		    time_growth = 0;
		    time_member = 0;
		}

		best_score = q.size;
		best_hash = pop[p0].hash;

		if (best_count < MAX_BEST && !bhash[pop[p0].hash & HMASK]) {
		    //for (int _=0;_<g->nnodes;_++)
		    //	putchar('0'+q.used[_]);
		    //putchar('\n');
		    bhash[pop[p0].hash & HMASK]++;
		    copy_clique(&best[best_count++], &pop[p0]);
		    time_member = 0;
		}
		//printf("best_count=%d time_member=%d\n", best_count, time_member);
	    }

//	    for (int j = 1; j < population; j++) {
//		printf("%d: %d %u %d\n", i, j, pop[j].hash, check_chash(g, &pop[j]));
//	    }
	}

#ifdef BAIL_OUT
	if (((time_member > BO_TMEMBER || best_count == MAX_BEST)
	     && time_growth > BO_TGROWTH)
	    || time_member > BO_TMEMBER2) {
	    fprintf(outfp, "Terminating with time_member %d and time_growth %d\n",
		    time_member, time_growth);
	    iter_end = i;
	    break;
	}
#endif

#if 1
	for (int i = 0; i < population/16; i++) {
	    int p = rand() % population;
//	    for (int j = 0; j < 10; j++) {
//		p = rand() % population;
//		if (pop[p].size <= best_score * 0.8)
//		    break;
//	    }
	    //printf("mutate %d (score %d)", p, pop[p].size);
	    mutate(g, &pop[p], 8, 1);
	    //mutate(g, &pop[p], pop[p].size/8+1, i==0);
	    //printf(" to %d (score %d)\n", p, pop[p].size);
	}
#endif

#if 1
	// TODO (try)
	// Split into two halfs.
	// 1st half: merge 1st + 2nd
	// 2nd half: random new cliques

	if (no_gain >= 3000) {
	    // Re-seed and try again
	    fprintf(stderr, "Reseeding at iter %d with score %d x %d.  "
		    "time mem/growth=%d/%d\n",
		    i, best_score, best_count,
		    time_member, time_growth);
	    memset(chash, 0, (1<<HSIZE) * sizeof(chash[0]));
	    for (int i = 0; i < population; i++) {
		do {
		    random_clique(g, &pop[i], 2);
		    improve(g, &pop[i], 9999, i < population/10);
		} while (!test_clique(g, &pop[i]));
		if (best_score < pop[i].size) {
		    best_hash = pop[i].hash;
		    bhash[pop[i].hash & HMASK]++;
		    best_count = 0;
		    copy_clique(&best[best_count++], &pop[i]);
		    fprintf(outfp, "New best clique of size %d time %f iter 0\n",
			    best_score,
			    (clock() - start_clock)/(double)CLOCKS_PER_SEC);
		    best_clock = clock();
		}
	    }
	    skip_repop = 1;
	    no_gain = 0;
	} else {
	    // Unimproved small mini cliques for purposes of merging in
	    for (int i = 0; i < population/16; i++) {
		int p;
		for (int j = 0; j < 10; j++) {
		    //p = (rand() % (population/2)) + population/2;
		    p = rand() % population;
		    if (pop[p].size <= best_score * 0.7)
			break;
		}
		//printf("random clique %d, old hash %u, score %d of %d\n", p, pop[p].hash, pop[p].size, best_score);
		random_clique(g, &pop[p], 2);
		//improve(g, &pop[p], 0);
	    }
	}
#endif

#if 0
	// Adjusted niter for 9.5s
	// if 0  23 23 23 25 25
	// P/32  22 25 25 25 25
	// P/16  22 23 23 25 25

	// For 2.1s
	// if 0  22 22 22 23 23 23 23 23 25 25
	// P/32  22 22 22 22 22 22 23 23 23 25
	for (int i = 0; i < population/32; i++) {
	    int p;
	    for (int j = 0; j < 10; j++) {
		//p = (rand() % (population/2)) + population/2;
		p = rand() % population;
		if (pop[p].size <= best_score * 0.7)
		    break;
	    }
	    random_clique(g, &pop[p], 2);
	    //improve(g, &pop[p], 9999, i==0);
	}
#endif
	// Keeping the best clique can cause us to end up in a trough.
	//printf("Update pop[0]\n");
	if (!skip_repop) {
	    for (int i = 0; i < 3; i++) {
		remove_chash(g, &pop[i]);
		copy_clique(&pop[i], &best[MAX(0,MIN(i,best_count-1))]);
		if (i > 0)
		    mutate(g, &pop[i], i*2, i==1);
		add_chash(g, &pop[i]);
	    }
//		//printf("Update pop[1]\n");
//		remove_chash(g, &pop[1]);
//		copy_clique(&pop[1], &best[MAX(0,MIN(1,best_count-1))]);
//		mutate(g, &pop[1], 2, 1);
//		add_chash(g, &pop[1]);
//
//		//printf("Update pop[2]\n");
//		remove_chash(g, &pop[2]);
//		copy_clique(&pop[2], &best[MAX(0,MIN(2,best_count-1))]);
//		mutate(g, &pop[2], 4, 0);
//		add_chash(g, &pop[2]);
	}
	skip_repop = 0;
	
	// p1 + p2 -> p0
    }

    *ncliques = best_count;
    *top_score = best_score;

    fprintf(outfp, "Finished iterations with time_member %d and time_growth %d\n",
	    time_member, time_growth);
    fprintf(outfp, "Found %4d cliques with best score %4d first_time %7.2f total time %7.2f, iter %d of %d\n",
	    best_count, best_score,
	    (best_clock - start_clock)/(double)CLOCKS_PER_SEC,
	    (clock() - start_clock)/(double)CLOCKS_PER_SEC,
	    iter_end, niter);

    return best;
}

#ifndef NO_ACLIQUE_MAIN
static void usage(FILE *fp) {
    fprintf(fp, "Usage: aclique filename.clq [-s seed] [-t timeout] [niter [population]]\n");
    exit(fp == stderr);
}

int main(int argc, char **argv) {
    //srand(0);
    srand(time(NULL) + clock());
    int v = rand();
    if (argc < 2)
	usage(stdout);

    char *fn = argv[1];
    if (strcmp(fn, "-h") == 0)
	usage(stdout);

    while (argc > 2 && argv[2][0] == '-') {
	if (argc > 3 && strcmp(argv[2], "-s") == 0) {
	    v = atoi(argv[3]);
	    argc-=2;
	    argv+=2;
	} else if (argc > 3 && strcmp(argv[2], "-t") == 0) {
	    time_limit = atoi(argv[3]);
	    argc-=2;
	    argv+=2;
	} else {
	    usage(stderr);
	}
    }
    printf("RAND %d\n", v);
    srand(v);
    srand48(rand());

    int niter = NITER;
    if (argc > 2)
	niter=atoi(argv[2]);

    int population = POPULATION;
    if (argc > 3)
	population=atoi(argv[3]);

    int fin_pop = 10;

    printf("Population %d, niter %d\n", population, niter);

    graph *g = load_graph(fn);
    if (!g)
	exit(1);

    int nclique = 0, best_score = 0;
    clique *clique = find_cliques(population, niter, time_limit,
				  /*fin_pop,*/ g, &nclique, &best_score);

    // TODO: keep an array of best cliques
    printf("Returned %d cliques of size %d\n", nclique, best_score);
    for (int b = 0; b < nclique; b++) {
	printf("Best clique %d of size %d:", b+1, best_score);
	for (int i = 0; i < g->nnodes; i++)
	    if (clique[b].used[i])
		printf(" %d", i);
	printf("\n");
    }

    free_graph(g);

    return 0;
}
#endif
