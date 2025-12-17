// TODO: switch to a bit-vector for clique
// TODO: Fix repeated identical hashes for e.g. gen400_p0.9_75.clq

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>

#define MAX_CLIQUE 4000

//#define POPULATION 256
//#define NITER 10000

#define MAX_POPULATION 1024

#ifndef POPULATION
#    define POPULATION 32
#endif

#ifndef NITER
#    define NITER 10000
#endif

#ifndef AGGRATE
#    define AGGRATE 0.01
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
} pair;

int pair_cmp(const void *vp1, const void *vp2) {
    const pair *p1 = (const pair *)vp1;
    const pair *p2 = (const pair *)vp2;
    return p2->i2 - p1->i2;
}

typedef struct {
    char **conn;
    int *weight;
    int nnodes, nfiltered;
    int max_possible;
    int *filtered;
    int *degree;
    int random_sz;
    int *node_map; // for picking nodes biased to degree
    int degree_order[MAX_CLIQUE];
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

// Find the number of connections to a node
int degree(graph *g, int n) {
    int d = 0;
    for (int i = 0; i < g->nnodes; i++) {
	if (i != n)
	    d += g->conn[n][i];
    }
    return d;
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
    }

    g->nnodes = n;
    g->max_possible = n;
    g->degree = (int *)malloc(n * sizeof(*g->degree));
    g->filtered = (int *)calloc(n, sizeof(*g->filtered));
    g->conn = (char **)calloc(n, sizeof(*g->conn));
    g->weight = (int *)malloc(n * sizeof(*g->weight));
    for (int i = 0; i < n; i++) {
	g->conn[i] = (char *)calloc(n, 1);
	g->weight[i] = 1;
	//g->conn[i][i] = 1; // why does this not work?
    }
    
    int n1, n2, d, w;
    while ((d = fscanf(fp, "e %d %d %d", &n1, &n2, &w)) >= 2) {
	if (d == 2)
	    w = 1;

	int c;
	if (d == 3) {
	    do {
		c = fgetc(fp);
	    } while (c != '\n' && c != EOF);
	}

	n1--;
	n2--;
	if (n1 < 0 || n1 >= n ||
	    n2 < 0 || n2 >= n) {
	    fprintf(stderr, "Node number out of range\n");
	    goto err;
	}

	g->conn[n1][n2] = 1;
	g->conn[n2][n1] = 1;
	g->weight[n1] = w;
    }

    for (int i = 0; i < n; i++) {
	g->degree[i] = degree(g, i);
    }

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
    pair p[MAX_CLIQUE];
    for (int i = 0; i < n; i++) {
	p[i].i1 = i;
	p[i].i2 = g->degree[i];
	total_degree += p[i].i2;
	//total_degree += pow(p[i].i2, 1.5);
    }

    g->random_sz = roundup(total_degree);
    g->node_map = (int *)malloc(g->random_sz * sizeof(g->node_map));

    qsort(p, n, sizeof(*p), pair_cmp);
    for (int i = 0; i < n; i++)
	g->degree_order[i] = p[i].i1;

    // map random value to node number such that high degrees are first.
    int k, i;
    for (k = i = 0; i < g->random_sz && k < n; k++) {
	for (int j = 0;
	     j < p[k].i2 * (g->random_sz/(double)total_degree)
		 && i < g->random_sz;
	     j++, i++) {
	    g->node_map[i] = p[k].i1;
	}
    }
    while (i < g->random_sz)
	g->node_map[i++] = p[0].i1;
//    for (int i = 0; i < g->random_sz; i++) {
//	printf("%d\t%d\t%d\n", i, g->node_map[i], g->degree[g->node_map[i]]);
//    }
    
    fclose(fp);
    return g;

 err:
    free(g);
    fclose(fp);
    return NULL;
}

int random_node(graph *g) {
    return g->node_map[random() % g->random_sz];
    //return rand() % g->nnodes;
}

void free_graph(graph *g) {
    if (g->conn) {
	for (int i = 0; i < g->nnodes; i++)
	    free(g->conn[i]);
	free(g->conn);
    }
    free(g->node_map);
    free(g->degree);
    free(g->filtered);
    free(g->weight);
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
		p1 = random_node(g);
	    } while (!out->used[p1]);
	    
	    p2 = random_node(g);
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

    while (!test_clique(g, q)) {
	int n;
	do {
	    // We want inverse of random_node here, so it's biased to
	    // low degree.  We don't have it so use a flat distrib instead.
	    n = rand() % g->nnodes;
	} while (!q->used[n]);

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

// Hill climb
int improve_(graph *g, clique *q, int cycles, double aggressive) {
    int n = g->nnodes;
    int cnodes[MAX_CLIQUE];
    int csize = 0;
    for (int i = 0; i < n; i++) {
	if (q->used[i])
	    cnodes[csize++] = i;
    }

    int p[MAX_CLIQUE];
    int psize = 0;
    for (int i = 0; i < n; i++) {
	if (!q->used[g->degree_order[i]]) {
	    p[psize] = g->degree_order[i];
	    psize++;
	}
    }

    // 46 46 46 47 47 47 48 48 48 48  base (p[psize] = i, not degree_order)
    // 48 48 49 49 49 50 50 50 50 50  n/2
    // 48 49 50 50 50 50 50 50 51 51  n/10
    // 49 49 49 49 50 50 50 50 50 50  n/20

    if (drand48() > aggressive) {
	for (int i = 0; i < n/10; i++) {
	    int p1 = i;
	    int p2 = rand() % psize;

	    int t = p[p1];
	    p[p1] = p[p2];
	    p[p2] = t;
	}
    }

    // Look for a connected node.
    // TODO: We should be able to AND the entire connection row against
    // the node numbers and see what's left?
    // Or at least as a check that *ANY* node connects as a filter.
    // However good for high csize, but bad for low csize.
    int add = -1;
    for (unsigned int i = 0; i < psize; i++) {
	int connected = 1;
	char *conn_row = g->conn[p[i]];
	for (unsigned int j = 0; j < csize; j++) {
	    if (!conn_row[cnodes[j]]) {
		connected = 0;
		break;
	    }
	}

	if (connected) {
	    q->used[p[i]] = 1;
	    if (q->hash == 666) printf("hash=%p %d\n", q, q->hash);
	    if (!check_chash(g, q)) {
		add = p[i];
		break;
	    } else {
		add = -1;
		q->used[p[i]] = 0;
	    }
	}
    }

    if (add >= 0) {
	q->used[add] = 1;
	q->size+=g->weight[add];

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
    int c = improve_(g, q, cycles, aggressive);
    //printf("improved cycles used=%d of %d => %d\n", cycles-c, cycles, q->size);
    add_chash(g, q);
}

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

#define MAX_BEST 1000
clique *find_cliques(int population, int niter,
		     graph *g, int *ncliques, int *top_score) {
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

	if (best_score <= pop[i].size) {
	    best_score = pop[i].size;
	    best_hash = pop[i].hash;
	    bhash[pop[i].hash & HMASK]++;
	    best_count = 0;
	    copy_clique(&best[best_count++], &pop[i]);
	    fprintf(stderr, "New best clique of size %d\n", best_score);
	}
    }
    fprintf(stderr, "Done\n");
//    for (int j = 1; j < population; j++) {
//	printf("INIT: %d %u %d\n", j, pop[j].hash, check_chash(g, &pop[j]));
//    }
    //filter_graph(g, best_score);

    // Merge random cliques
    int no_gain = 0;
    for (int i = 0; i < niter; i++) {
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

	// It's probably not a clique now, so remove things so it is.
	// 1100.1010: finds 270 in ~2min (try 1)
	//            finds 270 in ~1min (try 2)
	make_clique(g, &q);
	improve_(g, &q, 9999, drand48()<AGGRATE?1:0);

	//printf("merged %d(%d) and %d(%d) into %d(%d).  Is_clique=%d %s, present %d, hash %u\n",
	//       p1, pop[p1].size,
	//       p2, pop[p2].size,
	//       p0, q.size,
	//       test_clique(g, &q),
	//       "",//print_clique(g, &q),
	//       chash[q.hash & HMASK],
	//       q.hash);

	if (best_score < q.size)
	    no_gain = 0;
	else
	    no_gain++;

	if (q.is_clique && !check_chash(g, &q)) {
	    //improve(g, &q, 9999);
	    //printf("Improved to %d\n", q.size);
	    remove_chash(g, &pop[p0]);
	    copy_clique(&pop[p0], &q);
	    add_chash(g, &pop[p0]);

	    if (best_score <= q.size && best_hash != pop[p0].hash) {
		if (best_score < q.size) {
		    best_count = 0;
		    printf("New best clique of size %d, hash %u, idx %d\n",
			   q.size, pop[p0].hash, p0);
		    memset(bhash, 0, (HMASK+1)*sizeof(*bhash));
		}

		best_score = q.size;
		best_hash = pop[p0].hash;

		if (best_count < MAX_BEST && !bhash[pop[p0].hash & HMASK]) {
		    //for (int _=0;_<g->nnodes;_++)
		    //	putchar('0'+q.used[_]);
		    //putchar('\n');
		    bhash[pop[p0].hash & HMASK]++;
		    copy_clique(&best[best_count++], &pop[p0]);
		}
	    }

//	    for (int j = 1; j < population; j++) {
//		printf("%d: %d %u %d\n", i, j, pop[j].hash, check_chash(g, &pop[j]));
//	    }
	}

#if 1
	// For 2.1s
	// if 0  22 22 22 23 23 23 23 23 25 25
	// P/16  22 22 22 22 22 23 23 24 25 25  cnt=2
	// P/16  21 22 22 23 23 23 23 23 23 25  cnt=size/8+1
	// P/16  22 22 22 22 23 23 23 25 25 25  cnt=size/4
	// P/16  21 22 22 22 22 22 22 22 23 25  cnt=size/4+1
	// P/16  22 22 22 22 22 22 23 23 23 23  cnt=size/2
	// P/4   22 22 22 22 23 23 23 25 25 25  cnt=size/8+1
	// P/2   21 22 22 22 22 24 24 25 25 25  cnt=size/8+1
	//
	// P/4   23 23 23 23 23 24 24 24 24 24 (using degree order)
	for (int i = 0; i < population/4; i++) {
	    int p;
	    for (int j = 0; j < 10; j++) {
		p = rand() % population;
		if (pop[p].size <= best_score * 0.8)
		    break;
	    }
	    //printf("mutate %d (score %d)", p, pop[p].size);
	    mutate(g, &pop[p], pop[p].size/8+1, i==0);
	    //printf(" to %d (score %d)\n", p, pop[p].size);
	}
#endif
	
	// TODO (try)
	// Split into two halfs.
	// 1st half: merge 1st + 2nd
	// 2nd half: random new cliques

	if (no_gain >= 1000) {
	    // Re-seed and try again
	    fprintf(stderr, "Reseeding\n");
	    for (int i = 0; i < population; i++) {
		int p = rand() % population;
		random_clique(g, &pop[p], 2);
		improve(g, &pop[p], 9999, 0);
	    }
	    no_gain = 0;
	} else {
	    // Unimproved small mini cliques for purposes of merging in
	    for (int i = 0; i < population/2; i++) {
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
	    improve(g, &pop[p], 9999, i==0);
	}
#endif
	// Keeping the best clique can cause us to end up in a trough.
	//printf("Update pop[0]\n");
	remove_chash(g, &pop[0]);
	copy_clique(&pop[0], &best[0]);
	add_chash(g, &pop[0]);

	//printf("Update pop[1]\n");
	remove_chash(g, &pop[1]);
	copy_clique(&pop[1], &best[MAX(0,MIN(1,best_count-1))]);
	mutate(g, &pop[1], 2, 1);
	add_chash(g, &pop[1]);

	//printf("Update pop[2]\n");
	remove_chash(g, &pop[2]);
	copy_clique(&pop[2], &best[MAX(0,MIN(2,best_count-1))]);
	mutate(g, &pop[2], 4, 0);
	add_chash(g, &pop[2]);
	
	// p1 + p2 -> p0
    }

    *ncliques = best_count;
    *top_score = best_score;

    return best;
}

int main(int argc, char **argv) {
    //srand(0);
    srand(time(NULL) + clock());
    int v = rand();
    printf("RAND=%d\n", v);
    // RAND=1425238825: many hits on DIMACS_subset_ascii/brock200_2.clq
    // RAND=1121685628: impossible clique_size
    srand(v);
    srand48(rand());

    if (argc == 1) {
	fprintf(stderr, "Usage: aclique filename.clq\n");
	exit(1);
    }

    graph *g = load_graph(argv[1]);
    int nclique = 0, best_score = 0;
    clique *clique = find_cliques(POPULATION, NITER, g, &nclique, &best_score);

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
