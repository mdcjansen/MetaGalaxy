#ifndef SRC_IGRAPH_H_
#define SRC_IGRAPH_H_

#include <cstdint>
#include <assert.h>
#include <cstring>
#include <iostream>

typedef double igraph_real_t;
typedef int igraph_bool_t;
typedef uint_fast32_t node_t;
typedef uint_fast64_t edge_t;

#define igraph_Calloc(n,t)    (t*) calloc( (size_t)(n), sizeof(t) )
#define igraph_Realloc(p,n,t) (t*) realloc((void*)(p), (size_t)((n)*sizeof(t)))
#define igraph_Free(p)        (free( (void *)(p) ), (p) = NULL)
#define VECTOR(v) ((v).stor_begin)
#define IGRAPH_FROM(g,e) ((node_t)(VECTOR((g)->from)[(e)]))
#define IGRAPH_TO(g,e)   ((node_t)(VECTOR((g)->to)  [(e)]))
#define IGRAPH_OTHER(g,e,v) ((node_t)(IGRAPH_TO(g,(e))==(v) ? IGRAPH_FROM((g),(e)) : IGRAPH_TO((g),(e))))

#define N 624   /* Period parameters */
#define M 397

typedef struct {
	unsigned long mt[N];
	int mti;
} igraph_i_rng_mt19937_state_t;

int igraph_rng_mt19937_seed(void *vstate, unsigned long int seed) {
	igraph_i_rng_mt19937_state_t *state = (igraph_i_rng_mt19937_state_t *) vstate;
	int i;

	memset(state, 0, sizeof(igraph_i_rng_mt19937_state_t));

	if (seed == 0) {
		seed = 4357; /* the default seed is 4357 */
	}
	state->mt[0] = seed & 0xffffffffUL;

	for (i = 1; i < N; i++) {
		/* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
		 Ed. p.106 for multiplier. */
		state->mt[i] = (1812433253UL * (state->mt[i - 1] ^ (state->mt[i - 1] >> 30)) + (unsigned long) i);
		state->mt[i] &= 0xffffffffUL;
	}

	state->mti = i;
	return 0;
}

int igraph_rng_mt19937_init(void **state) {
	igraph_i_rng_mt19937_state_t *st;

	st = igraph_Calloc(1, igraph_i_rng_mt19937_state_t);
	if (!st) {
		exit(1);
	}
	(*state) = st;

	igraph_rng_mt19937_seed(st, 0);

	return 0;
}

void igraph_rng_mt19937_destroy(void *vstate) {
	igraph_i_rng_mt19937_state_t *state = (igraph_i_rng_mt19937_state_t*) vstate;
	igraph_Free(state);
}

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;
/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;

unsigned long int igraph_rng_mt19937_get(void *vstate) {
	igraph_i_rng_mt19937_state_t *state = (igraph_i_rng_mt19937_state_t*) vstate;

	unsigned long k;
	unsigned long int * const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfUL : 0)

	if (state->mti >= N) {
		/* generate N words at one time */
		int kk;

		for (kk = 0; kk < N - M; kk++) {
			unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
		}
		for (; kk < N - 1; kk++) {
			unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
		}

		{
			unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
			mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
		}

		state->mti = 0;
	}

#undef MAGIC

	/* Tempering */

	k = mt[state->mti];
	k ^= (k >> 11);
	k ^= (k << 7) & 0x9d2c5680UL;
	k ^= (k << 15) & 0xefc60000UL;
	k ^= (k >> 18);

	state->mti++;

	return k;
}

igraph_real_t igraph_rng_mt19937_get_real(void *vstate) {
	return igraph_rng_mt19937_get(vstate) / 4294967296.0;
}

typedef struct igraph_rng_type_t {
	const char *name;
	unsigned long int min;
	unsigned long int max;
	int (*init)(void **state);
	void (*destroy)(void *state);
	int (*seed)(void *state, unsigned long int seed);
	unsigned long int (*get)(void *state);
	igraph_real_t (*get_real)(void *state);
	igraph_real_t (*get_norm)(void *state);
	igraph_real_t (*get_geom)(void *state, igraph_real_t p);
	igraph_real_t (*get_binom)(void *state, long int n, igraph_real_t p);
	igraph_real_t (*get_exp)(void *state, igraph_real_t rate);
} igraph_rng_type_t;

typedef struct igraph_rng_t {
	const igraph_rng_type_t *type;
	void *state;
	int def;
} igraph_rng_t;

const igraph_rng_type_t igraph_rngtype_mt19937 = {
/* name= */"MT19937",
/* min=  */0,
/* max=  */0xffffffffUL,
/* init= */igraph_rng_mt19937_init,
/* destroy= */igraph_rng_mt19937_destroy,
/* seed= */igraph_rng_mt19937_seed,
/* get= */igraph_rng_mt19937_get,
/* get_real= */igraph_rng_mt19937_get_real,
/* get_norm= */0,
/* get_geom= */0,
/* get_binom= */0,
/* get_exp= */0 };

igraph_i_rng_mt19937_state_t igraph_i_rng_default_state;
#define addr(a) (&a)
igraph_rng_t igraph_i_rng_default = { addr(igraph_rngtype_mt19937), addr(igraph_i_rng_default_state),
/* def= */1 };

igraph_rng_t *igraph_rng_default() {
	return &igraph_i_rng_default;
}

int igraph_rng_seed(igraph_rng_t *rng, unsigned long int seed) {
	const igraph_rng_type_t *type = rng->type;
	rng->def = 0;
	type->seed(rng->state, seed);
	return 0;
}

long int igraph_rng_get_integer(igraph_rng_t *rng, long int l, long int h) {
	const igraph_rng_type_t *type = rng->type;
	if (type->get_real) {
		return (long int) (type->get_real(rng->state) * (h - l + 1) + l);
	} else if (type->get) {
		unsigned long int max = type->max;
		return (long int) (type->get(rng->state) / ((double) max + 1) * (h - l + 1) + l);
	}
	//IGRAPH_ERROR("Internal random generator error", IGRAPH_EINTERNAL);
	exit(1);
	return 0;
}

#define RNG_BEGIN()      if (igraph_rng_default()->def==1) {	\
  igraph_rng_seed(igraph_rng_default(), time(0));		\
  igraph_rng_default()->def=2;					\
  }
#define RNG_END()		/* do nothing */

#define RNG_INTEGER(l,h) (igraph_rng_get_integer(igraph_rng_default(),(l),(h)))

#undef N
#undef M

template <typename T>
struct igraph_vector_t {
	T* stor_begin; //storage begin
	T* stor_end;   //storage end
	T* end;        //current pointer
};

typedef igraph_vector_t<uint_least32_t> igraph_node_vector_t;
typedef igraph_vector_t<uint_least64_t> igraph_edge_vector_t;
typedef igraph_vector_t<float> igraph_weight_vector_t;

typedef struct igraph_s {
	node_t n;
	igraph_bool_t directed;
	igraph_node_vector_t from;
	igraph_node_vector_t to;
	igraph_edge_vector_t *incs; //incidence list which has edge id instead of node id (compared to adjacent list)
} igraph_t;

template <typename T>
int igraph_vector_init(igraph_vector_t<T>* v, uint_fast64_t size) {
	uint_fast64_t alloc_size = size > 0 ? size : 1;
	if (size < 0) {
		size = 0;
	}
	v->stor_begin = igraph_Calloc(alloc_size, T);
	if (v->stor_begin == 0) {
		std::cerr << "cannot init vector" << std::endl;
		exit(1);
	}
	v->stor_end = v->stor_begin + alloc_size;
	v->end = v->stor_begin + size;

	return 0;
}

template <typename T>
void igraph_vector_destroy(igraph_vector_t<T>* v) {
	assert(v != 0);
	if (v->stor_begin != 0) {
		igraph_Free(v->stor_begin);
		v->stor_begin = NULL;
	}
}

int igraph_destroy(igraph_t *graph) {
	igraph_vector_destroy(&graph->from);
	igraph_vector_destroy(&graph->to);
	return 0;
}

template <typename T>
uint_fast64_t igraph_vector_size(const igraph_vector_t<T>* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	return v->end - v->stor_begin;
}

template <typename T>
T igraph_vector_min(const igraph_vector_t<T>* v) {
	T min;
	T *ptr;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	min = *(v->stor_begin);
	ptr = v->stor_begin + 1;
	while (ptr < v->end) {
		if ((*ptr) < min) {
			min = *ptr;
		}
		ptr++;
	}
	return min;
}

node_t igraph_vcount(const igraph_t *graph) {
	return graph->n;
}

edge_t igraph_ecount(const igraph_t *graph) {
	return (edge_t) igraph_vector_size(&graph->from);
}

template <typename T>
void igraph_vector_null(igraph_vector_t<T>* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	if (igraph_vector_size(v) > 0) {
		memset(v->stor_begin, 0, sizeof(T) * igraph_vector_size(v));
	}
}

template <typename T>
int igraph_vector_reserve(igraph_vector_t<T>* v, uint_fast64_t size) {
	long int actual_size = igraph_vector_size(v);
	T *tmp;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	if (size <= igraph_vector_size(v)) {
		return 0;
	}

	tmp = igraph_Realloc(v->stor_begin, size, T);
	if (tmp == 0) {
		std::cerr << "cannot reserve space for vector" << std::endl;
		exit(1);
	}
	v->stor_begin = tmp;
	v->stor_end = v->stor_begin + size;
	v->end = v->stor_begin + actual_size;

	return 0;
}

template <typename T>
int igraph_vector_push_back(igraph_vector_t<T>* v, T e) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);

	/* full, allocate more storage */
	if (v->stor_end == v->end) {
		uint_fast64_t new_size = igraph_vector_size(v) * 2;
		if (new_size == 0) {
			new_size = 1;
		}
		igraph_vector_reserve(v, new_size);
	}

	*(v->end) = e;
	v->end += 1;

	return 0;
}

template <typename T>
int igraph_vector_resize(igraph_vector_t<T>* v, uint_fast64_t newsize) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	igraph_vector_reserve(v, newsize);
	v->end = v->stor_begin + newsize;
	return 0;
}

template <typename T>
int igraph_vector_resize_min(igraph_vector_t<T>*v) {
	uint_fast64_t size;
	T *tmp;
	if (v->stor_end == v->end) {
		return 0;
	}

	size = (uint_fast64_t) (v->end - v->stor_begin);
	if (size == 0)
		return 0;

	tmp = igraph_Realloc(v->stor_begin, size, T);
	if (tmp == 0) {
		std::cerr << "cannot resize vector" << std::endl;
		exit(1);
	} else {
		v->stor_begin = tmp;
		v->stor_end = v->end = v->stor_begin + size;
	}

	return 0;
}

template <typename T>
int igraph_vector_shuffle(igraph_vector_t<T> *v) {
	long int n = igraph_vector_size(v);
	long int k;
	T dummy;

	RNG_BEGIN();
	while (n > 1) {
		k = RNG_INTEGER(0, n - 1);
		n--;
		dummy = VECTOR(*v)[n];
		VECTOR(*v)[n] = VECTOR(*v)[k];
		VECTOR(*v)[k] = dummy;
	}RNG_END();

	return 0;
}

template <typename T>
void igraph_vector_clear(igraph_vector_t<T>* v) {
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	v->end = v->stor_begin;
}

template <typename T>
void igraph_vector_fill(igraph_vector_t<T>* v, T e) {
	T *ptr;
	assert(v != NULL);
	assert(v->stor_begin != NULL);
	for (ptr = v->stor_begin; ptr < v->end; ptr++) {
		*ptr = e;
	}
}

int igraph_add_vertices(igraph_t *graph, node_t nv, void *attr) {

	graph->n += nv;

	return 0;
}

int igraph_empty_attrs(igraph_t *graph, node_t n, igraph_bool_t directed, void* attr) {

	graph->n = 0;
	graph->directed = directed;
	igraph_vector_init(&graph->from, 0);
	igraph_vector_init(&graph->to, 0);

	igraph_add_vertices(graph, n, 0);

	return 0;
}

int igraph_empty(igraph_t *graph, node_t n, igraph_bool_t directed) {
	return igraph_empty_attrs(graph, n, directed, 0);
}

void igraph_inclist_destroy(igraph_t* g) {
	for (node_t i = 0; i < g->n; i++) {
		/* This works if some igraph_vector_t's are 0, because igraph_vector_destroy can
		 handle this. */
		igraph_vector_destroy(&g->incs[i]);
	}
	igraph_Free(g->incs);
}

#endif /* SRC_IGRAPH_H_ */
