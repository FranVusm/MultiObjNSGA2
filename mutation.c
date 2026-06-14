/* Mutation routines for perm/cuts/alpha chromosomes */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

static void swap_mutation (individual *ind);
static void insert_mutation (individual *ind);
static void cut_shift_mutation (individual *ind);
static void cut_reset_mutation (individual *ind);
static void alpha_mutation (individual *ind);
static double alpha_upper_bound (int slot);
static double clip_alpha (double value, double upper);

void mutation_pop (population *pop)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        mutation_ind(&(pop->ind[i]));
    }
    return;
}

void mutation_ind (individual *ind)
{
    if (randomperc() <= pmut_bin)
    {
        swap_mutation(ind);
    }
    if (randomperc() <= pmut_bin)
    {
        insert_mutation(ind);
    }
    if (randomperc() <= pmut_bin)
    {
        cut_shift_mutation(ind);
    }
    if (randomperc() <= 0.35)
    {
        cut_reset_mutation(ind);
    }
    if (randomperc() <= pmut_bin)
    {
        alpha_mutation(ind);
    }
    repair_cuts_capacity_aware(ind);
    return;
}

void bin_mutate_ind (individual *ind)
{
    mutation_ind(ind);
    return;
}

void real_mutate_ind (individual *ind)
{
    return;
}

static void swap_mutation (individual *ind)
{
    int i;
    int j;
    int temp;
    if (n_clients < 2) return;
    i = rnd(0, n_clients-1);
    j = rnd(0, n_clients-1);
    temp = ind->perm[i];
    ind->perm[i] = ind->perm[j];
    ind->perm[j] = temp;
    nbinmut++;
}

static void insert_mutation (individual *ind)
{
    int source;
    int destination;
    int customer;
    int i;
    if (n_clients < 2) return;
    source = rnd(0, n_clients-1);
    destination = rnd(0, n_clients-1);
    if (source == destination) return;
    customer = ind->perm[source];
    if (source < destination)
    {
        for (i=source; i<destination; i++)
        {
            ind->perm[i] = ind->perm[i+1];
        }
    }
    else
    {
        for (i=source; i>destination; i--)
        {
            ind->perm[i] = ind->perm[i-1];
        }
    }
    ind->perm[destination] = customer;
    nbinmut++;
}

static void cut_shift_mutation (individual *ind)
{
    int index;
    int delta;
    if (n_slots <= 1) return;
    index = rnd(1, n_slots-1);
    delta = (randomperc() < 0.5) ? -1 : 1;
    ind->cuts[index] += delta;
    nbinmut++;
}

static void cut_reset_mutation (individual *ind)
{
    int i;
    ind->cuts[0] = 0;
    for (i=1; i<n_slots; i++)
    {
        ind->cuts[i] = rnd(0, n_clients);
    }
    ind->cuts[n_slots] = n_clients;
    nbinmut++;
}

static void alpha_mutation (individual *ind)
{
    int slot;
    double upper;
    if (n_slots < 1) return;
    slot = rnd(0, n_slots-1);
    upper = alpha_upper_bound(slot);
    if (randomperc() < 0.1)
    {
        ind->alpha[slot] = rndreal(0.0, upper);
    }
    else
    {
        ind->alpha[slot] = clip_alpha(ind->alpha[slot] + rndreal(-1.0, 1.0), upper);
    }
    nbinmut++;
}

static double alpha_upper_bound (int slot)
{
    double min_li;
    double max_ls;
    double horizon;
    int i;
    min_li = problem.LI[0];
    max_ls = problem.LS[0];
    for (i=1; i<problem.nP; i++)
    {
        if (problem.LI[i] < min_li) min_li = problem.LI[i];
        if (problem.LS[i] > max_ls) max_ls = problem.LS[i];
    }
    horizon = max_ls - min_li;
    if (problem.slots[slot].trip == problem.V[0])
    {
        return minimum(1.0, horizon * 0.1);
    }
    return 0.25;
}

static double clip_alpha (double value, double upper)
{
    if (value < 0.0) return 0.0;
    if (value > upper) return upper;
    return value;
}
