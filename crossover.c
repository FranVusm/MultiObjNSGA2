/* Ordered crossover for perm/cuts/alpha chromosomes */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

static void ordered_crossover_perm (int *p1, int *p2, int *child, int start, int end);
static int contains_customer (int *values, int start, int end, int customer);
static double alpha_upper_bound (int slot);
static double clip_alpha (double value, double upper);

void crossover (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    int i;
    int start;
    int end;
    int temp;

    if (randomperc() <= pcross_bin && n_clients >= 2)
    {
        nbincross++;
        start = rnd(0, n_clients-1);
        end = rnd(0, n_clients-1);
        if (start > end)
        {
            temp = start;
            start = end;
            end = temp;
        }
        ordered_crossover_perm(parent1->perm, parent2->perm, child1->perm, start, end);
        ordered_crossover_perm(parent2->perm, parent1->perm, child2->perm, start, end);
    }
    else
    {
        for (i=0; i<n_clients; i++)
        {
            child1->perm[i] = parent1->perm[i];
            child2->perm[i] = parent2->perm[i];
        }
    }

    for (i=0; i<n_slots+1; i++)
    {
        child1->cuts[i] = (randomperc() < 0.5) ? parent1->cuts[i] : parent2->cuts[i];
        child2->cuts[i] = (randomperc() < 0.5) ? parent2->cuts[i] : parent1->cuts[i];
    }
    for (i=0; i<n_slots; i++)
    {
        double average;
        double upper;
        upper = alpha_upper_bound(i);
        average = (parent1->alpha[i] + parent2->alpha[i]) / 2.0;
        child1->alpha[i] = clip_alpha(average + rndreal(-0.5, 0.5), upper);
        child2->alpha[i] = clip_alpha(average + rndreal(-0.5, 0.5), upper);
    }
    repair_cuts_capacity_aware(child1);
    repair_cuts_capacity_aware(child2);
    return;
}

void realcross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    return;
}

void bincross (individual *parent1, individual *parent2, individual *child1, individual *child2)
{
    crossover(parent1, parent2, child1, child2);
    return;
}

static void ordered_crossover_perm (int *p1, int *p2, int *child, int start, int end)
{
    int i;
    int fill;
    for (i=0; i<n_clients; i++)
    {
        child[i] = -1;
    }
    for (i=start; i<=end; i++)
    {
        child[i] = p1[i];
    }
    fill = 0;
    for (i=0; i<n_clients; i++)
    {
        if (contains_customer(child, start, end, p2[i]))
        {
            continue;
        }
        while (fill < n_clients && child[fill] != -1)
        {
            fill++;
        }
        if (fill < n_clients)
        {
            child[fill] = p2[i];
        }
    }
}

static int contains_customer (int *values, int start, int end, int customer)
{
    int i;
    for (i=start; i<=end; i++)
    {
        if (values[i] == customer)
        {
            return 1;
        }
    }
    return 0;
}

static double alpha_upper_bound (int slot)
{
    double horizon;
    int i;
    double min_li;
    double max_ls;
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
