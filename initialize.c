/* Initialization and structural repair for perm/cuts/alpha chromosomes */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

static double planning_horizon_c (void);
static double max_ls (void);
static double min_li (void);
static void random_capacity_feasible_cuts (individual *ind);
static int cuts_capacity_valid (individual *ind);
static void normalize_cuts (individual *ind);

void initialize_pop (population *pop)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        initialize_ind (&(pop->ind[i]));
    }
    return;
}

void initialize_ind (individual *ind)
{
    int i;
    int j;
    double horizon;
    double first_upper;
    double upper;
    int temp;

    for (i=0; i<n_clients; i++)
    {
        ind->perm[i] = problem.customers[i];
    }
    for (i=0; i<n_clients; i++)
    {
        j = rnd(i, n_clients-1);
        temp = ind->perm[i];
        ind->perm[i] = ind->perm[j];
        ind->perm[j] = temp;
    }

    random_capacity_feasible_cuts(ind);

    horizon = planning_horizon_c();
    first_upper = minimum(1.0, horizon * 0.1);
    for (i=0; i<n_slots; i++)
    {
        upper = (problem.slots[i].trip == problem.V[0]) ? first_upper : 0.25;
        if (upper <= 0.0)
        {
            ind->alpha[i] = 0.0;
        }
        else
        {
            ind->alpha[i] = rndreal(0.0, upper);
        }
    }
    return;
}

void repair_cuts_capacity_aware (individual *ind)
{
    normalize_cuts(ind);
    if (cuts_capacity_valid(ind))
    {
        return;
    }
    random_capacity_feasible_cuts(ind);
}

int validate_individual_structure (individual *ind)
{
    int i;
    int j;
    int found;
    int first_empty;

    if (ind->cuts[0] != 0 || ind->cuts[n_slots] != n_clients)
    {
        return 0;
    }
    for (i=0; i<n_slots; i++)
    {
        if (ind->cuts[i] > ind->cuts[i+1] || ind->cuts[i] < 0 || ind->cuts[i] > n_clients)
        {
            return 0;
        }
    }
    for (i=0; i<n_clients; i++)
    {
        found = 0;
        for (j=0; j<n_clients; j++)
        {
            if (ind->perm[j] == problem.customers[i])
            {
                found++;
            }
        }
        if (found != 1)
        {
            return 0;
        }
    }
    for (j=0; j<problem.nK; j++)
    {
        first_empty = 0;
        for (i=0; i<n_slots; i++)
        {
            if (problem.slots[i].vehicle != problem.K[j])
            {
                continue;
            }
            if (ind->cuts[i] == ind->cuts[i+1])
            {
                first_empty = 1;
            }
            else if (first_empty)
            {
                return 0;
            }
        }
    }
    return 1;
}

static void random_capacity_feasible_cuts (individual *ind)
{
    int groups[1024];
    int group_count;
    int current_size;
    double current_demand;
    int i;
    int idx;
    int active_slots;
    int splittable_count;
    int target_group;
    int left;
    int pos;

    group_count = 0;
    current_size = 0;
    current_demand = 0.0;
    for (i=0; i<n_clients; i++)
    {
        double demand;
        demand = problem.d[node_index(&problem, ind->perm[i])];
        if (current_size > 0 && current_demand + demand > problem.q)
        {
            groups[group_count++] = current_size;
            current_size = 0;
            current_demand = 0.0;
        }
        current_size++;
        current_demand += demand;
    }
    if (current_size > 0)
    {
        groups[group_count++] = current_size;
    }
    if (group_count == 0)
    {
        group_count = 1;
        groups[0] = 0;
    }

    active_slots = group_count;
    if (n_clients > active_slots && n_slots > active_slots && randomperc() > 0.25)
    {
        active_slots = rnd(group_count, minimum(n_slots, n_clients));
    }
    while (group_count < active_slots)
    {
        splittable_count = 0;
        for (i=0; i<group_count; i++)
        {
            if (groups[i] > 1)
            {
                splittable_count++;
            }
        }
        if (splittable_count == 0)
        {
            break;
        }
        target_group = rnd(0, splittable_count-1);
        idx = -1;
        for (i=0; i<group_count; i++)
        {
            if (groups[i] > 1)
            {
                target_group--;
                if (target_group < 0)
                {
                    idx = i;
                    break;
                }
            }
        }
        if (idx < 0)
        {
            break;
        }
        left = rnd(1, groups[idx]-1);
        for (i=group_count; i>idx+1; i--)
        {
            groups[i] = groups[i-1];
        }
        groups[idx+1] = groups[idx] - left;
        groups[idx] = left;
        group_count++;
    }

    ind->cuts[0] = 0;
    pos = 0;
    for (i=0; i<group_count && i<n_slots; i++)
    {
        pos += groups[i];
        ind->cuts[i+1] = pos;
    }
    for (i=group_count+1; i<n_slots+1; i++)
    {
        ind->cuts[i] = n_clients;
    }
    ind->cuts[n_slots] = n_clients;
    normalize_cuts(ind);
}

static int cuts_capacity_valid (individual *ind)
{
    int slot;
    int pos;
    int customer;
    double load;

    for (slot=0; slot<n_slots; slot++)
    {
        load = 0.0;
        for (pos=ind->cuts[slot]; pos<ind->cuts[slot+1]; pos++)
        {
            customer = ind->perm[pos];
            load += problem.d[node_index(&problem, customer)];
        }
        if (load > problem.q)
        {
            return 0;
        }
    }
    return validate_individual_structure(ind);
}

static void normalize_cuts (individual *ind)
{
    int i;
    if (ind->cuts[0] != 0)
    {
        ind->cuts[0] = 0;
    }
    for (i=1; i<n_slots+1; i++)
    {
        if (ind->cuts[i] < ind->cuts[i-1])
        {
            ind->cuts[i] = ind->cuts[i-1];
        }
        if (ind->cuts[i] < 0)
        {
            ind->cuts[i] = 0;
        }
        if (ind->cuts[i] > n_clients)
        {
            ind->cuts[i] = n_clients;
        }
    }
    ind->cuts[n_slots] = n_clients;
}

static double planning_horizon_c (void)
{
    return max_ls() - min_li();
}

static double min_li (void)
{
    int i;
    double value;
    value = problem.LI[0];
    for (i=1; i<problem.nP; i++)
    {
        if (problem.LI[i] < value)
        {
            value = problem.LI[i];
        }
    }
    return value;
}

static double max_ls (void)
{
    int i;
    double value;
    value = problem.LS[0];
    for (i=1; i<problem.nP; i++)
    {
        if (problem.LS[i] > value)
        {
            value = problem.LS[i];
        }
    }
    return value;
}
