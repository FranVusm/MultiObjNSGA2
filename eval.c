/* Decoder and objective evaluation for perm/cuts/alpha chromosomes */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

static int locate_period_c (double time_value);
static double min_li_c (void);
static double max_ls_c (void);
static int first_trip (void);
static void active_trip_prefix_violations (individual *ind, int *violations);

void evaluate_pop (population *pop)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        evaluate_ind (&(pop->ind[i]));
    }
    return;
}

void evaluate_ind (individual *ind)
{
    int slot;
    int pos;
    int violations;
    double F1;
    double F2;
    double start_time;
    double planning_end;
    int vehicle_index;

    violations = 0;
    F1 = 0.0;
    F2 = 0.0;

    if (!validate_individual_structure(ind))
    {
        violations++;
        repair_cuts_capacity_aware(ind);
    }

    active_trip_prefix_violations(ind, &violations);

    for (slot=0; slot<n_slots; slot++)
    {
        double load;
        load = 0.0;
        for (pos=ind->cuts[slot]; pos<ind->cuts[slot+1]; pos++)
        {
            load += problem.d[node_index(&problem, ind->perm[pos])];
        }
        if (load > problem.q)
        {
            violations++;
        }
    }

    start_time = min_li_c();
    planning_end = max_ls_c();
    for (vehicle_index=0; vehicle_index<problem.nK; vehicle_index++)
    {
        double previous_vehicle_available_time;
        previous_vehicle_available_time = start_time;
        for (slot=0; slot<n_slots; slot++)
        {
            int from;
            int to;
            int route_pos;
            int route_size;
            int period;
            int period_idx;
            int route_completed;
            double current_time;

            if (problem.slots[slot].vehicle != problem.K[vehicle_index])
            {
                continue;
            }
            if (ind->cuts[slot] == ind->cuts[slot+1])
            {
                continue;
            }

            if (problem.slots[slot].trip == first_trip())
            {
                current_time = start_time + ind->alpha[slot];
            }
            else
            {
                current_time = previous_vehicle_available_time + ind->alpha[slot];
            }

            route_size = (ind->cuts[slot+1] - ind->cuts[slot]) + 2;
            route_completed = 1;
            for (route_pos=0; route_pos<route_size-1; route_pos++)
            {
                double depart_time;
                double emission;
                double cost;
                double travel_time;
                double service_time;

                if (route_pos == 0)
                {
                    from = problem.O;
                }
                else
                {
                    from = ind->perm[ind->cuts[slot] + route_pos - 1];
                }
                if (route_pos == route_size-2)
                {
                    to = problem.dummy_depot;
                }
                else
                {
                    to = ind->perm[ind->cuts[slot] + route_pos];
                }

                depart_time = current_time;
                period = locate_period_c(depart_time);
                if (period < 0)
                {
                    violations++;
                    route_completed = 0;
                    break;
                }
                period_idx = period_index(&problem, period);

                emission = data_3d(problem.e, &problem, from, to, period) + data_2d(problem.ee, &problem, to, period);
                cost = data_3d(problem.g, &problem, from, to, period) + data_2d(problem.gg, &problem, to, period);
                travel_time = data_3d(problem.T, &problem, from, to, period);
                service_time = data_2d(problem.tt, &problem, to, period);
                if (emission >= INF || cost >= INF || travel_time >= INF || service_time >= INF)
                {
                    violations++;
                    route_completed = 0;
                    break;
                }

                F1 += emission;
                F2 += cost;
                current_time = current_time + travel_time + service_time;
                if (current_time > problem.LS[period_idx])
                {
                    violations++;
                }
            }
            if (route_completed)
            {
                previous_vehicle_available_time = current_time;
                if (previous_vehicle_available_time > planning_end)
                {
                    violations++;
                }
            }
        }
    }

    ind->obj[0] = F1;
    ind->obj[1] = F2;
    ind->constr_violation = (violations == 0) ? 0.0 : -(double)violations;
    return;
}

static void active_trip_prefix_violations (individual *ind, int *violations)
{
    int vehicle_index;
    int slot;
    for (vehicle_index=0; vehicle_index<problem.nK; vehicle_index++)
    {
        int empty_seen;
        empty_seen = 0;
        for (slot=0; slot<n_slots; slot++)
        {
            if (problem.slots[slot].vehicle != problem.K[vehicle_index])
            {
                continue;
            }
            if (ind->cuts[slot] == ind->cuts[slot+1])
            {
                empty_seen = 1;
            }
            else if (empty_seen)
            {
                (*violations)++;
            }
        }
    }
}

static int locate_period_c (double time_value)
{
    int i;
    for (i=0; i<problem.nP; i++)
    {
        if (problem.LI[i] <= time_value && time_value <= problem.LS[i])
        {
            return problem.P[i];
        }
    }
    return -1;
}

static int first_trip (void)
{
    return problem.V[0];
}

static double min_li_c (void)
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

static double max_ls_c (void)
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
