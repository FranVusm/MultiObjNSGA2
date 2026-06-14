/* NSGA-II routine (implementation of the 'main' function) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include <ctype.h>
# include <errno.h>
# ifdef _WIN32
# include <direct.h>
# define MKDIR(path) _mkdir(path)
# define PATH_SEP "\\"
# else
# include <sys/stat.h>
# define MKDIR(path) mkdir(path, 0777)
# define PATH_SEP "/"
# endif

# include "global.h"
# include "rand.h"

int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1;
int obj2;
int obj3;
int angle1;
int angle2;
problem_data problem;
int n_clients;
int n_slots;

# define DEFAULT_MULTI_RUNS 10
# define RESULTS_ROOT "resultados"

static char output_dir[1024] = ".";

static void sanitize_model_name (char *instance_route, char *model_name, int model_name_size)
{
    char *base;
    char *slash_back;
    char *slash_forward;
    char *dot;
    int i;
    int j;
    int length;

    slash_back = strrchr(instance_route, '\\');
    slash_forward = strrchr(instance_route, '/');
    base = instance_route;
    if (slash_back != NULL && slash_back + 1 > base)
    {
        base = slash_back + 1;
    }
    if (slash_forward != NULL && slash_forward + 1 > base)
    {
        base = slash_forward + 1;
    }
    dot = strrchr(base, '.');
    length = (dot != NULL && dot > base) ? (int)(dot - base) : (int)strlen(base);
    if (length <= 0)
    {
        strcpy(model_name, "modelo");
        return;
    }

    j = 0;
    for (i=0; i<length && j<model_name_size-1; i++)
    {
        unsigned char ch;
        ch = (unsigned char)base[i];
        if (isalnum(ch) || ch == '_' || ch == '-')
        {
            model_name[j++] = (char)ch;
        }
        else
        {
            model_name[j++] = '_';
        }
    }
    model_name[j] = '\0';
    if (j == 0)
    {
        strcpy(model_name, "modelo");
    }
}

static void create_output_directory (char *instance_route)
{
    char model_name[256];
    char timestamp[64];
    time_t now;
    struct tm *local_time;
    int attempt;

    sanitize_model_name(instance_route, model_name, sizeof(model_name));
    now = time(NULL);
    local_time = localtime(&now);
    if (local_time == NULL)
    {
        strcpy(timestamp, "unknown_time");
    }
    else
    {
        strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", local_time);
    }

    if (MKDIR(RESULTS_ROOT) != 0 && errno != EEXIST)
    {
        printf("\n Could not create results directory %s, hence exiting \n", RESULTS_ROOT);
        exit(1);
    }

    for (attempt=0; attempt<100; attempt++)
    {
        if (attempt == 0)
        {
            sprintf(output_dir, "%s%s%s_%s", RESULTS_ROOT, PATH_SEP, model_name, timestamp);
        }
        else
        {
            sprintf(output_dir, "%s%s%s_%s_%02d", RESULTS_ROOT, PATH_SEP, model_name, timestamp, attempt + 1);
        }
        if (MKDIR(output_dir) == 0)
        {
            printf("\n Output directory: %s\n", output_dir);
            return;
        }
        if (errno != EEXIST)
        {
            printf("\n Could not create output directory %s, hence exiting \n", output_dir);
            exit(1);
        }
    }
    printf("\n Could not create a unique output directory under %s, hence exiting \n", RESULTS_ROOT);
    exit(1);
}

static FILE* open_output_file (char *filename)
{
    FILE *file;
    char path[1536];
    char fallback[1600];
    sprintf(path, "%s%s%s", output_dir, PATH_SEP, filename);
    file = fopen(path, "w");
    if (file != NULL)
    {
        return file;
    }
    sprintf(fallback, "%s.new", path);
    printf("\n Warning: could not open %s, writing %s instead", path, fallback);
    return fopen(fallback, "w");
}

static double seed_for_run (double base_seed, int run_index)
{
    double value;
    value = fmod(base_seed + (0.097 * (double)run_index), 1.0);
    if (value <= 0.0)
    {
        value += 0.001;
    }
    if (value >= 1.0)
    {
        value = 0.999999;
    }
    return value;
}

static void write_population_headers (FILE *fpt1, FILE *fpt2, FILE *fpt3, FILE *fpt4)
{
    fprintf(fpt1,"# This file contains the data of initial population\n");
    fprintf(fpt2,"# This file contains the data of final population\n");
    fprintf(fpt3,"# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4,"# This file contains the data of all generations\n");
    fprintf(fpt1,"# F1, F2, perm[%d], cuts[%d], alpha[%d], constr_violation, rank, crowding_distance\n",n_clients,n_slots+1,n_slots);
    fprintf(fpt2,"# F1, F2, perm[%d], cuts[%d], alpha[%d], constr_violation, rank, crowding_distance\n",n_clients,n_slots+1,n_slots);
    fprintf(fpt3,"# F1, F2, perm[%d], cuts[%d], alpha[%d], constr_violation, rank, crowding_distance\n",n_clients,n_slots+1,n_slots);
    fprintf(fpt4,"# F1, F2, perm[%d], cuts[%d], alpha[%d], constr_violation, rank, crowding_distance\n",n_clients,n_slots+1,n_slots);
}

static int collect_final_candidates (population *pop, population *archive, int *archive_count, int archive_capacity, int *source_run, double *source_seed, int run_number, double run_seed)
{
    int i;
    int stored;
    stored = 0;
    for (i=0; i<popsize; i++)
    {
        if (pop->ind[i].rank == 1)
        {
            if (*archive_count < archive_capacity)
            {
                copy_ind(&(pop->ind[i]), &(archive->ind[*archive_count]));
                source_run[*archive_count] = run_number;
                source_seed[*archive_count] = run_seed;
                (*archive_count)++;
                stored++;
            }
        }
    }
    return stored;
}

static int run_nsga2_once (int run_number, double run_seed, population *archive, int *archive_count, int archive_capacity, int *source_run, double *source_seed, FILE *run_summary)
{
    int i;
    int stored;
    double elapsed;
    clock_t start_clock;
    clock_t end_clock;
    char filename[128];
    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

    sprintf(filename, "seed_%02d_initial_pop.out", run_number);
    fpt1 = open_output_file(filename);
    sprintf(filename, "seed_%02d_final_pop.out", run_number);
    fpt2 = open_output_file(filename);
    sprintf(filename, "seed_%02d_best_pop.out", run_number);
    fpt3 = open_output_file(filename);
    sprintf(filename, "seed_%02d_all_pop.out", run_number);
    fpt4 = open_output_file(filename);
    if (fpt1==NULL || fpt2==NULL || fpt3==NULL || fpt4==NULL)
    {
        printf("\n Could not open one or more seed output files, hence exiting \n");
        exit(1);
    }
    write_population_headers(fpt1, fpt2, fpt3, fpt4);

    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;
    seed = run_seed;
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);

    start_clock = clock();
    randomize();
    initialize_pop (parent_pop);
    printf("\n seed_%02d initialization done, now performing first generation", run_number);
    decode_pop(parent_pop);
    evaluate_pop (parent_pop);
    assign_rank_and_crowding_distance (parent_pop);
    report_pop (parent_pop, fpt1);
    fprintf(fpt4,"# gen = 1\n");
    report_pop(parent_pop,fpt4);
    printf("\n seed_%02d gen = 1", run_number);
    fflush(stdout);
    fflush(fpt1);
    fflush(fpt4);
    for (i=2; i<=ngen; i++)
    {
        selection (parent_pop, child_pop);
        mutation_pop (child_pop);
        decode_pop(child_pop);
        evaluate_pop(child_pop);
        merge (parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort (mixed_pop, parent_pop);
        fprintf(fpt4,"# gen = %d\n",i);
        report_pop(parent_pop,fpt4);
        fflush(fpt4);
        printf("\n seed_%02d gen = %d", run_number, i);
        fflush(stdout);
    }
    printf("\n seed_%02d generations finished, now reporting solutions\n", run_number);
    fflush(stdout);
    report_pop(parent_pop,fpt2);
    report_feasible(parent_pop,fpt3);
    stored = collect_final_candidates(parent_pop, archive, archive_count, archive_capacity, source_run, source_seed, run_number, run_seed);
    end_clock = clock();
    elapsed = (double)(end_clock - start_clock) / (double)CLOCKS_PER_SEC;
    fprintf(run_summary,"%d,%.12f,%d,%.6f\n", run_number, run_seed, stored, elapsed);

    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    fclose(fpt4);
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2*popsize);
    free (parent_pop);
    free (child_pop);
    free (mixed_pop);
    return stored;
}

static int same_final_representation (individual *a, individual *b)
{
    int i;
    for (i=0; i<n_clients; i++)
    {
        if (a->perm[i] != b->perm[i])
        {
            return 0;
        }
    }
    for (i=0; i<n_slots+1; i++)
    {
        if (a->cuts[i] != b->cuts[i])
        {
            return 0;
        }
    }
    for (i=0; i<n_slots; i++)
    {
        if (fabs(a->alpha[i] - b->alpha[i]) > 1.0e-9)
        {
            return 0;
        }
    }
    return 1;
}

static int dominates_objectives (individual *a, individual *b)
{
    int j;
    int better;
    better = 0;
    for (j=0; j<nobj; j++)
    {
        if (a->obj[j] > b->obj[j] + EPS)
        {
            return 0;
        }
        if (a->obj[j] + EPS < b->obj[j])
        {
            better = 1;
        }
    }
    return better;
}

static void mark_unique_and_pareto (population *archive, int count, int *is_unique, int *is_pareto)
{
    int i;
    int j;
    int has_feasible;
    has_feasible = 0;
    for (i=0; i<count; i++)
    {
        is_unique[i] = 1;
        is_pareto[i] = 0;
    }
    for (i=0; i<count; i++)
    {
        for (j=0; j<i; j++)
        {
            if (is_unique[j] && same_final_representation(&(archive->ind[i]), &(archive->ind[j])))
            {
                is_unique[i] = 0;
                break;
            }
        }
    }
    for (i=0; i<count; i++)
    {
        if (is_unique[i] && archive->ind[i].constr_violation == 0.0)
        {
            has_feasible = 1;
            break;
        }
    }
    for (i=0; i<count; i++)
    {
        if (!is_unique[i])
        {
            continue;
        }
        if (has_feasible && archive->ind[i].constr_violation != 0.0)
        {
            continue;
        }
        is_pareto[i] = 1;
        for (j=0; j<count; j++)
        {
            if (i == j || !is_unique[j])
            {
                continue;
            }
            if (has_feasible && archive->ind[j].constr_violation != 0.0)
            {
                continue;
            }
            if (dominates_objectives(&(archive->ind[j]), &(archive->ind[i])))
            {
                is_pareto[i] = 0;
                break;
            }
        }
    }
}

static void write_int_array (FILE *fpt, int *values, int count)
{
    int i;
    fprintf(fpt, "[");
    for (i=0; i<count; i++)
    {
        if (i > 0)
        {
            fprintf(fpt, ", ");
        }
        fprintf(fpt, "%d", values[i]);
    }
    fprintf(fpt, "]");
}

static void write_alpha_dict (FILE *fpt, individual *ind)
{
    int i;
    fprintf(fpt, "{");
    for (i=0; i<n_slots; i++)
    {
        if (i > 0)
        {
            fprintf(fpt, ", ");
        }
        fprintf(fpt, "(%d, %d): %.15g", problem.slots[i].trip, problem.slots[i].vehicle, ind->alpha[i]);
    }
    fprintf(fpt, "}");
}

static double slot_demand (individual *ind, int slot)
{
    int pos;
    double load;
    load = 0.0;
    for (pos=ind->cuts[slot]; pos<ind->cuts[slot+1]; pos++)
    {
        int idx;
        idx = node_index(&problem, ind->perm[pos]);
        if (idx >= 0)
        {
            load += problem.d[idx];
        }
    }
    return load;
}

static void write_customer_slice (FILE *fpt, individual *ind, int slot)
{
    int pos;
    fprintf(fpt, "[");
    for (pos=ind->cuts[slot]; pos<ind->cuts[slot+1]; pos++)
    {
        if (pos > ind->cuts[slot])
        {
            fprintf(fpt, ", ");
        }
        fprintf(fpt, "%d", ind->perm[pos]);
    }
    fprintf(fpt, "]");
}

static void write_active_trips_and_routes (FILE *fpt, individual *ind)
{
    int slot;
    int printed_empty;
    fprintf(fpt, "Active trips:\n");
    for (slot=0; slot<n_slots; slot++)
    {
        if (ind->cuts[slot] != ind->cuts[slot+1])
        {
            fprintf(fpt, "  (%d, %d): ", problem.slots[slot].trip, problem.slots[slot].vehicle);
            write_customer_slice(fpt, ind, slot);
            fprintf(fpt, " demand=%.1f\n", slot_demand(ind, slot));
        }
    }
    fprintf(fpt, "Empty slots: [");
    printed_empty = 0;
    for (slot=0; slot<n_slots; slot++)
    {
        if (ind->cuts[slot] == ind->cuts[slot+1])
        {
            if (printed_empty)
            {
                fprintf(fpt, ", ");
            }
            fprintf(fpt, "(%d, %d)", problem.slots[slot].trip, problem.slots[slot].vehicle);
            printed_empty = 1;
        }
    }
    fprintf(fpt, "]\n");
    fprintf(fpt, "Active routes:\n");
    for (slot=0; slot<n_slots; slot++)
    {
        int pos;
        if (ind->cuts[slot] != ind->cuts[slot+1])
        {
            fprintf(fpt, "  (%d, %d): [%d", problem.slots[slot].trip, problem.slots[slot].vehicle, problem.O);
            for (pos=ind->cuts[slot]; pos<ind->cuts[slot+1]; pos++)
            {
                fprintf(fpt, ", %d", ind->perm[pos]);
            }
            fprintf(fpt, ", %d]\n", problem.dummy_depot);
        }
    }
}

static void write_solution_detail (FILE *fpt, char *title, individual *ind, int source_run, double source_seed)
{
    fprintf(fpt, "--- %s ---\n", title);
    fprintf(fpt, "Source run: seed_%02d, seed=%.12f\n", source_run, source_seed);
    fprintf(fpt, "Final representation:\n");
    fprintf(fpt, "Raw F1: %.15g\n", ind->obj[0]);
    fprintf(fpt, "Raw F2: %.15g\n", ind->obj[1]);
    fprintf(fpt, "Penalized F1: %.15g\n", ind->obj[0]);
    fprintf(fpt, "Penalized F2: %.15g\n", ind->obj[1]);
    fprintf(fpt, "Feasible: %s\n", ind->constr_violation == 0.0 ? "True" : "False");
    fprintf(fpt, "Violations:\n");
    if (ind->constr_violation == 0.0)
    {
        fprintf(fpt, "  none\n");
    }
    else
    {
        fprintf(fpt, "  %.0f\n", -ind->constr_violation);
    }
    fprintf(fpt, "Permutation: ");
    write_int_array(fpt, ind->perm, n_clients);
    fprintf(fpt, "\nCuts: ");
    write_int_array(fpt, ind->cuts, n_slots+1);
    fprintf(fpt, "\nAlpha: ");
    write_alpha_dict(fpt, ind);
    fprintf(fpt, "\n");
    write_active_trips_and_routes(fpt, ind);
    fprintf(fpt, "\n");
}

static int best_pareto_index (population *archive, int count, int *is_pareto, int objective)
{
    int i;
    int best;
    best = -1;
    for (i=0; i<count; i++)
    {
        if (!is_pareto[i])
        {
            continue;
        }
        if (best < 0 || archive->ind[i].obj[objective] < archive->ind[best].obj[objective])
        {
            best = i;
        }
    }
    return best;
}

static int balanced_pareto_index (population *archive, int count, int *is_pareto)
{
    int i;
    int best;
    double min_f1;
    double max_f1;
    double min_f2;
    double max_f2;
    double best_score;
    min_f1 = INF;
    min_f2 = INF;
    max_f1 = -INF;
    max_f2 = -INF;
    for (i=0; i<count; i++)
    {
        if (!is_pareto[i])
        {
            continue;
        }
        if (archive->ind[i].obj[0] < min_f1) min_f1 = archive->ind[i].obj[0];
        if (archive->ind[i].obj[0] > max_f1) max_f1 = archive->ind[i].obj[0];
        if (archive->ind[i].obj[1] < min_f2) min_f2 = archive->ind[i].obj[1];
        if (archive->ind[i].obj[1] > max_f2) max_f2 = archive->ind[i].obj[1];
    }
    best = -1;
    best_score = INF;
    for (i=0; i<count; i++)
    {
        double nf1;
        double nf2;
        double score;
        if (!is_pareto[i])
        {
            continue;
        }
        nf1 = (max_f1 > min_f1) ? (archive->ind[i].obj[0] - min_f1) / (max_f1 - min_f1) : 0.0;
        nf2 = (max_f2 > min_f2) ? (archive->ind[i].obj[1] - min_f2) / (max_f2 - min_f2) : 0.0;
        score = sqrt((nf1 * nf1) + (nf2 * nf2));
        if (score < best_score)
        {
            best_score = score;
            best = i;
        }
    }
    return best;
}

static void write_pareto_points_csv (population *archive, int count, int *is_pareto, int *source_run, double *source_seed)
{
    int idx;
    int written;
    int *already_written;
    FILE *fpt;
    fpt = open_output_file("aggregate_pareto_points.csv");
    if (fpt == NULL)
    {
        return;
    }
    already_written = (int *)malloc(count * sizeof(int));
    if (already_written == NULL)
    {
        fclose(fpt);
        return;
    }
    for (idx=0; idx<count; idx++)
    {
        already_written[idx] = 0;
    }
    fprintf(fpt, "index,F1,F2,source_run,seed,feasible,perm,cuts,alpha\n");
    written = 0;
    for (;;)
    {
        int best;
        best = -1;
        for (idx=0; idx<count; idx++)
        {
            if (!is_pareto[idx])
            {
                continue;
            }
            if (already_written[idx])
            {
                continue;
            }
            if (best < 0 || archive->ind[idx].obj[0] < archive->ind[best].obj[0])
            {
                best = idx;
            }
        }
        if (best < 0)
        {
            break;
        }
        already_written[best] = 1;
        written++;
        fprintf(fpt, "%d,%.15g,%.15g,%d,%.12f,%s,", written, archive->ind[best].obj[0], archive->ind[best].obj[1], source_run[best], source_seed[best], archive->ind[best].constr_violation == 0.0 ? "True" : "False");
        write_int_array(fpt, archive->ind[best].perm, n_clients);
        fprintf(fpt, ",");
        write_int_array(fpt, archive->ind[best].cuts, n_slots+1);
        fprintf(fpt, ",");
        write_alpha_dict(fpt, &(archive->ind[best]));
        fprintf(fpt, "\n");
    }
    free(already_written);
    fclose(fpt);
}

static void write_aggregate_report (population *archive, int count, int *is_unique, int *is_pareto, int *source_run, double *source_seed, double total_elapsed)
{
    int i;
    int unique_count;
    int feasible_count;
    int pareto_count;
    int best_f1;
    int best_f2;
    int balanced;
    int printed;
    FILE *fpt;
    unique_count = 0;
    feasible_count = 0;
    pareto_count = 0;
    for (i=0; i<count; i++)
    {
        if (is_unique[i]) unique_count++;
        if (archive->ind[i].constr_violation == 0.0) feasible_count++;
        if (is_pareto[i]) pareto_count++;
    }
    fpt = open_output_file("aggregate_pareto_solutions.txt");
    if (fpt == NULL)
    {
        return;
    }
    fprintf(fpt, "TD-MT-GVRP Multi-Seed Aggregate Pareto Solutions\n\n");
    fprintf(fpt, "Aggregate summary\n");
    fprintf(fpt, "Stored candidate solutions: %d\n", count);
    fprintf(fpt, "Feasible candidate solutions: %d\n", feasible_count);
    fprintf(fpt, "Unique final representations: %d\n", unique_count);
    fprintf(fpt, "Global non-dominated solutions: %d\n", pareto_count);
    fprintf(fpt, "Total elapsed seconds: %.6f\n\n", total_elapsed);
    fprintf(fpt, "Candidate source summary\n");
    fprintf(fpt, "run,seed,source_file_prefix\n");
    for (i=0; i<count; i++)
    {
        if (is_unique[i])
        {
            fprintf(fpt, "%d,%.12f,seed_%02d\n", source_run[i], source_seed[i], source_run[i]);
        }
    }
    fprintf(fpt, "\nSelected solutions\n\n");
    best_f1 = best_pareto_index(archive, count, is_pareto, 0);
    best_f2 = best_pareto_index(archive, count, is_pareto, 1);
    balanced = balanced_pareto_index(archive, count, is_pareto);
    if (best_f1 >= 0)
    {
        write_solution_detail(fpt, "Best environmental solution (min F1)", &(archive->ind[best_f1]), source_run[best_f1], source_seed[best_f1]);
    }
    if (best_f2 >= 0)
    {
        write_solution_detail(fpt, "Best time/cost solution (min F2)", &(archive->ind[best_f2]), source_run[best_f2], source_seed[best_f2]);
    }
    if (balanced >= 0)
    {
        write_solution_detail(fpt, "Most balanced solution", &(archive->ind[balanced]), source_run[balanced], source_seed[balanced]);
    }
    fprintf(fpt, "Global Pareto front - non-dominated final representations\n\n");
    printed = 0;
    for (;;)
    {
        int best;
        best = -1;
        for (i=0; i<count; i++)
        {
            if (!is_pareto[i] || source_run[i] < 0)
            {
                continue;
            }
            if (best < 0 || archive->ind[i].obj[0] < archive->ind[best].obj[0])
            {
                best = i;
            }
        }
        if (best < 0)
        {
            break;
        }
        source_run[best] = -source_run[best];
        printed++;
        {
            char title[64];
            sprintf(title, "Pareto solution %d", printed);
            write_solution_detail(fpt, title, &(archive->ind[best]), -source_run[best], source_seed[best]);
        }
    }
    for (i=0; i<count; i++)
    {
        if (source_run[i] < 0)
        {
            source_run[i] = -source_run[i];
        }
    }
    fclose(fpt);
    write_pareto_points_csv(archive, count, is_pareto, source_run, source_seed);
}

static void generate_pareto_plot (void)
{
    char command[2300];
    int status;
    sprintf(command, "python plot_pareto.py \"%s\"", output_dir);
    printf("\n Generating Pareto plot\n");
    fflush(stdout);
    status = system(command);
    if (status != 0)
    {
        printf("\n Warning: could not generate pareto_front.png automatically. Run: %s", command);
    }
}

int main (int argc, char **argv)
{
    int i;
    FILE *fpt5;
    FILE *run_summary;

    char * instance_route;
    population *archive_pop;
    int *source_run;
    double *source_seed;
    int *is_unique;
    int *is_pareto;
    int archive_count;
    int archive_capacity;
    int runs;
    int stored_total;
    double base_seed;
    double total_elapsed;
    clock_t all_start;
    clock_t all_end;
    if (argc<7)
    {
        printf("\n Usage: ./nsga2r seed instance_route popsize ngen pcross pmut [runs]\n");
        printf("\n Example: ./nsga2r 0.123 ../../data/modelo_intermedio.dat 40 100 0.9 0.2\n");
        exit(1);
    }
    base_seed = (double)atof(argv[1]);
    if (base_seed<=0.0 || base_seed>=1.0){
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    runs = DEFAULT_MULTI_RUNS;
    if (argc >= 8)
    {
        runs = atoi(argv[7]);
        if (runs < 1)
        {
            printf("\n Number of runs must be positive, hence exiting \n");
            exit(1);
        }
    }
    instance_route = argv[2];
    create_output_directory(instance_route);
    fpt5 = open_output_file("params.out");
    run_summary = open_output_file("multi_seed_run_summary.csv");
    if (fpt5==NULL || run_summary==NULL)
    {
        if (fpt5==NULL) printf("\n Could not open params.out");
        if (run_summary==NULL) printf("\n Could not open multi_seed_run_summary.csv");
        printf("\n Could not open one or more output files, hence exiting \n");
        exit(1);
    }
    fprintf(fpt5,"# This file contains information about inputs as read by the program\n");
    fprintf(run_summary,"run,seed,rank1_candidates,elapsed_seconds\n");

    if (!read_ampl_data(instance_route, &problem))
    {
        printf("\n Could not load AMPL instance, hence exiting \n");
        exit(1);
    }
    n_clients = problem.nCustomers;
    n_slots = problem.nSlots;
    nreal = 0;
    nbin = 0;
    ncon = 0;
    nobj = 2;

    popsize = atoi(argv[3]);
    if (popsize<4 || (popsize%4)!= 0){
        printf("\n population size read is : %d",popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }
    ngen = atoi(argv[4]);
    if (ngen<1){
        printf("\n number of generations read is : %d",ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }
    if (n_clients<1 || n_slots<1){
        printf("\n Instance must contain positive-demand customers and trip slots, hence exiting \n");
        exit (1);
    }
    /* Setear en lectura de instancia --
    printf("\n Enter the number of binary variables : ");
    scanf("%d",&nbin);
    if (nbin<0)
    {
        printf ("\n number of binary variables entered is : %d",nbin);
        printf ("\n Wrong number of binary variables entered, hence exiting \n");
        exit(1);
    }
    if (nbin != 0)
    {
        nbits = (int *)malloc(nbin*sizeof(int));
        min_binvar = (double *)malloc(nbin*sizeof(double));
        max_binvar = (double *)malloc(nbin*sizeof(double));
        for (i=0; i<nbin; i++)
        {
            printf ("\n Enter the number of bits for binary variable %d : ",i+1);
            scanf ("%d",&nbits[i]);
            if (nbits[i] < 1)
            {
                printf("\n Wrong number of bits for binary variable entered, hence exiting");
                exit(1);
            }
            printf ("\n Enter the lower limit of binary variable %d : ",i+1);
            scanf ("%lf",&min_binvar[i]);
            printf ("\n Enter the upper limit of binary variable %d : ",i+1);
            scanf ("%lf",&max_binvar[i]);
            if (max_binvar[i] <= min_binvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
                exit(1);
            }
        }
    */
    pcross_bin = atof (argv[5]);
    if (pcross_bin<0.0 || pcross_bin>1.0){
        printf("\n Probability of crossover entered is : %e",pcross_bin);
        printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
        exit (1);
    }
    pmut_bin = atof (argv[6]);
    if (pmut_bin<0.0 || pmut_bin>1.0){
        printf("\n Probability of mutation entered is : %e",pmut_bin);
        printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
        exit (1);
    }

    printf("\n Input data successfully entered, now performing initialization \n");
    fprintf(fpt5,"\n Population size = %d",popsize);
    fprintf(fpt5,"\n Number of generations = %d",ngen);
    fprintf(fpt5,"\n Number of objective functions = %d",nobj);
    fprintf(fpt5,"\n Instance route = %s", instance_route);
    fprintf(fpt5,"\n Output directory = %s", output_dir);
    fprintf(fpt5,"\n Number of nodes = %d", problem.nN);
    fprintf(fpt5,"\n Number of customers = %d", n_clients);
    fprintf(fpt5,"\n Number of periods = %d", problem.nP);
    fprintf(fpt5,"\n Number of vehicles = %d", problem.nK);
    fprintf(fpt5,"\n Number of trips = %d", problem.nV);
    fprintf(fpt5,"\n Number of trip slots = %d", n_slots);
    fprintf(fpt5,"\n Origin = %d", problem.O);
    fprintf(fpt5,"\n Dummy depot = %d", problem.dummy_depot);
    fprintf(fpt5,"\n Vehicle capacity = %e", problem.q);
    /*fprintf(fpt5,"\n Number of constraints = %d",ncon);
    fprintf(fpt5,"\n Number of real variables = %d",nreal);
    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,min_realvar[i]);
            fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,max_realvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of real variable = %e",pcross_real);
        fprintf(fpt5,"\n Probability of mutation of real variable = %e",pmut_real);
        fprintf(fpt5,"\n Distribution index for crossover = %e",eta_c);
        fprintf(fpt5,"\n Distribution index for mutation = %e",eta_m);
    }*/
    fprintf(fpt5,"\n Probability of permutation crossover = %e",pcross_bin);
    fprintf(fpt5,"\n Probability of chromosome mutation = %e",pmut_bin);
    fprintf(fpt5,"\n Base seed for random number generator = %e",base_seed);
    fprintf(fpt5,"\n Multi-seed runs = %d",runs);
    bitlength = 0;
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            bitlength += nbits[i];
        }
    }
    archive_capacity = runs * popsize;
    archive_count = 0;
    stored_total = 0;
    archive_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop(archive_pop, archive_capacity);
    source_run = (int *)malloc(archive_capacity * sizeof(int));
    source_seed = (double *)malloc(archive_capacity * sizeof(double));
    is_unique = (int *)malloc(archive_capacity * sizeof(int));
    is_pareto = (int *)malloc(archive_capacity * sizeof(int));
    all_start = clock();
    for (i=0; i<runs; i++)
    {
        double run_seed;
        int stored;
        run_seed = seed_for_run(base_seed, i);
        printf("\n\n Multi-seed run %d/%d with seed %.12f", i+1, runs, run_seed);
        stored = run_nsga2_once(i+1, run_seed, archive_pop, &archive_count, archive_capacity, source_run, source_seed, run_summary);
        stored_total += stored;
    }
    all_end = clock();
    total_elapsed = (double)(all_end - all_start) / (double)CLOCKS_PER_SEC;
    mark_unique_and_pareto(archive_pop, archive_count, is_unique, is_pareto);
    write_aggregate_report(archive_pop, archive_count, is_unique, is_pareto, source_run, source_seed, total_elapsed);
    generate_pareto_plot();
    if (nreal!=0)
    {
        fprintf(fpt5,"\n Number of crossover of real variable = %d",nrealcross);
        fprintf(fpt5,"\n Number of mutation of real variable = %d",nrealmut);
    }
    if (nbin!=0)
    {
        fprintf(fpt5,"\n Number of crossover of binary variable = %d",nbincross);
        fprintf(fpt5,"\n Number of mutation of binary variable = %d",nbinmut);
    }
    fprintf(fpt5,"\n Stored feasible rank-1 candidates = %d",stored_total);
    fflush(stdout);
    fflush(fpt5);
    fflush(run_summary);
    fclose(fpt5);
    fclose(run_summary);
    if (nbin!=0)
    {
        free (min_binvar);
        free (max_binvar);
        free (nbits);
    }
    deallocate_memory_pop (archive_pop, archive_capacity);
    free (archive_pop);
    free (source_run);
    free (source_seed);
    free (is_unique);
    free (is_pareto);
    free_problem_data(&problem);
    printf("\n Routine successfully exited \n");
    return (0);
}
