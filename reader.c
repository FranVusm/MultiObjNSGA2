/* AMPL .dat reader for the TD-MT-GVRP NSGA-II representation */

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>

# include "global.h"

static char* read_clean_file (char *file_path);
static char* find_declaration (char *text, char *kind, char *name);
static char* declaration_body (char *decl);
static char* declaration_end (char *decl);
static int parse_int_set (char *text, char *name, int **values, int *count);
static int parse_scalar_int (char *text, char *name, int *value);
static int parse_scalar_double (char *text, char *name, double *value);
static int parse_1d_param (char *text, char *name, problem_data *data, double *values, int by_period);
static int parse_2d_pairs (char *text, char *name, problem_data *data, double *values);
static int parse_3d_matrix (char *text, char *name, problem_data *data, double *values);
static int parse_numbers (char *start, char *end, double **values, int *count);
static int next_number (char **cursor, char *end, double *value);
static int compare_ints (const void *a, const void *b);
static void build_customers_and_slots (problem_data *data);
static void zero_problem_data (problem_data *data);
static double* allocate_double_array (int size);
static int index_3d (problem_data *data, int from_index, int to_index, int period_index_value);
static int index_2d (problem_data *data, int node_index_value, int period_index_value);

int read_ampl_data (char *file_path, problem_data *data)
{
    char *text;
    int total_3d;
    int total_2d;

    zero_problem_data(data);
    text = read_clean_file(file_path);
    if (text == NULL)
    {
        printf("Could not read AMPL data file: %s\n", file_path);
        return 0;
    }

    if (!parse_int_set(text, "N", &data->N, &data->nN)) return 0;
    if (!parse_int_set(text, "P", &data->P, &data->nP)) return 0;
    if (!parse_int_set(text, "K", &data->K, &data->nK)) return 0;
    if (!parse_int_set(text, "V", &data->V, &data->nV)) return 0;

    qsort(data->N, data->nN, sizeof(int), compare_ints);
    qsort(data->P, data->nP, sizeof(int), compare_ints);
    qsort(data->K, data->nK, sizeof(int), compare_ints);
    qsort(data->V, data->nV, sizeof(int), compare_ints);

    if (!parse_scalar_int(text, "O", &data->O)) return 0;
    if (!parse_scalar_double(text, "q", &data->q)) return 0;
    data->dummy_depot = data->N[data->nN-1];

    data->d = allocate_double_array(data->nN);
    data->LI = allocate_double_array(data->nP);
    data->LS = allocate_double_array(data->nP);
    total_3d = data->nN * data->nN * data->nP;
    total_2d = data->nN * data->nP;
    data->e = allocate_double_array(total_3d);
    data->g = allocate_double_array(total_3d);
    data->T = allocate_double_array(total_3d);
    data->ee = allocate_double_array(total_2d);
    data->gg = allocate_double_array(total_2d);
    data->tt = allocate_double_array(total_2d);

    if (!parse_1d_param(text, "d", data, data->d, 0)) return 0;
    if (!parse_1d_param(text, "LI", data, data->LI, 1)) return 0;
    if (!parse_1d_param(text, "LS", data, data->LS, 1)) return 0;
    if (!parse_3d_matrix(text, "e", data, data->e)) return 0;
    if (!parse_3d_matrix(text, "g", data, data->g)) return 0;
    if (!parse_3d_matrix(text, "T", data, data->T)) return 0;
    if (!parse_2d_pairs(text, "ee", data, data->ee)) return 0;
    if (!parse_2d_pairs(text, "gg", data, data->gg)) return 0;
    if (!parse_2d_pairs(text, "tt", data, data->tt)) return 0;

    build_customers_and_slots(data);
    free(text);
    return 1;
}

void free_problem_data (problem_data *data)
{
    free(data->N);
    free(data->P);
    free(data->K);
    free(data->V);
    free(data->d);
    free(data->LI);
    free(data->LS);
    free(data->e);
    free(data->ee);
    free(data->g);
    free(data->gg);
    free(data->T);
    free(data->tt);
    free(data->customers);
    free(data->slots);
    zero_problem_data(data);
}

int node_index (problem_data *data, int node)
{
    int i;
    for (i=0; i<data->nN; i++)
    {
        if (data->N[i] == node)
        {
            return i;
        }
    }
    return -1;
}

int period_index (problem_data *data, int period)
{
    int i;
    for (i=0; i<data->nP; i++)
    {
        if (data->P[i] == period)
        {
            return i;
        }
    }
    return -1;
}

double data_3d (double *values, problem_data *data, int from, int to, int period)
{
    int from_idx;
    int to_idx;
    int p_idx;
    from_idx = node_index(data, from);
    to_idx = node_index(data, to);
    p_idx = period_index(data, period);
    if (from_idx < 0 || to_idx < 0 || p_idx < 0)
    {
        return INF;
    }
    return values[index_3d(data, from_idx, to_idx, p_idx)];
}

double data_2d (double *values, problem_data *data, int node, int period)
{
    int node_idx;
    int p_idx;
    node_idx = node_index(data, node);
    p_idx = period_index(data, period);
    if (node_idx < 0 || p_idx < 0)
    {
        return INF;
    }
    return values[index_2d(data, node_idx, p_idx)];
}

static char* read_clean_file (char *file_path)
{
    FILE *file;
    long size;
    char *raw;
    char *clean;
    long i;
    long j;
    int in_comment;

    file = fopen(file_path, "rb");
    if (file == NULL)
    {
        return NULL;
    }
    fseek(file, 0, SEEK_END);
    size = ftell(file);
    rewind(file);
    raw = (char *)malloc((size + 1) * sizeof(char));
    clean = (char *)malloc((size + 1) * sizeof(char));
    if (raw == NULL || clean == NULL)
    {
        fclose(file);
        free(raw);
        free(clean);
        return NULL;
    }
    fread(raw, 1, size, file);
    raw[size] = '\0';
    fclose(file);

    in_comment = 0;
    j = 0;
    for (i=0; i<size; i++)
    {
        if (raw[i] == '#')
        {
            in_comment = 1;
        }
        if (raw[i] == '\n' || raw[i] == '\r')
        {
            in_comment = 0;
            clean[j++] = raw[i];
        }
        else if (!in_comment)
        {
            clean[j++] = raw[i];
        }
    }
    clean[j] = '\0';
    free(raw);
    return clean;
}

static char* find_declaration (char *text, char *kind, char *name)
{
    char *cursor;
    int kind_len;
    int name_len;

    cursor = text;
    kind_len = strlen(kind);
    name_len = strlen(name);
    while ((cursor = strstr(cursor, kind)) != NULL)
    {
        char *p;
        p = cursor + kind_len;
        if (cursor != text && (isalnum((unsigned char)cursor[-1]) || cursor[-1] == '_'))
        {
            cursor++;
            continue;
        }
        while (*p && isspace((unsigned char)*p)) p++;
        if (strncmp(p, name, name_len) == 0 && !(isalnum((unsigned char)p[name_len]) || p[name_len] == '_'))
        {
            return cursor;
        }
        cursor++;
    }
    return NULL;
}

static char* declaration_body (char *decl)
{
    char *body;
    body = strstr(decl, ":=");
    if (body == NULL)
    {
        return NULL;
    }
    return body + 2;
}

static char* declaration_end (char *decl)
{
    return strchr(decl, ';');
}

static int parse_int_set (char *text, char *name, int **values, int *count)
{
    char *decl;
    char *body;
    char *end;
    double *numbers;
    int number_count;
    int i;

    decl = find_declaration(text, "set", name);
    if (decl == NULL) return 0;
    body = declaration_body(decl);
    end = declaration_end(decl);
    if (body == NULL || end == NULL) return 0;
    if (!parse_numbers(body, end, &numbers, &number_count)) return 0;
    *values = (int *)malloc(number_count * sizeof(int));
    *count = number_count;
    for (i=0; i<number_count; i++)
    {
        (*values)[i] = (int)numbers[i];
    }
    free(numbers);
    return 1;
}

static int parse_scalar_int (char *text, char *name, int *value)
{
    double scalar;
    if (!parse_scalar_double(text, name, &scalar)) return 0;
    *value = (int)scalar;
    return 1;
}

static int parse_scalar_double (char *text, char *name, double *value)
{
    char *decl;
    char *body;
    char *end;
    double parsed;

    decl = find_declaration(text, "param", name);
    if (decl == NULL) return 0;
    body = declaration_body(decl);
    end = declaration_end(decl);
    if (body == NULL || end == NULL) return 0;
    if (!next_number(&body, end, &parsed)) return 0;
    *value = parsed;
    return 1;
}

static int parse_1d_param (char *text, char *name, problem_data *data, double *values, int by_period)
{
    char *decl;
    char *body;
    char *end;
    double *numbers;
    int number_count;
    int i;

    decl = find_declaration(text, "param", name);
    if (decl == NULL) return 0;
    body = declaration_body(decl);
    end = declaration_end(decl);
    if (body == NULL || end == NULL) return 0;
    if (!parse_numbers(body, end, &numbers, &number_count)) return 0;
    for (i=0; i+1<number_count; i+=2)
    {
        int key;
        int idx;
        key = (int)numbers[i];
        idx = by_period ? period_index(data, key) : node_index(data, key);
        if (idx >= 0)
        {
            values[idx] = numbers[i+1];
        }
    }
    free(numbers);
    return 1;
}

static int parse_2d_pairs (char *text, char *name, problem_data *data, double *values)
{
    char *decl;
    char *body;
    char *end;
    double *numbers;
    int number_count;
    int i;

    decl = find_declaration(text, "param", name);
    if (decl == NULL) return 0;
    body = declaration_body(decl);
    end = declaration_end(decl);
    if (body == NULL || end == NULL) return 0;
    if (!parse_numbers(body, end, &numbers, &number_count)) return 0;
    for (i=0; i+2<number_count; i+=3)
    {
        int node;
        int period;
        int n_idx;
        int p_idx;
        node = (int)numbers[i];
        period = (int)numbers[i+1];
        n_idx = node_index(data, node);
        p_idx = period_index(data, period);
        if (n_idx >= 0 && p_idx >= 0)
        {
            values[index_2d(data, n_idx, p_idx)] = numbers[i+2];
        }
    }
    free(numbers);
    return 1;
}

static int parse_3d_matrix (char *text, char *name, problem_data *data, double *values)
{
    char *decl;
    char *end_decl;
    char *cursor;

    decl = find_declaration(text, "param", name);
    if (decl == NULL) return 0;
    end_decl = declaration_end(decl);
    if (end_decl == NULL) return 0;
    cursor = decl;

    while ((cursor = strstr(cursor, "[*,*,")) != NULL && cursor < end_decl)
    {
        char *period_start;
        char *header_start;
        char *assign;
        char *block_start;
        char *next_block;
        char *block_end;
        double period_value;
        double *columns;
        int column_count;
        int p_idx;
        char *row_cursor;

        period_start = cursor + 5;
        if (!next_number(&period_start, end_decl, &period_value)) return 0;
        p_idx = period_index(data, (int)period_value);
        header_start = strchr(cursor, ':');
        assign = strstr(cursor, ":=");
        if (header_start == NULL || assign == NULL || assign > end_decl) return 0;
        if (!parse_numbers(header_start + 1, assign, &columns, &column_count)) return 0;
        block_start = assign + 2;
        next_block = strstr(block_start, "[*,*,");
        block_end = (next_block != NULL && next_block < end_decl) ? next_block : end_decl;
        row_cursor = block_start;
        while (row_cursor < block_end)
        {
            char *line_end;
            double *row_numbers;
            int row_count;
            int row_node;
            int row_idx;
            int j;

            line_end = strchr(row_cursor, '\n');
            if (line_end == NULL || line_end > block_end)
            {
                line_end = block_end;
            }
            if (parse_numbers(row_cursor, line_end, &row_numbers, &row_count) && row_count > 1)
            {
                row_node = (int)row_numbers[0];
                row_idx = node_index(data, row_node);
                for (j=1; j<row_count && j<=column_count; j++)
                {
                    int col_idx;
                    col_idx = node_index(data, (int)columns[j-1]);
                    if (row_idx >= 0 && col_idx >= 0 && p_idx >= 0)
                    {
                        values[index_3d(data, row_idx, col_idx, p_idx)] = row_numbers[j];
                    }
                }
                free(row_numbers);
            }
            row_cursor = line_end + 1;
        }
        free(columns);
        cursor = block_end;
    }
    return 1;
}

static int parse_numbers (char *start, char *end, double **values, int *count)
{
    int capacity;
    char *cursor;
    double value;

    capacity = 16;
    *count = 0;
    *values = (double *)malloc(capacity * sizeof(double));
    if (*values == NULL) return 0;
    cursor = start;
    while (next_number(&cursor, end, &value))
    {
        if (*count >= capacity)
        {
            double *resized;
            capacity *= 2;
            resized = (double *)realloc(*values, capacity * sizeof(double));
            if (resized == NULL)
            {
                free(*values);
                *values = NULL;
                *count = 0;
                return 0;
            }
            *values = resized;
        }
        (*values)[*count] = value;
        (*count)++;
    }
    return 1;
}

static int next_number (char **cursor, char *end, double *value)
{
    char *p;
    char *after;
    p = *cursor;
    while (p < end)
    {
        if (isdigit((unsigned char)*p) || *p == '-' || *p == '+' || *p == '.')
        {
            *value = strtod(p, &after);
            if (after != p)
            {
                *cursor = after;
                return 1;
            }
        }
        p++;
    }
    *cursor = p;
    return 0;
}

static int compare_ints (const void *a, const void *b)
{
    int left;
    int right;
    left = *((const int *)a);
    right = *((const int *)b);
    return (left > right) - (left < right);
}

static void build_customers_and_slots (problem_data *data)
{
    int i;
    int k;
    int v;
    int count;

    data->customers = (int *)malloc(data->nN * sizeof(int));
    data->nCustomers = 0;
    for (i=0; i<data->nN; i++)
    {
        if (data->d[i] > 0.0)
        {
            data->customers[data->nCustomers] = data->N[i];
            data->nCustomers++;
        }
    }
    qsort(data->customers, data->nCustomers, sizeof(int), compare_ints);

    data->nSlots = data->nK * data->nV;
    data->slots = (trip_slot *)malloc(data->nSlots * sizeof(trip_slot));
    count = 0;
    for (v=0; v<data->nV; v++)
    {
        for (k=0; k<data->nK; k++)
        {
            data->slots[count].trip = data->V[v];
            data->slots[count].vehicle = data->K[k];
            count++;
        }
    }
}

static void zero_problem_data (problem_data *data)
{
    memset(data, 0, sizeof(problem_data));
}

static double* allocate_double_array (int size)
{
    double *values;
    int i;
    values = (double *)malloc(size * sizeof(double));
    if (values == NULL)
    {
        return NULL;
    }
    for (i=0; i<size; i++)
    {
        values[i] = 0.0;
    }
    return values;
}

static int index_3d (problem_data *data, int from_index, int to_index, int period_index_value)
{
    return (period_index_value * data->nN * data->nN) + (from_index * data->nN) + to_index;
}

static int index_2d (problem_data *data, int node_index_value, int period_index_value)
{
    return (period_index_value * data->nN) + node_index_value;
}
