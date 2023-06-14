#define _CRT_SECURE_NO_WARNINGS
#include <mpi.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <random>

#define city_total 280
#define testcase_name "a280.txt"
#define population_size 2000
#define iter_total 10000
#define iter_one 200
double city_location[city_total][2];
double distance_matrix[city_total][city_total];

bool caculation_end = false;

#define city_distance(idx, idy) distance_matrix[(idx)][(idy)]

struct chromosome {
    int path[city_total];
    double path_length;
};
chromosome *population;
chromosome *buffer;

int rand_between(int start, int end, std::mt19937& rng) {
    return rng() % (end - start) + start;
}
double genrand_real(std::mt19937& rng) { return rng() * (1.0 / 4294967296.0); }
void citys_init() {
    FILE* fp = fopen(testcase_name, "r");
    for (int i = 0; i < city_total; i++) {
        fscanf(fp, "%lf %lf", &city_location[i][0], &city_location[i][1]);
    }
    fclose(fp);
    for (int i = 0; i < city_total; i++) {
        for (int j = 0; j <= i; j++) {
            if (i == j) {
                distance_matrix[i][j] = 0;
            }
            double x = city_location[i][0] - city_location[j][0];
            double y = city_location[i][1] - city_location[j][1];
            double t = sqrt(x * x + y * y);
            city_distance(i, j) = city_distance(j, i) = t;
        }
    }
}

void path_random_init(int* path, std::mt19937& rng) {
    for (int i = 0; i < city_total; i++) {
        path[i] = i;
    }
    std::shuffle(path, path + city_total, rng);
}

void path_duplicate(int* to, const int* from) {
    memcpy(to, from, sizeof(int) * city_total);
}
double path_length(int* path) {
    double length = city_distance(path[0], path[city_total - 1]);
    for (int i = 0; i < city_total - 1; i++) {
        length += city_distance(path[i], path[i + 1]);
    }
    return length;
}

// compare by fitness (path length)
bool chromosome_cmp_func(const chromosome& a, const chromosome& b) {
    return a.path_length < b.path_length;
}

void population_init(chromosome* p, std::mt19937& rng) {
    for (int i = 0; i < population_size; i++) {
        path_random_init(p[i].path, rng);
        p[i].path_length = path_length(p[i].path);
    }
}

void population_sort(chromosome* p) {
    std::sort(p, p + population_size, chromosome_cmp_func);
}

void population_select(chromosome* p, chromosome* buffer, int elite_size, std::mt19937& rng) {
    int random_city_idx[population_size];
    for (int i = 0; i < population_size; i++) {
        random_city_idx[i] = i;
    }
    std::shuffle(random_city_idx, random_city_idx + population_size, rng);
    for (int i = elite_size; i < population_size; i++) {
        buffer[i] = p[random_city_idx[i]];
    }
    memcpy(&p[elite_size], &buffer[elite_size],
        sizeof(struct chromosome) * (population_size - elite_size));
}

void population_breed(chromosome* p, chromosome* buffer, int elite_size, std::mt19937& rng) {
    int cross_path[city_total];
    int diff_path[city_total];
    int new_path1[city_total], new_path2[city_total];
    int path_in_cross_path_flag[city_total];

    for (int i = elite_size; i < population_size; i++) {
        int idx, idy;
        do {
            idx = rand_between(0, population_size, rng);
            idy = rand_between(0, population_size, rng);
        } while (idx == idy);
        memset(path_in_cross_path_flag, 0, sizeof(path_in_cross_path_flag));
        int x, y;
        int cross_length;
        do {
            x = rand_between(0, city_total, rng);
            y = rand_between(0, city_total, rng);
        } while (x >= y);
        cross_length = y - x;
        memcpy(cross_path, p[idx].path, sizeof(int) * cross_length);
        for (int i = 0; i < cross_length; i++) {
            path_in_cross_path_flag[cross_path[i]] = 1;
        }
        int diff_length = 0;
        for (int i = 0; i < city_total; i++) {
            if (!path_in_cross_path_flag[p[idy].path[i]]) {
                diff_path[diff_length++] = p[idy].path[i];
            }
        }
        memcpy(new_path1, cross_path, sizeof(int) * cross_length);
        memcpy(new_path1 + cross_length, diff_path, sizeof(int) * diff_length);
        memcpy(new_path2, diff_path, sizeof(int) * diff_length);
        memcpy(new_path2 + diff_length, cross_path, sizeof(int) * cross_length);
        double pl1 = path_length(new_path1);
        double pl2 = path_length(new_path2);
        int* np;
        double nl;
        if (pl1 < pl2) {
            np = new_path1;
            nl = pl1;
        }
        else {
            np = new_path2;
            nl = pl2;
        }
        path_duplicate(buffer[i].path, np);
        buffer[i].path_length = nl;
    }
    memcpy(&p[elite_size], &buffer[elite_size], sizeof(chromosome) * (population_size - elite_size));
}

void population_mutate(chromosome* p, double rate, std::mt19937& rng) {
    for (int i = 1; i < population_size; i++) {
        if (genrand_real(rng) > rate) {
            continue;
        }
        int x = rand_between(0, city_total, rng);
        int y = rand_between(0, city_total, rng);
        int t = p[i].path[x];
        p[i].path[x] = p[i].path[y];
        p[i].path[y] = t;
        p[i].path_length = path_length(p[i].path);
    }
}

void population_update(chromosome* p, chromosome* buffer, std::mt19937& rng) {
    population_sort(p);
    population_select(p, buffer, population_size / 5, rng);
    population_breed(p, buffer, population_size / 5, rng);
    population_mutate(p, 0.1, rng);
}
#define TSP_GA_MASTER_PROCESS 0
int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
        processor_name, world_rank, world_size);
    fflush(stdout);

    std::mt19937 rng = std::mt19937(std::random_device{}());
    population = new chromosome[world_size * population_size];
    buffer = new chromosome[world_size * population_size];

    if (world_rank == TSP_GA_MASTER_PROCESS)
    {
        citys_init();
        for (int i = 0; i < world_size; i++)
        {
            population_init(&population[i * population_size], rng);
        }
    }

    MPI_Bcast(distance_matrix, sizeof(distance_matrix), MPI_BYTE, TSP_GA_MASTER_PROCESS, MPI_COMM_WORLD);
    MPI_Bcast(population, world_size * population_size * sizeof(chromosome), MPI_BYTE, TSP_GA_MASTER_PROCESS, MPI_COMM_WORLD);

    for (int it = 0; it < iter_total / iter_one; it++) {
        for (int i = 0; i < iter_one; i++) {
            population_update(&population[world_rank * population_size], &buffer[world_rank * population_size], rng);
        }
        printf("Process %d iter %d best path length %f\n", world_rank, (it + 1) * iter_one,
            population[world_rank * population_size].path_length);
        fflush(stdout);
        if (world_rank == TSP_GA_MASTER_PROCESS)
        {
            for (int i = 1; i < world_size; i++)
            {
                MPI_Recv(&population[i * population_size], population_size * sizeof(chromosome), MPI_BYTE, i, 0, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
                // printf("Process MASTER received population from process %d\n", i);
            }
            std::shuffle(population, &population[world_size * population_size], rng);
        }
        else
        {
            MPI_Send(&population[world_rank * population_size], population_size * sizeof(chromosome), MPI_BYTE, TSP_GA_MASTER_PROCESS, 0,
                MPI_COMM_WORLD);
            // printf("Process %d sent population to MASTER\n", world_rank);
        }
        MPI_Bcast(population, world_size * population_size * sizeof(chromosome), MPI_BYTE, TSP_GA_MASTER_PROCESS, MPI_COMM_WORLD);
        // printf("Process MASTER broadcast population\n");
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}