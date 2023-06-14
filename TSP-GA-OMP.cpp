#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <random>
#include <omp.h>

#define thread_count 4
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
chromosome population[population_size * thread_count];
chromosome buffer[population_size * thread_count];

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


int main() {
    omp_set_num_threads(thread_count);
    citys_init();
    std::mt19937 rng[thread_count];
    for (int i = 0; i < thread_count; i++)
    {
        rng[i] = std::mt19937(std::random_device{}());
        population_init(&population[i * population_size], rng[i]);
    }

    for (int it = 0; it < iter_total / iter_one; it++) {
#pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            for (int i = 0; i < iter_one; i++) {
                population_update(&population[thread_id * population_size], &buffer[thread_id * population_size], rng[thread_id]);
            }
            printf("iter %d thread %d best path length %f\n", (it + 1) * iter_one, thread_id,
                population[thread_id * population_size].path_length);
        }
        std::shuffle(population, &population[thread_count * population_size], rng[0]);
    }
}
