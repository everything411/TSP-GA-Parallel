#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define city_total 280
#define testcase_name "a280.txt"
#define population_size 2000
#define iter_total 10000
double city_location[city_total][2];
double distance_matrix[city_total][city_total];
#define city_distance(idx, idy) distance_matrix[(idx)][(idy)]

struct chromosome {
  int path[city_total];
  double path_length;
};
struct population {
  struct chromosome individuals[population_size];
  struct chromosome buffer[population_size];
};

int rand_between(int start, int end) { return rand() % (end - start) + start; }
double genrand_real(void) { return rand() * (1.0 / 4294967296.0); }
void citys_init() {
  FILE *fp = fopen(testcase_name, "r");
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

void path_random_init(int *path) {
  for (int i = 0; i < city_total; i++) {
    path[i] = i;
  }
  for (int i = 0; i < city_total; i++) {
    int p = rand() % city_total;
    int t = path[i];
    path[i] = path[p];
    path[p] = t;
  }
}

void path_duplicate(int *to, const int *from) {
  memcpy(to, from, sizeof(int) * city_total);
}
double path_length(int *path) {
  double length = city_distance(path[0], path[city_total - 1]);
  for (int i = 0; i < city_total - 1; i++) {
    length += city_distance(path[i], path[i + 1]);
  }
  return length;
}

int sgn(double val) { return (0 < val) - (val < 0); }

// compare by fitness (path length)
int chromosome_cmp_func(const void *a, const void *b) {
  double fit = ((const struct chromosome *)a)->path_length -
               ((const struct chromosome *)b)->path_length;
  return sgn(fit);
}

struct population *population_create() {
  return calloc(sizeof(struct population), 1);
}

void population_destory(struct population *p) { free(p); }

void population_init(struct population *p) {
  for (int i = 0; i < population_size; i++) {
    path_random_init(p->individuals[i].path);
    p->individuals[i].path_length = path_length(p->individuals[i].path);
  }
}

void population_sort(struct population *p) {
  qsort(p->individuals, population_size, sizeof(struct chromosome),
        chromosome_cmp_func);
}

void population_select(struct population *p, int elite_size) {
  int random_city_idx[population_size];
  for (int i = 0; i < population_size; i++) {
    random_city_idx[i] = i;
  }
  for (int i = 0; i < population_size; i++) {
    int p = rand() % population_size;
    int t = random_city_idx[i];
    random_city_idx[i] = random_city_idx[p];
    random_city_idx[p] = t;
  }
  for (int i = elite_size; i < population_size; i++) {
    p->buffer[i] = p->individuals[random_city_idx[i]];
  }
  memcpy(&p->individuals[elite_size], &p->buffer[elite_size],
         sizeof(struct chromosome) * (population_size - elite_size));
}

void population_breed(struct population *p, int elite_size) {
  int cross_path[city_total];
  int diff_path[city_total];
  int new_path1[city_total], new_path2[city_total];
  int path_in_cross_path_flag[city_total];

  for (int i = elite_size; i < population_size; i++) {
    int idx, idy;
    do {
      idx = rand_between(0, population_size);
      idy = rand_between(0, population_size);
    } while (idx == idy);
    memset(path_in_cross_path_flag, 0, sizeof(path_in_cross_path_flag));
    int x, y;
    int cross_length;
    do {
      x = rand_between(0, city_total);
      y = rand_between(0, city_total);
    } while (x >= y);
    cross_length = y - x;
    memcpy(cross_path, p->individuals[idx].path, sizeof(int) * cross_length);
    for (int i = 0; i < cross_length; i++) {
      path_in_cross_path_flag[cross_path[i]] = 1;
    }
    int diff_length = 0;
    for (int i = 0; i < city_total; i++) {
      if (!path_in_cross_path_flag[p->individuals[idy].path[i]]) {
        diff_path[diff_length++] = p->individuals[idy].path[i];
      }
    }
    memcpy(new_path1, cross_path, sizeof(int) * cross_length);
    memcpy(new_path1 + cross_length, diff_path, sizeof(int) * diff_length);
    memcpy(new_path2, diff_path, sizeof(int) * diff_length);
    memcpy(new_path2 + diff_length, cross_path, sizeof(int) * cross_length);
    double pl1 = path_length(new_path1);
    double pl2 = path_length(new_path2);
    int *np;
    double nl;
    if (pl1 < pl2) {
      np = new_path1;
      nl = pl1;
    } else {
      np = new_path2;
      nl = pl2;
    }
    path_duplicate(p->buffer[i].path, np);
    p->buffer[i].path_length = nl;
  }
  memcpy(&p->individuals[elite_size], &p->buffer[elite_size],
         sizeof(struct chromosome) * (population_size - elite_size));
}

void population_mutate(struct population *p, double rate) {
  for (int i = 0; i < population_size; i++) {
    if (genrand_real() > rate) {
      continue;
    }
    struct chromosome *pp = &p->individuals[i];
    int x = rand_between(0, city_total);
    int y = rand_between(0, city_total);
    int t = pp->path[x];
    pp->path[x] = pp->path[y];
    pp->path[y] = t;
    pp->path_length = path_length(pp->path);
  }
}

double population_update(struct population *p) {
  population_sort(p);
  double best_path_length = p->individuals[0].path_length;
  population_select(p, population_size / 5);
  population_breed(p, population_size / 5);
  population_mutate(p, 0.1);
  return best_path_length;
}

int main() {
  srand(time(NULL));
  citys_init();
  struct population *population = population_create();
  double curr_best;
  population_init(population);
  for (int i = 0; i < iter_total; i++) {
    curr_best = population_update(population);
    if (i % 100 == 0) {
      printf("best path length %f\n", curr_best);
    }
  }
  population_destory(population);
}
