#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <random>
#include <semaphore>
#include <thread>

#define city_total 280
#define testcase_name "a280.txt"
#define population_size 2000
#define iter_total 10000
#define iter_one 200
double city_location[city_total][2];
double distance_matrix[city_total][city_total];

bool caculation_end = false;
std::counting_semaphore sem_calc_thread(0);
std::counting_semaphore sem_main_thread(0);

#define city_distance(idx, idy) distance_matrix[(idx)][(idy)]

struct chromosome {
  int path[city_total];
  double path_length;
};
struct population {
  chromosome individuals[population_size];
  chromosome buffer[population_size];
};

int rand_between(int start, int end, std::mt19937 &rng) {
  return rng() % (end - start) + start;
}
double genrand_real(std::mt19937 &rng) { return rng() * (1.0 / 4294967296.0); }
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

void path_random_init(int *path, std::mt19937 &rng) {
  for (int i = 0; i < city_total; i++) {
    path[i] = i;
  }
  std::shuffle(path, path + city_total, rng);
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
bool chromosome_cmp_func(const chromosome &a, const chromosome &b) {
  return a.path_length < b.path_length;
}

population *population_create() {
  return (population *)calloc(sizeof(struct population), 1);
}

void population_destory(population *p) { free(p); }

void population_init(population *p, std::mt19937 &rng) {
  for (int i = 0; i < population_size; i++) {
    path_random_init(p->individuals[i].path, rng);
    p->individuals[i].path_length = path_length(p->individuals[i].path);
  }
}

void population_sort(population *p) {
  std::sort(p->individuals, p->individuals + population_size,
            chromosome_cmp_func);
}

void population_select(population *p, int elite_size, std::mt19937 &rng) {
  int random_city_idx[population_size];
  for (int i = 0; i < population_size; i++) {
    random_city_idx[i] = i;
  }
  std::shuffle(random_city_idx, random_city_idx + population_size, rng);
  for (int i = elite_size; i < population_size; i++) {
    p->buffer[i] = p->individuals[random_city_idx[i]];
  }
  memcpy(&p->individuals[elite_size], &p->buffer[elite_size],
         sizeof(struct chromosome) * (population_size - elite_size));
}

void population_breed(struct population *p, int elite_size, std::mt19937 &rng) {
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

void population_mutate(struct population *p, double rate, std::mt19937 &rng) {
  for (int i = 2; i < population_size; i++) {
    if (genrand_real(rng) > rate) {
      continue;
    }
    struct chromosome *pp = &p->individuals[i];
    int x = rand_between(0, city_total, rng);
    int y = rand_between(0, city_total, rng);
    int t = pp->path[x];
    pp->path[x] = pp->path[y];
    pp->path[y] = t;
    pp->path_length = path_length(pp->path);
  }
}

double population_update(struct population *p, std::mt19937 &rng) {
  population_sort(p);
  double best_path_length = p->individuals[0].path_length;
  population_select(p, population_size / 5, rng);
  population_breed(p, population_size / 5, rng);
  population_mutate(p, 0.1, rng);
  return best_path_length;
}

void thread_main(population *p) {
  std::mt19937 thread_local rng(std::random_device{}());
  double curr_best;
  population_init(p, rng);
  while (!caculation_end) {
    sem_calc_thread.acquire();
    for (int i = 0; i < iter_one; i++) {
      curr_best = population_update(p, rng);
    }
    sem_main_thread.release();
  }
}

int main() {
  citys_init();
  std::mt19937 rng(std::random_device{}());
  std::thread v[4];
  population *p[4];
  chromosome *main_chromosome_buffer = new chromosome[4 * population_size];
  for (int i = 0; i < 4; i++) {
    p[i] = population_create();
    v[i] = std::thread(thread_main, p[i]);
  }
  for (int i = 0; i < iter_total / iter_one; i++) {
    sem_calc_thread.release(4);
    for (int i = 0; i < 4; i++) {
      sem_main_thread.acquire();
    }
    for (int i = 0; i < 4; i++) {
      printf("thread %i best path length %f\n", i,
             p[i]->individuals->path_length);
      memcpy(&main_chromosome_buffer[i * population_size], p[i]->individuals,
             sizeof(chromosome) * population_size);
    }
    std::shuffle(main_chromosome_buffer,
                 &main_chromosome_buffer[4 * population_size], rng);
    for (int i = 0; i < 4; i++) {
      memcpy(p[i]->individuals, &main_chromosome_buffer[i * population_size],
             sizeof(chromosome) * population_size);
    }
  }
  caculation_end = true;
  sem_calc_thread.release(4);
  for (int i = 0; i < 4; i++) {
    v[i].join();
    population_destory(p[i]);
  }
}
