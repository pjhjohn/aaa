#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits>
#include<fstream>

/* PROJECT 2 ESSENTIALS */
#define MAXN    600 // Maximum value of N ( of project 2 )
#define PSIZE   10  // Size of the population

double X[MAXN], Y[MAXN];    // (X[i], Y[i]) := the i-th location
int N;                      // The number of locations
long long TimeLimit;        // Time Limitation
double Dist[MAXN][MAXN];    // Distance Map

typedef struct {            // Note that the representation is currently order-based.
    int ch[MAXN];           // chromosome : a permutation of (0, 1, ..., N-1).
    double f;               // fitness
} SOL;
SOL population[PSIZE];      // Population
SOL record;                 // Best solution ( eval() automatically updates this )

double eval(SOL *s) {       // calculate the fitness of s and store it into s->f
    s->f = 0;
    for (int i = 0; i < N; i++) s->f += Dist[s->ch[i]][s->ch[(i+1)%N]];
    if (record.f > s->f) {
        record.f = s->f;
        for (int i = 0; i < N; i++) record.ch[i] = s->ch[i];
    } return s->f;
}
double fitness(SOL *s) {
    double f = 0;
    for (int i = 0; i < N; i++) f += Dist[s->ch[i]][s->ch[(i+1)%N]];
    return f;
}

void gen_rand_sol(SOL *s) { // generate a random order-based solution at s
    for (int i = 0; i < N; i++) s->ch[i] = i;
    for (int i = 0; i < N; i++) {
        int r = i + rand() % (N - i);   // r is a random number in [i..N-i)
        int tmp = s->ch[i]; s->ch[i] = s->ch[r]; s->ch[r] = tmp; // swap
    } eval(s);
}


void pprint(const char* TAG, const SOL *s) {
  printf("<%s>\n", TAG);
  for(int i = 0; i < N; i ++) {
      printf("%3d ", s->ch[i]);
      if(i%10==9) printf("\n");
  } printf("%lf\n", s->f);
}



/* Selection Algorithms */
/* Roulette - Wheel Selection */
void RouletteWheelSelection(SOL **p) {
  int i;
  // find min, max, total fitness
  double f_min = std::numeric_limits<double>::max();
  double f_max = std::numeric_limits<double>::min();
  for(i = 0; i < PSIZE; i ++) {
    double f_current = population[i].f;
    if(f_min > f_current) f_min = f_current;
    if(f_max < f_current) f_max = f_current;
  }
  /* Start shooting Roulette */
  double K = 3;
  double f_total = 0;
  for(i = 0; i < PSIZE; i ++) f_total += (f_min - population[i].f) + (f_min - f_max)/(K-1);
  double point = static_cast <float> (rand()) / static_cast <float> (f_total);
  double sum = 0;
  for(i = 1; i < PSIZE; i++) {
    double f_current = population[i].f;
    double region = (f_min - f_current) + (f_min - f_max)/(K-1);
    sum += region;
    if(point < sum) *p = &population[i];
  }
}

/* Random Selection */
void randomSelection(SOL **p) {
    *p = &population[rand() % PSIZE];
}

// choose one solution from the population
void selection(SOL **p) {
  randomSelection(p);
}

/* CROSSOVER */
int s1, s2; // slice index
void setRandomSliceIndexes() {
  s1 = rand()%N;
  s2 = s1;
  while(s1==s2) {
    s2 = rand()%N;
  }
  if (s1 > s2) {
    int temp = s2;
    s2 = s1;
    s1 = temp;
  }
}
int hasDuplication(const int index, const int *ch, int from, int to) {
  int i;
  for(i = from; i <= to; i ++) {
    if(index == ch[i]) return 1;
  } return 0;
}
int isOkayToAppend(const int index, const int *ch) {
  return hasDuplication(index, ch, 0, N-1) == 0 ? 1 : 0;
}
int isOffspringCompleted(const SOL *offspring) {
  int sum = 0;
  for(int i = 0; i < N; i ++) {
      if(offspring->ch[i] == -1) return 0;
      else sum += offspring->ch[i];
  }
  return sum == N * (N-1) / 2;
}
void PMXCrossover(const SOL *p1, const SOL *p2, SOL *c, const int stackTrace) {
  setRandomSliceIndexes();
  /* Copy Parent1 in Slice btw s1 and s2 */
  if(stackTrace) {pprint("p1", p1); pprint("p2", p2); printf("slice from %d to %d\n", s1, s2);}
  int i, remains = N - (s2 - s1 + 1);
  for(i =  0; i <  N; i ++) c->ch[i] = -1;
  if(stackTrace) pprint("after init", c);
  for(i = s1; i <= s2; i ++) c->ch[i] = p1->ch[i];
  if(stackTrace) pprint("after copy", c);
  int current = (s2+1)%N;
  int j = s1;
  while (remains > 0) {
    if(isOkayToAppend(p2->ch[current], c->ch)) {
      c->ch[current] = p2->ch[current];
      current = (current+1)%N;
      remains --;
      if(stackTrace) pprint("after safe", c);
    } else {
      int appended = 0;
      while(!appended) {
        if(hasDuplication(p2->ch[j], p1->ch, s1, s2)) {
            if(stackTrace) printf("collision. passing @ %d\n", j);
            j = (j+1)%N;
        }
        else {
          c->ch[current] = p2->ch[j];
          current = (current+1)%N;
          appended = 1;
          remains --;
          if(stackTrace) { printf("takes # @ %d\n", j); pprint("after notsafe", c); }
          j = (j+1)%N;
        }
      }
    }
  }
  if (isOffspringCompleted(c)) {
    return;
  } else {
    if(stackTrace) exit(1);
    else {
        PMXCrossover(p1, p2, c, 1);
        exit(1);
    }
  }
}
int duplicated(const int val, const int *ch, int from, int to) {
    int i;
    for(i = from; i < to; i ++) if (val == ch[i]) return 1;
    return 0;
}
void OrderCrossover(const SOL *p1, const SOL *p2, SOL *c) {
    setRandomSliceIndexes();
    int i;
    for(i = 0; i < N; i++) c->ch[i] = -1;
    for(i = s1; i < s2; i++) c->ch[i] = p1->ch[i];

    int iC = 0, iP2 = 0;
    while (iP2 < N) {
        if(duplicated(p2->ch[iP2], p1->ch, s1, s2)) {
            iP2 ++;
        } else {
            if(iC == s1) {
                iC += s2 - s1;
            } else {
                c->ch[iC] = p2->ch[iP2];
                iC ++;
                iP2 ++;
            }
        }
    }
}
void directFromFirst(const SOL *p1, const SOL *p2, SOL *c) {
    for (int i = 0; i < N; i++) {
        c->ch[i] = p1->ch[i];
    }
}
void crossover(const SOL *p1, const SOL *p2, SOL *c) {
  OrderCrossover(p1, p2, c);
  eval(c);
}

/* MUTATION */
void inversionMutation(SOL *s) {
    setRandomSliceIndexes();
    int i, temp;
    for(i = 0; i < (s2 - s1) / 2; i ++) {
        temp = s->ch[s1+i];
        s->ch[s1+i] = s->ch[s2-i];
        s->ch[s2-i] = temp;
    }
}
void mutation(SOL *s) {
    inversionMutation(s);
    eval(s);
}

/* LOCAL OPTIMIZATION */
void inverse(SOL *s, int from, int to) {
    int i, temp;
    for(i = 0; i < (to - from + 1) / 2; i ++) {
        temp = s->ch[from+i];
        s->ch[from+i] = s->ch[to-i];
        s->ch[to-i] = temp;
    }
}
void copy(SOL *target, SOL *result) {
    result->f = target->f;
    for(int i = 0; i < N; i ++) result->ch[i] = target->ch[i];
}
void _2_opt_(SOL *s) {
    int from, to;
    SOL sol, inversed;
    copy(s, &sol);
    while(true) {
        double solf = fitness(&sol);
        int better = false;
        for(from = 0; from < N - 1; from ++) {
            for(to = from + 1; to < N; to ++) {
                copy(&sol, &inversed);
                inverse(&inversed, from, to);
                double invf = fitness(&inversed);
                if (solf > invf) {
                    copy(&inversed, &sol);
                    better = true;
                    break;
                } else continue;
            } if(better) break;
        }
        if (!better) break;
    }
    copy(&sol, s);
}
void local_optimization(SOL *s) {
    _2_opt_(s);
    eval(s);
}

/* REPLACEMENT */
int getIndexOfWorst() {
  double f_max = std::numeric_limits<double>::min();
  int index = -1;
  for(int i = 0; i < PSIZE; i ++) {
    double f_current = population[i].f;
    if(f_max < f_current){
      f_max = f_current;
      index = i;
    }
  } return index;
}
void replacement(const SOL *offspr) {
    int i, p = getIndexOfWorst();
    population[p].f = offspr->f;
    for (i = 0; i < N; i++) {
        population[p].ch[i] = offspr->ch[i];
    }
}

int GENERATION;
// a "steady-state" GA
void GA() {
    int i;
    SOL c;
    time_t begin = time(NULL);

    for (i = 0; i < PSIZE; i++) {
        gen_rand_sol(&population[i]);
    }

    while (1) {
        if(time(NULL) - begin >= TimeLimit - 1) return; // end condition
        SOL *p1, *p2;
        selection(&p1); selection(&p2);
        crossover(p1, p2, &c);
        mutation(&c);
        local_optimization(&c);
        replacement(&c);
        GENERATION ++;
    }
}


// read the test case from stdin
// and initialize some values such as record.f and Dist
void init() {
    FILE *pf = fopen("../input/cycle.in.200", "r");
    int i, j, tmp;
    double time_limit;

    tmp = fscanf(pf, "%d", &N);
    for (i = 0; i < N; i++) {
        tmp = fscanf(pf, "%lf %lf", &X[i], &Y[i]);
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            double dx = X[i] - X[j], dy = Y[i] - Y[j];
            Dist[i][j] = sqrt(dx*dx + dy*dy);
        }
    }
    tmp = fscanf(pf, "%lf", &time_limit);
    TimeLimit = (long long) time_limit;

    record.f = 1e100;
    GENERATION = 0;
}


// print the best solution found to stdout
void answer() {
    int i;

    for (i = 0; i < N; i++) {
        if (i > 0) printf(" ");
        printf("%d", record.ch[i]+1);
    }
    printf(" %lf", record.f);
    printf(" GENERATION : %d", GENERATION);
    printf("\n");
}

int main() {
    srand(time(NULL));
  int nLoop = 1;
  for(int loop = 0; loop < nLoop; loop++) {
    init();
    GA();
    answer();
  }
    return 0;
}