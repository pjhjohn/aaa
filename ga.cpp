#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits>

#define MAXN	318 // Maximum value of N
#define PSIZE	10 // Size of the population

/*****************************************************************
	Input variables
*****************************************************************/
// (X[i], Y[i]) := the i-th location
double X[MAXN], Y[MAXN];
// The number of locations
int N;

// Time limit for the test case
long long TimeLimit;


// Dist[i][j] := the distance between (X[i], Y[i]) and (X[j], Y[j])
// will be automatically calculated
double Dist[MAXN][MAXN];



/*****************************************************************
	GA variables and functions
	Note that the representation is currently order-based.
	A chromosome will be a permutation of (0, 1, ..., N-1).
*****************************************************************/
typedef struct {
	int ch[MAXN];	// chromosome
	double f;		// fitness
} SOL;

// population of solutions
SOL population[PSIZE];

// best (found) solution
// eval() updates this
SOL record;

// calculate the fitness of s and store it into s->f
double eval(SOL *s) {
	int i;

	s->f = 0;
	for (i = 0; i < N; i++) {
		s->f += Dist[s->ch[i]][s->ch[(i+1)%N]];
	}

	if (record.f > s->f) {
		record.f = s->f;
		for (i = 0; i < N; i++) {
			record.ch[i] = s->ch[i];
		}
	}

	return s->f;
}


// generate a random order-based solution at s
void gen_rand_sol(SOL *s) {
	int i;
	for (i = 0; i < N; i++) {
		s->ch[i] = i;
	}
	for (i = 0; i < N; i++) {
		int r = i + rand() % (N - i);	// r is a random number in [i..N-i)
		int tmp = s->ch[i]; s->ch[i] = s->ch[r]; s->ch[r] = tmp; // swap
	}
	
	// calculate the fitness
	eval(s);
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
void PMXCrossover(const SOL *p1, const SOL *p2, SOL *c) {
  setRandomSliceIndexes();
  /* Copy Parent1 in Slice btw s1 and s2 */
  int i, remains = N - (s2 - s1 + 1);
  for(i =  0; i <  N; i ++) c->ch[i] = -1;
  for(i = s1; i <= s2; i ++) c->ch[i] = p1->ch[i];
  int current = (s2+1)%N;
  int j = s1;
  while (remains > 0) {
    if(isOkayToAppend(p2->ch[current], c->ch)) {
      c->ch[current] = p2->ch[current];
      current = (current+1)%N;
      remains --;
    } else {
      int appended = 0;
      while(!appended) {
        if(hasDuplication(p2->ch[j], p1->ch, s1, s2)) {
			j = (j+1)%N;
		}
        else {
          c->ch[current] = p2->ch[j];
          current = (current+1)%N;
          appended = 1;
          remains --;
          j = (j+1)%N;
        }
      }
    }
  }
  if (isOffspringCompleted(c)) {
    return;
  } 
}

void directFromFirst(const SOL *p1, const SOL *p2, SOL *c) {
	for (int i = 0; i < N; i++) {
		c->ch[i] = p1->ch[i];
	}
}
void crossover(const SOL *p1, const SOL *p2, SOL *c) {
  PMXCrossover(p1, p2, c);
  eval(c);
}


// mutate the solution s
// currently this operator does nothing
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


// replace one solution from the population with the new offspring
// currently any random solution can be replaced
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
double getAverageFitness() {
	double total = 0.0f;
	for(int i = 0; i < PSIZE; i++) {
		total += population[i].f;
	} return total / PSIZE;
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

	FILE *pf = fopen("output", "w");
	while (1) {
		if(time(NULL) - begin >= TimeLimit - 1) return; // end condition
		SOL *p1, *p2;
		selection(&p1); selection(&p2);
		crossover(p1, p2, &c);
		mutation(&c);
		replacement(&c);
		GENERATION ++;
	} fclose(pf);
}


// read the test case from stdin
// and initialize some values such as record.f and Dist
void init() {
	int i, j, tmp;
	double time_limit;

	tmp = scanf("%d", &N);
	for (i = 0; i < N; i++) {
		tmp = scanf("%lf %lf", &X[i], &Y[i]);
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			double dx = X[i] - X[j], dy = Y[i] - Y[j];
			Dist[i][j] = sqrt(dx*dx + dy*dy);
		}
	}
	tmp = scanf("%lf", &time_limit);
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
	printf("\n");
}


int main() {
	srand(time(NULL));
  	init();
  	GA();
  	answer();
	return 0;
}


