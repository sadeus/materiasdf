
/*
    Simulación de un sistema 2D Ising por medio de un algoritmo
    Metropolis Monte Carlo NVT (canónico)

    Para compilar
    >>> gcc -o3 ising ising.c -lm


    Parametros
    ----------------
    -L (int) tamaño de la red LxL
    -T (double) inverso de la temperatura adimensional (beta = J/(k*T))
    -n (int) cantidad de muestras a tomar
    -fs (int) cantidad de pasos entre muestra y muestra
    -nT (int) cantidad de pasos hasta termalizar
    -o (string) archivo donde agrega los datos finalmente

    Resultados
    ---------------
    Imprime en pantalla la temperatura, magnetización y la energía,
    además lo guarda en el archivo indicado (por defecto en ./med)

    Transfondo
    --------------
    El Hamiltoneano del modelo de Ising es

    H = -J sum_<ij> s_i * s_j,

    donde "sum_<ij>" es la suma sobre los primeros vecinos de la
    posición (i,j). Esta suma depende de la forma de la red, en el
    caso 2D se consideran los primeros vecinos a los spines arriba,
    abajo y a los costados. La proyección de spin  utilizada
    corresponde a s_i = +/- 1, y J es un parámetro arbitrario de
    acomplamiento.

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//Parametros del programa
struct Parametros {
    int L;
    double beta;
    int n_samp;
    int n_term;
    int f_samp;
    int seed;
};


//Computa el cambio de energía debido al giro del spin i,j,
//que corresponde a la energía asociada a los primeros vecinos
//De esa forma DeltaE puede valor -8,-4,0,4,8 (sin unidades)
int DeltaE(char *F, int L, int i) {
	//Cantidad total de elementos. Matriz de N,N.
	//Para portar a una red genérica se podría usar una clase
  int N = L * L;
  if (i < 0 || i > N) {
    printf("Error en los índices");
    return;
	}
	//Condiciones de contorno periodica en ambas direcciones
  int i_izq =  i % L == 0 ? i + (L - 1) : i - 1;
  int i_der = (i + 1) % L == 0 ? i + 1 - L : i + 1;
  int i_up = i - L < 0 ? i + N - L : i - L;
  int i_down = i + L + 1 > N ? i + L - N : i + L;
  int de = -2 * *(F + i) * (*(F + i_izq) + *(F + i_der) + *( F + i_down) + *(F + i_up));
  return de;

}

//Función que mide la magnetización y la energía por sitio
//e y m se pasan por referencia, energía y magnetización
void sampleo(char * F, struct Parametros param, long * m, long * e){
  int i, i_der, i_down;
  int N = param.L*param.L;
  for (i = 0; i < N; i++) {
      *m += *(F + i);
      i_der = (i + 1) % param.L == 0 ? i + 1 - param.L : i + 1;
      i_down = i + param.L + 1 > N ? i + param.L - N : i + param.L;
      *e -= *(F + i) * (*(F + i_der) + *(F + i_down));
  }
}


//Genera el array de forma aleatoria
//Usa un entero como un bitfield, si es 1 => 1, 0 => -1
void initRandom(char *F, int L) {
  int i;
  int N = L*L;
  while (N) {
      int r = rand();
      for (i = 0; N && i < sizeof(int); ++i) {
          F[--N] = r >> i & 1 ? -1 : 1;
      }
  }
}

//Genera un array de forma ordenada, todos 1
void initOrd(char *F, int L) {
  int i;
  int N = L * L;
  for (i = 0; i < N; i++) {
      if (i % 2)
          F[i] = 1;
      else
          F[i] = -1;
  }
}


void parseParametros(int argc, char * argv[], struct Parametros *param) {
  int i = 0;
  for (i = 0; i < argc; i++)
  {
    if (!strcmp(argv[i],"-L"))
        param->L = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T"))
        param->beta = atof(argv[++i]);
    else if (!strcmp(argv[i],"-n"))
        param->n_samp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nT"))
        param->n_term = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-fs"))
        param->f_samp = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-s"))
        param->seed = (unsigned long)atoi(argv[++i]);
  }
}

void mmc(char *F, struct Parametros param) {
    int de = 0; long E = 0; double e_s?um = 0.0; //Energía
    long M = 0; double m_sum = 0.0; //Magentización
    int N = param.L * param.L;
    double p[2]; //Factor de aceptación
    double r; //Número aleatorio, para comprar con la aceptación

    int n_samp = 1;  //Cantidad de mediciones
    int n_iter = param.f_samp * param.n_samp + param.n_term;
    int i, j, k;?
    sampleo(F, param, &M, &E);
    printf("%f %f\n", (float)M/N, (float)E/N);
    p[0] = exp(-4.0 * param.beta); //Prealoco las exponenciales.
    p[1] = exp(-8.0 * param.beta);
    for (j = 0; j < n_iter; j++) {
        //Por cada paso de MC hace N giros de spin.
        for (k = 0; k < N; k++) {
          i = rand() % N; //posición al azar del array
          de = DeltaE(F, param.L, i); //Delta de Enegía al posiblemente cambiar el estado.
          if (de < 0) {
              *(F + i) *= -1;
              M += *(F + i);
              E += de;
          }
          if (de != 0) {
              r = (double)rand()/RAND_MAX * 1.0; //Número aleatorio en [0,1].
              if (p[de == 4 ? 0 : 1] > r){
                  printf("Entró en transiciónn %d %d\n", j, k);
                  *(F + i) *= -1; //Acepta el cambio => Invierte el spin.
                  M += *(F + i);
                  E += de;
              }
          }
        }
        //Samplea n_samp veces, después de que haya termalizado (n_term pasos).
        if ((j % param.f_samp) == 0 && j > param.n_term) {
            //sampleo(F, param, &m, &e);
            m_sum += (double)M/N;
            e_sum += (double)E/N;
            n_samp++;
        }
    }
    m_sum = fabs(m_sum) / n_samp;
    e_sum /= n_samp;
    printf("%f %f %f\n", param.beta, m_sum, e_sum);
}


int main(int argc, char * argv[]) {
    struct Parametros param; //Parametros
    param.L = 16;
    param.n_samp = 200;
    param.f_samp = 100;
    param.n_term = 500;
    param.seed = time(NULL);
    int N = param.L * param.L;
    //Inicialización
    parseParametros(argc, argv, &param); //Parsea los parámetros.
    N = param.L * param.L;
    srand(param.seed);
    char F[N];
    initRandom(F, param.L);
    //Código
    mmc(F, param);
    return 0;
}
