// nombre de commandes
int n = ...;
range N = 0..n;
range O = 1..n;

// nombre de periodes de tarrifs TOU
int m = ...;


// nombre de periodes d'emission de CO2 par jour
int h = ...;

/* parametres */
// due dates
int d[0..n] = ...;
// deadlines
int db[0..n] = ...;
// profits
int e[0..n] = ...;
// processing time
int p[0..n] = ...;
// arrival time
int r[0..n] = ...;

//penalite de retard sur la commande i
float w[0..n] = ...;

// cout de setup entre les commandes i et j 
int s[0..n][0..n] = ...;

// debut de la periode k d'emission CO2 
int g[0..h-1] = ...;

// emission CO2/kWh durant la periode k
float q[1..h] = ...;

// taxe sur l'emission CO2 /kg
float Tax = ...;

// debut de la periode k de TOU
int b[0..m-1] = ...;

// cout de l'electricite sur la periode k
float EC[1..m] = ...;

// puissance (kW) de consommation de la commande i
int Omega[1..n] = ...;


/* time indexed */
// horizon
int T = max(i in O) (db[i]) +1;
int TT = max(i in O)(db[i]);
int pmax = max(i in 0..n) p[i];

// f[j][t] = profit généré a l'instant de debut
float f[j in 1..n][t in 0..T] = (t>=r[j])&&(t<=db[j])?(e[j] - w[j]*maxl(t-d[j],0) ):0;


// cout TOU & emission CO2 a l'instant t pour la commande j
/*float c[j in 0..n][t in 0..T] = ((t > T) ? 0 : ((j == 0) ? 0 : (1/60)*Omega[j]*sum(k in 1..m-1)(((t>=b[k-1])&&(t<b[k]))*EC[k])
  +(1/60)*Tax*Omega[j]*sum(l in 1..h-1) (((t>=g[l-1])&&(t<g[l]))*q[l])));*/
 float c[j in O][t in 0..T] = (t>=r[j]) && (t<=db[j])?
(1/60)*Omega[j]*sum(k in 1..m-1)(((t>=b[k-1])&&(t<b[k]))*EC[k])
  +(1/60)*Tax*Omega[j]*sum(l in 1..h-1) (((t>=g[l-1])&&(t<g[l]))*q[l]):0;

  // convenient
int pp[i in 0..n][j in 0..n] = (( i == j && i != 0)?1:p[i]);

int smin = min(i in 0..n, j in 1..n: i!=j) s[i][j];
int smax = max(i in 0..n, j in 1..n: i!=j) s[i][j];
int Smax[j in 1..n] = max(i in 0..n) s[i][j];
int Smin[j in 1..n] = min(i in 0..n: i!=j) s[i][j];

float C[i in 0..n][j in 1..n][t in 0..T] = ((t>=r[j]+s[i][j]+1) && (t<=db[j]-p[j]+1) )? sum(tp in t-s[i][j]..t-1) c[j][tp] : 0;

