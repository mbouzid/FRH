
// nombre de commandes
int n = ...;
range N = 0..n;
range O = 1..n;


// nombre de périodes de tarrifs TOU
int m = ...;


// nombre de périodes d'emission de CO2 par jour
int h = ...;

/* parametres */
// due dates
int d[N] = ...;
// deadlines
int db[N] = ...;
// profits
int e[N] = ...;
// processing time
int p[N] = ...;
// arrival time
int r[N] = ...;

//pénalité de retard sur la commande i
float w[N] = ...;

// cout de setup entre les commandes i et j 
int s[N][N] = ...;

// début de la période k d'émission CO2 
int g[0..h-1] = ...;

// emission CO2/kWh durant la période k
float q[1..h] = ...;

// taxe sur l'emission CO2 /kg
float Tax = ...;

// début de la période k de TOU
int b[0..m-1] = ...;

// cout de l'électricité sur la période k
float EC[1..m] = ...;

// puissance (kW) de consommation de la commande i
int Omega[O] = ...;


/* time indexed */
// horizon
int T = max(i in O) (db[i]+1) ;


int SMax[j in 1..n] = max(i in 0..n:i!=j) s[i][j];

// f[j][t] = profit généré a l'instant de debut
float f[j in 1..n][t in 0..T] = (t>=r[j]+p[j])&&(t<=db[j])?(e[j] - w[j]*maxl(t-d[j],0) ):0;


// c[j][t] = cout energetique a la periode t
float c[j in O][t in 0..T] = (t>=r[j]) && (t<=db[j]) ?
(1/60)*Omega[j]*sum(k in 1..m)(((t>=b[k-1])&&(t<b[k]))*EC[k])
  +(1/60)*Tax*Omega[j]*sum(l in 1..h) (((t>=g[l-1])&&(t<g[l]))*q[l]):0;



/*variable de decisions*/
// x[j][t]= 1 ssi commande j commence à l'instant t
dvar boolean x[j in 0..n][t in 0..T];
// u[i][j] = 1 ssi commande i precede directement commande j
dvar boolean u[i in 0..n][j in 0..n]; 

dvar boolean isProcessingApprox[i in 1..n];

int smin = min(i in 1..n, j in 1..n: i!=j) s[i][j];
int smax = max(i in 1..n, j in 1..n: i!=j) s[i][j];


maximize sum(j in 1..n) 
sum(t in r[j]..db[j]-p[j]+1) 
( 

	x[j][t]*(
	
	f[j][t+p[j]-1]
	- (
		
		sum(tp in 0..p[j]-1) c[j][t+tp] 
		+
		sum(i in 0..n: i!=j)
			u[i][j]* 
			sum(tp in 1..s[i][j]: t-tp>=0) c[j][t-tp]
		))
)
;

subject to 
{

  	succ: forall(i in 0..n)
  	{
  		sum(j in 1..n: j!=i) u[i][j] <= sum(t in r[i]..db[i]-p[i]+1) x[i][t];  	
  	}
  	  
  	pred: forall(i in 1..n)
  	{ 
  		sum(j in 0..n:j!=i) u[j][i] == sum(t in r[i]..db[i]-p[i]+1) x[i][t];
  	}	

	capacity: forall(t in 1..T)
	{
		sum(j in 1..n) (x[j][t]) <= 1;
	}
	
	capacity1: forall(j in 0..n)
	{
		sum(t in r[j]..db[j]-p[j]+1)	 x[j][t] <= 1;
	}
	

	prec:forall(j in 1..n, i in 0..n: i!=j)
	{
	 	sum(t in r[i]..db[i]-p[i]+1) t*x[i][t] + (p[i]+s[i][j])*u[i][j] - (db[i]-p[i]+1)*(1-u[i][j])
		<= sum(t in r[j]+s[i][j]+1..db[j]-p[j]+1) t*x[j][t];
	}
	

	def1: forall(j in 1..n)
	{
		sum(t in 0..r[j]-1) x[j][t] == 0;
	}
	
	def2: forall(j in 1..n)
	{
		sum(t in ((db[j]-p[j])+1)+1..T) x[j][t]== 0;	
	}	
	
	def3:
	{
		x[0][0]==1;		
	}

		

}; 

