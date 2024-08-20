/*********************************************
 * OPL 12.6.0.0 Model
 * Author: Marcos
 * Creation Date: 25/03/2014 at 09:32:41
 *********************************************/
int NPeriods = ...;
int NGrammage = ...;//Number of grammage
int NSubPeriods = ...; //number of subperiods per period
int NSpeed = ...; //number of speed
int NPMachines = ...;//number of paper machines

range Periods = 1..NPeriods;
range Grammage = 1..NGrammage;
range SubPeriods = 1..NPeriods*NSubPeriods;
range PMachines = 1..NPMachines;

float price[Grammage]=...;


{int} S[t in Periods]={s | s in NSubPeriods*(t-1)+1..NSubPeriods*t};

range Speed = 1..NSpeed;
int speed_Inicial=...;//12.0
int grammage_Inicial[PMachines]=...;
int Wini[PMachines] = ...;

float  alfa=...; //parameter transformation 
float  vmax_dig=...; //maximum speed
float  vmin_dig=...; //minimum speed

float  Speeds[Speed]=...;//possible speeds

float Nmin=...;//Minimum size of a sub-period
float  hour[Periods]=...;//number of hours in each period

//Pulp towers
float  I0_virg=...; //initial inventory of virgin pulp
float  Imin_virg=...; //Minimum inventory of virgin pulp
float  Imax_virg=...; //maximum inventory of virgin pulp

//Paper machine
float  sl[Grammage][Grammage][PMachines]=...; //setup of changeover grammage k to j in ton
float  st[Grammage][Grammage][PMachines]=...; // in min transf hor setup of time changeover grammage k to j in 
float  b_virg[Grammage][PMachines]=...; //percentage of virgin pulp in gammage j
float  b[Grammage][PMachines]=...; //time to produce one ton of grammage j
float  M[Grammage][PMachines]=...; //upper bound on the production amounts
float  m_g[Grammage][PMachines]=...; //lower bound on the production amounts
//int Zmax=...;//Maximum number of changeover in each paper machine for each period

string NomesGrammage[Grammage] =...;
int wG[Grammage] = ...;
string tG[Grammage] =...;
float  d[Grammage][Periods]=...; //demand
float  Atrasos[Grammage]=...;//Initial Backlog

float demandaTotal = sum(i in Grammage)(Atrasos[i]+sum(t in Periods)(d[i][t]));

//Power Plant
float  delta = ...; //transformation pulp to black weak liquor (production proportion)
float I0_liquor=...;//Initial inventory of black liquor
float  Imax_liquor=...; //maximum inventory of black liquor
float  Imin_liquor=...; //maximum inventory of black liquor
float  mii = ...; //transform black weak liquor into steam
float fator = ...;
float fator2 = ...;
float I0_dry = ...;
float Imax_dry = ...;
float Imin_dry = ...;
float zeta = ...;
float  f_heat = ...; //maximum steam production
float gama = ...;

//digestor and recycled pulp mill
int teta = ...;
dvar float+ Xdig[SubPeriods];
dvar float+ Nhour[SubPeriods];
dvar int+ Yv[Speed][SubPeriods] in 0..1;//escolhas velocidades
dvar float+ Yv0[Speed] in 0..1;//gramagens iniciais 
dvar float+ Nh_s[Speed][SubPeriods];

//Pulp Towers
dvar float+ Ivirg[SubPeriods];
dvar float+ Ovirg[SubPeriods][PMachines];

//Paper machine
dvar float+ Z[Grammage][Grammage][SubPeriods][PMachines] in 0..1;
dvar int+ Y[Grammage][SubPeriods][PMachines] in 0..1;
dvar float+ X[Grammage][SubPeriods][PMachines];
dvar float+ IG[Grammage][0..NPeriods];
dvar float+ IGBacklog[Grammage][0..NPeriods];
dvar float+ Y0[Grammage][PMachines] in 0..1;//gramagens iniciais 
{int} NTop[PMachines] = ...;
{int} NBot[PMachines] = ...;
{int} A[PMachines] = ...;

//Power Plant
dvar float+ Xliquor[SubPeriods];
dvar float+ Iliquor[SubPeriods];
dvar float+ Oliquor[SubPeriods];
dvar float+ Xdry[SubPeriods];
dvar float+ Idry[SubPeriods];
dvar float+ Odry[SubPeriods];
dvar float+ XOdry[SubPeriods];
dvar float+ Oheat[SubPeriods];

//objective function
int NParc=4;
dvar float Parc[0..NParc];
dvar float+ dV[SubPeriods];

dvar float f;

minimize
10;//f;


subject to {

f == Parc[0] + Parc[1] + Parc[2] - Parc[3] + Parc[4];  


Parc[0] == sum(j in Grammage, t in Periods) 0.25*price[j]/365*IG[j][t];
Parc[1] == sum(j in Grammage, t in Periods) 2.5*price[j]/365*IGBacklog[j][t]; 
Parc[2] == sum(j in Grammage, k in Grammage, s in SubPeriods, m in PMachines) price[j]/5*sl[k][j][m]*Z[k][j][s][m];
Parc[3] == sum(s in SubPeriods) 0.1*Oheat[s];
Parc[4] == sum(s in SubPeriods) 0.01*dV[s];

/*******************************************/
//Disgester constraints
/*******************************************/
c1:
forall(v in Speed:Speeds[v]!=speed_Inicial)
  Yv0[v] == 0;
c1b:
forall(v in Speed:Speeds[v]==speed_Inicial)
  Yv0[v] == 1;
c2: 
forall(s in SubPeriods) 
  sum(v in Speed) Yv[v][s]==1;

cadd2z:  
forall(v in Speed)  
   Nhour[1]>=(Yv[v][1]-Yv0[v])*Nmin;
cadd2:
forall(v in Speed, s in 2..NPeriods*NSubPeriods)  
   Nhour[s]>=(Yv[v][s]-Yv[v][s-1])*Nmin;

c3:    
forall(t in Periods, s in S[t], v in Speed) 
  Nh_s[v][s]<=hour[t]*Yv[v][s]; 
c4:
forall(s in SubPeriods) 
  Xdig[s] == alfa*sum(v in Speed) Speeds[v]* Nh_s[v][s];

c5: 
forall(s in SubPeriods) 
  Nhour[s]==sum(v in Speed) Nh_s[v][s];

c6:    
forall(t in Periods)    
  sum(s in S[t]) Nhour[s]==hour[t];

c7z:        
forall(v in Speed)
  Yv[v][1]<=sum(v1 in v-teta..v+teta: v1>=1 && v1<=NSpeed)Yv0[v1];
c7:
forall(s in 2..NSubPeriods*NPeriods, v in Speed)
  Yv[v][s]<=sum(v1 in v-teta..v+teta: v1>=1 && v1<=NSpeed)Yv[v1][s-1];
c8:
forall(s in 2..NPeriods*NSubPeriods, v in Speed){
  dV[s] >= Yv[v][s]-Yv[v][s-1];
  dV[s] >= Yv[v][s-1]-Yv[v][s];
}  
c9:
forall(v in Speed)
  dV[1] >= Yv[v][1]-Yv0[v];
c10:
forall(v in Speed)
  dV[1] >= Yv0[v]-Yv[v][1];
/*******************************************/
//Virgin pulp inventory
/*******************************************/
c12z:    
Xdig[1]+I0_virg == sum(m in PMachines) Ovirg[1][m]+Ivirg[1];
c12:
forall (s in 2..NSubPeriods*NPeriods)
  Xdig[s]+Ivirg[s-1] == sum(m in PMachines) Ovirg[s][m]+Ivirg[s];  
c14:   	
forall (s in SubPeriods) 
  Imin_virg<=Ivirg[s]<=Imax_virg;
/*******************************************/
//Paper machines constraints
/*******************************************/
c15: 
forall(s in SubPeriods, m in PMachines) 
  sum(j in Grammage) (b_virg[j][m]*X[j][s][m])+
  sum(k in Grammage, j in Grammage) (b_virg[j][m]*sl[k][j][m]*Z[k][j][s][m])==Ovirg[s][m];
c17: 
forall (s in SubPeriods, m in PMachines) 
  sum  (j in Grammage) (b[j][m]*X[j][s][m])+
  sum (k in Grammage, j in Grammage) ((1/60)*st[k][j][m]*Z[k][j][s][m])==Nhour[s];
c18:
forall(m in PMachines) 
	sum(i in Grammage) Y0[i][m]==1;
c19:
forall(m in PMachines)
  Y0[grammage_Inicial[m]][m]==1;  
c20:    
 forall (s in SubPeriods, m in PMachines) 
       sum  (i in Grammage) Y[i][s][m]==1;    
c21: 
forall (s in SubPeriods, m in PMachines) 
    sum  (i,j in Grammage) Z[i][j][s][m]==1;
c22:
forall(k in Grammage, m in PMachines)
   sum(j in Grammage) Z[k][j][1][m] == Y0[k][m];

c23:
forall(k in Grammage, s in 2..NPeriods*NSubPeriods, m in PMachines)
   sum(j in Grammage) Z[k][j][s][m] == Y[k][s-1][m];
c24:
forall(j in Grammage, s in SubPeriods, m in PMachines)
   sum(k in Grammage) Z[k][j][s][m] == Y[j][s][m];

c25:
forall (j in Grammage, t in Periods, s in S[t], m in PMachines)
   X[j][s][m]<=M[j][m]*Y[j][s][m];

c25a:
forall (j in Grammage, t in Periods, s in S[t], m in PMachines)
   X[j][s][m]<=hour[t]/b[j][m]*Y[j][s][m];


c26z:   	
forall (j in Grammage, m in PMachines)
  X[j][1][m]>=m_g[j][m]*(Y[j][1][m]-Y0[j][m]);
c26:   	
forall (j in Grammage, s in 2..NSubPeriods*NPeriods, m in PMachines)
  X[j][s][m]>=m_g[j][m]*(Y[j][s][m]-Y[j][s-1][m]);

c27:
forall (t in Periods,j in Grammage) 
   sum  (s in S[t], m in PMachines) X[j][s][m]-d[j][t]+ IG[j][t-1]+IGBacklog[j][t]-IGBacklog[j][t-1]==IG[j][t];
c28:
forall (j in Grammage)
	IGBacklog[j][0]==Atrasos[j];
c29:
forall (j in Grammage)
	IG[j][0]==0;
/*******************************************/
//Recovery constraints
/*******************************************/

c30:
   forall (s in SubPeriods) 
      Xliquor[s] == delta*sum(v in Speed) Speeds[v]* Nh_s[v][s];

c31z:    
  Xliquor[1]+I0_liquor == Oliquor[1]+Iliquor[1];

c31:	
forall (s in 2..NSubPeriods*NPeriods)
  Xliquor[s]+Iliquor[s-1] == Oliquor[s]+Iliquor[s]; 

forall(s in SubPeriods)
  Imin_liquor <= Iliquor[s] <= Imax_liquor;

forall (s in SubPeriods) 
  Oliquor[s]<=  mii*Nhour[s];
              
forall (s in SubPeriods) 
  Xdry[s]== gama*Oliquor[s];

Xdry[1]+I0_dry == Odry[1]+Idry[1];
forall (s in 2..NSubPeriods*NPeriods)
  Xdry[s]+Idry[s-1] == Odry[s]+Idry[s]; 
      
forall (s in SubPeriods) 
  Imin_dry <= Idry[s]<= Imax_dry;

forall(s in SubPeriods)
  XOdry[s]==fator*Odry[s];

forall (s in SubPeriods) 
  XOdry[s] <= fator2*Nhour[s]; 
 
forall (s in SubPeriods) 
  Oheat[s]==zeta*XOdry[s]; 
        
forall (s in SubPeriods) 
  Oheat[s] <= f_heat*Nhour[s];   

/*******************************************/
//Aditional Cuts
/*******************************************/
}
tuple tripla{
key int i;
key int m;
float value;
}

{tripla} b_aux = {<i,m,b[i][m]>|i in Grammage, m in PMachines};
{tripla} bv_aux = {<i,m,b_virg[i][m]>|i in Grammage, m in PMachines};

execute{
	writeln(b_aux)
	writeln(bv_aux)
}
