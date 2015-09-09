#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

/* scip includes */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* user defined includes */
#include "pricer_CapacityExp.h"


/* namespace usage */
using namespace std;
using namespace scip;

void GenerateDemand (vector<int> &x , vector<double> &y);
double CalProb (int n, double mean);
int CalMin (int  a, int b);
double mean=5;
double K=2;


int main ()
{

  vector<double> purchcost = {1,4,7};
  vector<int >vialsize = {1,5,10};
  double ConfRate=0.75;
  double Shortbeta=1;
  double Wastbeta=1;
  int I=3;
  double M=100;
  int Nper_to_spoil=1;
  int T=3; //T is the number of stages
  int NScen=pow(K, T-1);
  int Nnodes=0;
  for (int i=0; i<T ; i++){
	Nnodes=Nnodes+pow(K,i);}

  vector<double>Dprob;
  vector<int>Demand;
  vector<double>Nodedemand;
  vector<int>ScenDemand;
  vector<double>Scenprob;
  vector<double>Scenprob2;
  vector<double>Nodeprobe;

  GenerateDemand(Demand, Dprob);

  //generating the nodes matrix
  vector< vector<int> > nodemat;
  vector< vector<int> > jmat;

	int n=1;
	for (int t=1 ; t<T+1 ; t++){
		vector<int>row; //creat an empty vector
                vector<int>jrow;
		int len=NScen/pow(K,t-1);
		for (int i=0 ; i<pow(K,t-1) ; i++){
			for (int j=0 ; j<len ; j++){
				row.push_back(n);
                                jrow.push_back(i);}
		n++;}
		nodemat.push_back(row);
                jmat.push_back(jrow);}


	for (int i=0 ; i<pow(K,T-1) ; i++){
		Scenprob2.push_back(0);
	    Scenprob.push_back(0);}


	for (int j=0 ; j<K ; j++){
		Nodedemand.push_back(Demand[j]);
		Nodeprobe.push_back(Dprob[j]);
		Scenprob[j]=Dprob[j];}

	for (int t=3; t<T+1 ; t++){
		int n=0;
		for (int i=0 ; i<pow(K,t-2) ; i++){
			for (int j=0 ; j < K ; j++){
				Nodedemand.push_back(Demand[j]);
				Nodeprobe.push_back(Dprob[j]);
				Scenprob2[n]=Scenprob[i]* Dprob[j];
				n++;}}
		int ssize= (int)Scenprob2.size();
		for(int ii=0 ; ii<ssize ; ii++){
			Scenprob[ii]=Scenprob2[ii];}
	}

	//Generating ScenDemand
	for (int i=0 ; i<pow(K,T-2) ; i++){
		for (int j=0; j<K ; j++){
			ScenDemand.push_back(Demand[j]);}}
  SCIP* scip = NULL;

  static const char* vaccine_PRICER_NAME = "vaccine_Pricer";

  SCIP_CALL( SCIPcreate(&scip) );

   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* set verbosity parameter */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetBoolParam(scip, "display/lpinfo", TRUE) ); 

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip, "vaccine", 0, 0, 0, 0, 0, 0, 0) );
  
  /*Add variables*/
  char var_name[255];
  vector< SCIP_VAR* > z_var(NScen);
  for (int i=0; i<NScen; i++)
  {
     SCIP_VAR* var;
     SCIPsnprintf(var_name, 255, "Z%d", i );
     SCIP_CALL( SCIPcreateVar(scip,
                     &var,                   // returns new index
                     var_name,               // name
                     0.0,                    // lower bound
                     1.0,                    // upper bound
                     0,                      // objective
                     SCIP_VARTYPE_CONTINUOUS,   // variable type
                     true,                   // initial
                     false,                  // forget the rest ...
                     0, 0, 0, 0, 0) );
    SCIP_CALL( SCIPaddVar(scip, var) );
    z_var[i] = var;
   }

  /* Add w variables*/
  vector< vector< SCIP_VAR * > > wj_var(T-1);
  for (int t=2; t<T+1 ; t++)
  {
   wj_var[t-2].resize(pow(K,t-1));
   int lenght  = NScen/ pow(K,t-1);
   for (int j=0 ; j<pow(K,t-1); j++)
   { SCIP_Real objj = Nodeprobe[nodemat[t-1][j*lenght]-2]*purchcost[1]*1;
    SCIP_VAR * var;
    SCIPsnprintf(var_name, 255, "wt%d_%d", t, j );
    SCIP_CALL( SCIPcreateVar(scip,
                     &var,                   // returns new index
                     var_name,               // name
                     0.0,                    // lower bound
                     SCIPinfinity(scip),     // upper bound
                     objj,                      // objective
                     SCIP_VARTYPE_CONTINUOUS,   // variable type
                     true,                   // initial
                     false,                  // forget the rest ...
                     0, 0, 0, 0, 0) );
    SCIP_CALL( SCIPaddVar(scip, var) );
    wj_var[t-2][j] = var;
    }
  }
  /* Add w variables*/
  /*vector< vector< SCIP_VAR * > > wj_var2(T-1);
  for (int t=2; t<T+1 ; t++)
  {
   wj_var2[t-2].resize(pow(K,t-1));
   int lenght  = NScen/ pow(K,t-1);
   for (int j=0 ; j<pow(K,t-1); j++)
   { SCIP_Real objj = Nodeprobe[nodemat[t-1][j*lenght]-2]*purchcost[0]*Nodedemand[nodemat[t-1][j*lenght]-2];
    SCIP_VAR * var;
    SCIPsnprintf(var_name, 255, "wt%d_%d", t, j );
    SCIP_CALL( SCIPcreateVar(scip,
                     &var,                   // returns new index
                     var_name,               // name
                     0.0,                    // lower bound
                     SCIPinfinity(scip),    // upper bound
                     objj,                      // objective
                     SCIP_VARTYPE_CONTINUOUS,   // variable type
                     true,                   // initial
                     false,                  // forget the rest ...
                     0, 0, 0, 0, 0) );
    SCIP_CALL( SCIPaddVar(scip, var) );
    wj_var2[t-2][j] = var;
    }
  }*/
  
  /*Add constraint 5b */
  vector< vector < vector <SCIP_CONS* > > > sh_con(T-1);
 
 
  char con_name[255];
  for(int t=2; t<T+1; t++)
  {
    sh_con[t-2].resize(pow(K,t-1));
    int lenght=NScen/pow(K,t-1); 
    for (int j=0 ; j<pow(K,t-1) ; j++)
      {
         sh_con[t-2][j].resize(pow(K,T-t));
         for (int q=0 ; q<pow(K,T-t) ; q++)
            {
                SCIP_Real mwj =Nodedemand[nodemat[t-1][j*lenght]-2] ; 
                SCIP_Real mwj2 = 0;  
                SCIP_CONS* con;
                SCIPsnprintf(con_name, 255, "sh%d_%d_%d", t, j, q);
                SCIP_VAR* index = z_var[j * pow(K,T-t) + q];
                SCIP_Real coeff = Nodedemand[nodemat[t-1][j*lenght]-2];
                SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 1, &index, &coeff,
                     Nodedemand[nodemat[t-1][j*lenght]-2]-Shortbeta,    /* lhs */
                     SCIPinfinity(scip),     /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     true,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */
                SCIP_CALL( SCIPaddCons(scip, con) );
                sh_con[t-2][j][q] = con;
                if (t==T && j==0)
                   SCIP_CALL( SCIPaddCoefLinear(scip, sh_con[t-2][j][q], wj_var[t-2][j], mwj2) );
                else
                   SCIP_CALL( SCIPaddCoefLinear(scip, sh_con[t-2][j][q], wj_var[t-2][j], mwj) );
               // SCIP_CALL( SCIPaddCoefLinear(scip, sh_con[t-2][j][q], wj_var2[t-2][j], mwj) );
             }
       }
   }
   /*Add constraint 5c*/
   vector<vector<SCIP_CONS*> >z_con (NScen); 
   SCIP_Real mwjj = 0;
   SCIP_Real mwjj2 = -1;
   //char con_name[255];  
   for (int i=0; i<NScen; i++)
   {  
      z_con[i].resize(T-1);
      for (int j=2; j<T+1 ; j++)
       {
          SCIP_CONS* con;
          SCIPsnprintf(con_name, 255, "z%d_%d", i, j);
          SCIP_VAR* index = z_var[i];
          SCIP_Real coeff = 1;
          SCIP_CALL( SCIPcreateConsLinear(scip, &con, con_name, 1, &index, &coeff,
                     0,                      /* lhs */
                     SCIPinfinity(scip),     /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     true,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */
                SCIP_CALL( SCIPaddCons(scip, con) );
                z_con[i][j-2] = con;
                int jjjmat= jmat [j-1][i];
                if (j==T && i==0)
                   SCIP_CALL( SCIPaddCoefLinear(scip, z_con[i][j-2], wj_var[j-2][jjjmat], mwjj2) );
                else
                   SCIP_CALL( SCIPaddCoefLinear(scip, z_con[i][j-2], wj_var[j-2][jjjmat], mwjj) );           
                //SCIP_CALL( SCIPaddCoefLinear(scip, z_con[i][j-2], wj_var2[j-2][jjjmat], mwjj) );
        }
    }

   /*Add constraint 5d */
   //char con_name[255];
   SCIP_CONS*  alpha_con = NULL;
   SCIP_CALL( SCIPcreateConsLinear(scip, & alpha_con,"zn" , 0, NULL, NULL,
                     -SCIPinfinity(scip),    /* lhs */
                     1-ConfRate,                    /* rhs */
                     true,                   /* initial */
                     false,                  /* separate */
                     true,                   /* enforce */
                     true,                   /* check */
                     true,                   /* propagate */
                     false,                  /* local */
                     true,                   /* modifiable */
                     false,                  /* dynamic */
                     false,                  /* removable */
                     false) );               /* stickingatnode */

    for (int i = 0; i < NScen; ++i)
	 SCIP_CALL( SCIPaddCoefLinear(scip, alpha_con, z_var[i], Scenprob[i]) );
    SCIP_CALL( SCIPaddCons(scip, alpha_con) );

   /*Add constraint 5e */
   //char con_name[255]; 
   SCIP_Real wje =1;
   vector< vector<SCIP_CONS* > > w_con(T-1);
   for (int t=2; t< T+1 ; t++)
   {
    w_con[t-2].resize(pow (K, t-1));
    for (int j=0; j < pow (K, t-1); j++)
     {
      SCIP_CONS* con = NULL;
      SCIPsnprintf(con_name, 255, "w%d_%d", t, j);
      SCIP_CALL( SCIPcreateConsLinear( scip, &con, con_name, 0, NULL, NULL,
                                       1.0,   /* lhs */
                                       1.0,   /* rhs */
                                       true,  /* initial */
                                       false, /* separate */
                                       true,  /* enforce */
                                       true,  /* check */
                                       true,  /* propagate */
                                       false, /* local */
                                       true,  /* modifiable */
                                       false, /* dynamic */
                                       false, /* removable */
                                       false  /* stickingatnode */ ) );
      SCIP_CALL( SCIPaddCons(scip, con) );
      w_con[t-2][j] = con;
      SCIP_CALL( SCIPaddCoefLinear(scip, w_con[t-2][j], wj_var[t-2][j],wje ) );
      //SCIP_CALL( SCIPaddCoefLinear(scip, w_con[t-2][j], wj_var2[t-2][j],wje ) );
     }
   }
         
   /* include pricer */
   ObjPricervaccine* vaccine_pricer_ptr = new ObjPricervaccine(scip, vaccine_PRICER_NAME, purchcost, vialsize, Nper_to_spoil,        
 Wastbeta, T, NScen, Nnodes, K, I, Nodeprobe, nodemat, jmat, z_var, wj_var, sh_con, z_con, alpha_con, w_con);



   SCIP_CALL( SCIPincludeObjPricer(scip, vaccine_pricer_ptr, true) );

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, vaccine_PRICER_NAME)) );

   SCIP_CALL( SCIPwriteOrigProblem(scip, "vCapacityExp_init.lp", "lp", FALSE) );

   /*************
    *  Solve    *
    *************/

  SCIP_CALL( SCIPsolve(scip) );


   /**************
    * Statistics *
    *************/
   //FILE* ffp =fopen("/home/bahar/SCIPprpject/vaccine/resultt","w");
   FILE* ffp2 =fopen("/home/bahar/SCIPproject/vaccine/result2","w");
   SCIP_CALL( SCIPprintStatistics(scip, ffp2) );
  

   SCIP_CALL( SCIPprintBestSol(scip, ffp2, FALSE) );



   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return 0;
}//for main 

void GenerateDemand (vector<int> &x , vector<double> &y){
	int a, b;
        a = CalMin(mean,ceil((K-1)/2));
        if (mean > ceil ((K-1)/2))
           b= floor((K-1)/2);
        else
           b= K-mean-1;
        double prob = 0;
        for (int i= 0; i< mean-a +1 ; i++){
          prob = prob +CalProb (i , mean);
          }
        double probb = prob;
        x.push_back(mean-a);
        y.push_back(prob);
	for (int i=1 ; i<a+b ; i++){
           y.push_back(CalProb(mean-a+i , mean)); 
	   x.push_back(mean-a+i);
           probb = probb +y[i];}
        x.push_back (mean+b);
        y.push_back(1-probb);
}
double CalProb (int n, double mean) 
 {
  double prob=exp(-1*mean);
  if (n > 0)
  {
    for (int i=1 ; i<n+1 ; i++)
    {
     prob = prob * mean / i ;
    }
   } 
  return (prob);
 }

int CalMin (int  a, int b)
{
   int min;
   if (a < b)
     min=a;
   else
     min=b;
   return (min);
}
