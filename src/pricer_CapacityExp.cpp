#include "pricer_CapacityExp.h"
#include "scip/scipdefplugins.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "scip/cons_linear.h"


using namespace std;
using namespace scip;

#define PRICER_NAME            "vaccine"
#define PRICER_DESC            "variable pricer template"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */


/** Constructs the pricer object with the data needed
 *
 *  An alternative is to have a problem data class which allows to access the data.
 */
ObjPricervaccine::ObjPricervaccine(
		SCIP*                                            scip,          /**< SCIP pointer */
		const char*                                      p_name,        /**< name of pricer */
		vector<double> &                                 p_purchcost,
		vector<int > &                                   p_vialsize,
		int                                              p_Nper_to_spoil,
		double                                           p_Wastbeta,
		int                                              p_T,
		int                                              p_NScen,
		int                                              p_Nnodes,
		int                                              p_K,
		int                                              p_I,
		vector<double> &                                 p_Nodeprobe,
		vector< vector<int> > &                          p_nodemat,
		vector< vector<int> > &                          p_jmat,
		vector< SCIP_VAR* > &                            p_z_var,
		vector< vector< SCIP_VAR * > > &                 p_wj_var,
              //  vector< vector< SCIP_VAR * > > &                 p_wj_var2,
		vector< vector < vector <SCIP_CONS* > > > &      p_sh_con,
		vector<vector<SCIP_CONS* > > &                   p_z_con,
		SCIP_CONS*                                     & p_alpha_con,
		vector< vector<SCIP_CONS* > > &                   p_w_con
):
				   ObjPricer(scip, p_name, "Finds tour with negative reduced cost.", 0, TRUE),
				   _purchcost(p_purchcost),
				   _vialsize(p_vialsize),
				   _Nper_to_spoil(p_Nper_to_spoil),
				   _Wastbeta(p_Wastbeta),
				   _T(p_T),
				   _NScen(p_NScen),
				   _Nnodes(p_Nnodes),
				   _K(p_K),
				   _I(p_I),
				   _Nodeprobe(p_Nodeprobe),
				   _nodemat(p_nodemat),
				   _jmat(p_jmat),
				   _z_var(p_z_var),
				   _wj_var(p_wj_var),
				//   _wj_var2(p_wj_var2),
				   _sh_con(p_sh_con),
				   _z_con(p_z_con),
				   _alpha_con(p_alpha_con),
				   _w_con(p_w_con)
{}


/** Destructs the pricer object. */
ObjPricervaccine::~ObjPricervaccine()
{}

SCIP_DECL_PRICERINIT(ObjPricervaccine::scip_init)
{
	for (int i=0; i<_NScen; i++)
	{  //gets the transformed var and z_con
		SCIP_CALL( SCIPgetTransformedVar(scip, _z_var[i], &_z_var[i]) );
		for (int j=2; j<_T+1 ; j++)
		{
			SCIP_CALL( SCIPgetTransformedCons(scip, _z_con[i][j-2], &_z_con[i][j-2]) );
		}
	}


	for(int t=2; t<_T+1; t++)
	{
		for (int j=0 ; j<pow(_K,t-1) ; j++)
		{
			SCIP_CALL( SCIPgetTransformedCons(scip, _w_con[t-2][j], &_w_con[t-2][j]));
			SCIP_CALL( SCIPgetTransformedVar(scip, _wj_var[t-2][j], &_wj_var[t-2][j]) );
            //SCIP_CALL( SCIPgetTransformedVar(scip, _wj_var2[t-2][j], &_wj_var2[t-2][j]) );
			for (int q=0 ; q<pow(_K,_T-t) ; q++)
			{
				SCIP_CALL( SCIPgetTransformedCons(scip, _sh_con[t-2][j][q], &_sh_con[t-2][j][q]) );
			}
		}
	}
	//gets transformed for alpha_con
	SCIP_CALL( SCIPgetTransformedCons(scip, _alpha_con, &_alpha_con));

	return SCIP_OKAY;
}


/** perform pricing*/
SCIP_RETCODE ObjPricervaccine::pricing(SCIP*   scip)               /**< SCIP data structure */
{ //Get the duals of s_sh_con
	vector< vector < vector <SCIP_Real > > > pi(_T-1);
	for(int t=2; t<_T+1; t++)
	{
		pi[t-2].resize(pow(_K,t-1));
		for (int j=0 ; j<pow(_K,t-1) ; j++)
		{
			pi[t-2][j].resize(pow(_K,_T-t));
			for (int q=0 ; q<pow(_K,_T-t) ; q++)
			{
				pi[t-2][j][q] = SCIPgetDualsolLinear(scip, _sh_con[t-2][j][q]);
				cout << "pi" << pi[t-2][j][q] << endl;
			}
		}
	}

	//Get the duals of 5c*/
	vector<vector<SCIP_Real>> gamma(_NScen);

	for (int i=0; i<_NScen; i++)
	{
		gamma[i].resize(_T-1);
		for (int j=2; j<_T+1 ; j++)
		{
			gamma[i][j-2] = SCIPgetDualsolLinear(scip, _z_con[i][j-2]);
			cout << "gamma" << gamma[i][j-2] << endl;
		}
	}
	//Gets the dual of w_con
	vector < vector < SCIP_Real > > mu(_T-1);
	for (int t=2; t< _T+1 ; t++)
	{
		mu[t-2].resize(pow(_K,t-1));
		for (int j=0; j < pow (_K, t-1); j++)
		{
			mu[t-2][j]= SCIPgetDualsolLinear(scip, _w_con[t-2][j]);
			cout << "mu" << mu[t-2][j] << endl;
		}
	}


	vector< vector <vector< SCIP_Real > > > x_newcolumn (_I);
	vector< vector <vector< vector< vector< SCIP_Real > > > > > y_newcolumn (_I);
	vector< vector < SCIP_Real > > z_newcolumn (_T-1);
	vector< vector < SCIP_Real > > objval (_T-1);
	vector< vector < SCIP_Real > > redcost(_T-1);
	for (int i=0; i<_I ; i++)
	{
		x_newcolumn[i].resize(_T-1);
		y_newcolumn[i].resize(_T-1);
		for (int t=2; t<_T+1; t++)
		{
			x_newcolumn[i][t-2].resize(pow(_K, t-1));
			y_newcolumn[i][t-2].resize(pow(_K, t-1));
			for (int j=0; j<pow(_K, t-1); j++)
			{
				int YNN;
				if (_T-t < _Nper_to_spoil)
					YNN = _T-t;
				else
					YNN = _Nper_to_spoil;
				y_newcolumn[i][t-2][j].resize(YNN+1);
				for (int kk=0; kk < YNN+1 ; kk++)
				{
					y_newcolumn[i][t-2][j][kk].resize(pow(_K,kk));
				}
			}
		}//t
	}// y[i][t][j][kk][jj]
	for (int t=2; t<_T+1; t++)
	{
		z_newcolumn[t-2].resize(pow(_K, t-1));
		objval[t-2].resize(pow(_K,t-1));
		redcost[t-2].resize(pow(_K,t-1));
	}

	SCIP_CALL( find_new_column(pi , gamma,  x_newcolumn, y_newcolumn, z_newcolumn, objval));
	int n_added_var =0;
	vector< SCIP_Real > xcolumn (_I);
	vector< vector< vector < SCIP_Real > > > ycolumn (_I);
	for (int t=2; t<_T+1; t++)
	{
		for (int j=0 ; j<pow(_K,t-1); j++)
		{
			int NN;
			if (_T-t < _Nper_to_spoil)
				NN = _T-t;
			else
				NN = _Nper_to_spoil;

			for (int i=0 ; i<_I; i++)
			{
				ycolumn[i].resize(NN+1);
				xcolumn[i] = x_newcolumn[i][t-2][j];
				for(int kk=0 ; kk<NN+1; kk++)
				{
					ycolumn[i][kk].resize(pow(_K, kk));
					for(int jj=0 ;jj<pow(_K, kk); jj++)
					{
						ycolumn[i][kk][jj] = y_newcolumn[i][t-2][j][kk][jj];
					}
				}
			}

			redcost[t-2][j] = objval[t-2][j] - mu[t-2][j];
			cout << "objval" << objval[t-2][j] << endl;
			if ( SCIPisNegative(scip, redcost[t-2][j]) )
			{
				SCIP_CALL(add_newcolumn_variable(scip, xcolumn, ycolumn, z_newcolumn[t-2][j], t, j));
				n_added_var++;
			}
		}
	}
	if (n_added_var > 0 )
	{
		return SCIP_OKAY;
	}


	//see if reduced cost is negative
	/* add tour variable */

	SCIP_CALL( SCIPwriteTransProblem(scip, "CapacityExp.lp", "lp", FALSE) );


	return SCIP_OKAY;
}

SCIP_DECL_PRICERREDCOST(ObjPricervaccine::scip_redcost)
{
	SCIPdebugMessage("call scip_redcost ...\n");

	/* set result pointer, see above */
	*result = SCIP_SUCCESS;

	/* call pricing routine */
	SCIP_CALL( pricing(scip) );

	return SCIP_OKAY;
}

/** add new variable to problem */
SCIP_RETCODE ObjPricervaccine::add_newcolumn_variable(
		SCIP*                                                           scip,
		vector< SCIP_Real >                                           & x_newcolumn,
		vector< vector< vector < SCIP_Real > > >                      & y_newcolumn,
		SCIP_Real                                                      z_newcolumn,
		int                                                             t,
		int                                                             j )        /**< list of nodes in tour */
{
	SCIP_Real wobjcoef =0 ;
	char var_name[255];
	SCIP_VAR * w_var;
	int lenght = pow(_K,_T-1)/pow(_K, t-1);
	for (int i=0 ; i<_I ; i++)
	{
		//cout << x_newcolumn[i] << endl;
		wobjcoef = wobjcoef +  x_newcolumn[i] * _purchcost[i] * _Nodeprobe[_nodemat[t-1][j*lenght]-2];
	}
	//cout << "wobjcoef" << wobjcoef << endl;
	SCIPsnprintf(var_name, 255, "Alpha%d_%d", t, j);
	SCIP_CALL( SCIPcreateVar(scip, & w_var, var_name,
			0.0,                      // lower bound
			SCIPinfinity(scip),      // upper bound
			wobjcoef,               // objective
			SCIP_VARTYPE_CONTINUOUS, // variable type
			false, false, 0, 0, 0, 0, 0) );
	SCIPdebugMessage("new variable <%s>\n", var_name);

	/* add new variable to the list of variables to price into LP (score: leave 1 here) */
	SCIP_CALL( SCIPaddPricedVar(scip, w_var, 1.0) );

	// Add the new variable to the constraint 5b
	int NNN;
	if (_T-t < _Nper_to_spoil)
		NNN = _T-t;
	else
		NNN = _Nper_to_spoil;
	for (int tt=t; tt< t+NNN+1 ; tt++)
	{
		int lenn= pow(_K,tt-1)/pow(_K,t-1);
		for (int jj=0 ; jj <lenn ; jj++)
		{
			for (int q=0 ; q< pow(_K, _T-tt) ; q++)
			{
				SCIP_Real coefsh = 0;
				for (int i=0 ; i<_I; i++)
				{
					coefsh = coefsh + y_newcolumn [i][tt-t][jj];
				}
				//cout << "coefsh" << coefsh << endl;
				cout << _sh_con[tt-2][j*lenn+jj][q] << endl;
				SCIP_CALL( SCIPaddCoefLinear(scip, _sh_con[tt-2][j*lenn+jj][q], w_var, coefsh ) );
			}
		}
	}

	// Add to constraint 5c
	for (int jj=0 ; jj< pow(_K,_T-1)/pow(_K,t-1); jj++)
	{
		SCIP_CALL( SCIPaddCoefLinear(scip, _z_con[j*pow(_K,_T-1)/pow(_K,t-1)+jj][t-2], w_var, -1*z_newcolumn ) );
	}

	//Add to constraint 5e
	SCIP_CALL( SCIPaddCoefLinear(scip, _w_con[t-2][j], w_var, 1 ) );

	//cleanup

	SCIP_CALL( SCIPreleaseVar(scip, &w_var) );

	return SCIP_OKAY;
}

//Finds new columns
SCIP_RETCODE ObjPricervaccine::find_new_column(
		vector< vector < vector <SCIP_Real > > > pi,      /**< dual variable value*/
		vector< vector < SCIP_Real > > gamma,
		vector< vector <vector< SCIP_Real > > > & x_newcolumn,
		vector< vector <vector< vector< vector< SCIP_Real > > > > > & y_newcolumn,
		vector< vector < SCIP_Real > > & z_newcolumn,
		vector< vector < SCIP_Real > > & objval )
{

	vector < SCIP_Real > x_ub = {50,10,5};
	for (int t=2 ; t<_T+1; t++)
	{
		int lenght = _NScen/pow(_K,t-1);
		for (int j=0; j<pow(_K,t-1) ; j++)
		{
			SCIP* subscip = NULL;
			SCIP_CALL( SCIPcreate(&subscip) );
			SCIPprintVersion(subscip, NULL);
			SCIPinfoMessage(subscip, NULL, "\n");

			/* include default plugins */
			SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

			/* set verbosity parameter */
			SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
			//SCIP_CALL( SCIPsetBoolParabm(subscip, "display/lpinfo", TRUE) );

			/* create empty problem */
			SCIP_CALL( SCIPcreateProb(subscip, "p_vaccine", 0, 0, 0, 0, 0, 0, 0) );

			/*Add xi variables*/
			char var_name[255];
			vector< SCIP_VAR* > x_var(_I);
			for (int i = 0; i <_I; ++i)
			{
				SCIP_VAR* xvar;
				SCIPsnprintf(var_name, 255, "x%d", i );
				SCIP_CALL( SCIPcreateVar(subscip,
						&xvar,                                                      // returns new index
						var_name,                                                   // name
						0.0,                                                        // lower bound
						x_ub[i],                                      // upper bound
						_Nodeprobe[_nodemat[t-1][j*lenght]-2] * _purchcost[i],      // objective
						SCIP_VARTYPE_CONTINUOUS,                                       // variable type
						true,                                                       // initial
						false,                                                      // forget the rest ...
						0, 0, 0, 0, 0) );
				SCIP_CALL( SCIPaddVar(subscip, xvar) );
				x_var[i] = xvar;
			}
			SCIPinfoMessage(subscip, NULL, "Number of x variables is sunbproblem is" , _I);

			vector< vector< vector <SCIP_Real> > > shcoeff(_I);
			/*Add variables Yi(m)k */
			//char var_name[255];
			vector< vector< vector <SCIP_VAR*> > > y_var(_I);
			for (int i=0 ; i<_I ; i++)
			{
				int NNN;
				if (_T-t < _Nper_to_spoil)
					NNN = _T-t;
				else
					NNN = _Nper_to_spoil;
				shcoeff[i].resize(NNN+1);
				y_var[i].resize(NNN+1);
				for (int kk=0; kk <NNN+1; kk++)
				{
					shcoeff[i][kk].resize(pow(_K,kk));
					y_var[i][kk].resize(pow(_K,kk));
					int len = pow(_K,t+kk-1)/pow(_K,t-1);
					for(int jj=0; jj < pow(_K,kk) ; jj++)
					{
						for (int q=0; q<pow(_K,_T-t-kk) ; q++)
						{
							shcoeff[i][kk][jj] -= pi[t+kk-2][j*len + jj][q];
						}
						//Add
						SCIP_VAR* yvar;
						SCIPsnprintf(var_name, 255, "y%d", i );
						SCIP_CALL( SCIPcreateVar(subscip,
								&yvar,                                                      // returns new index
								var_name,                                                   // name
								0.0,                                                        // lower bound
								50,                                      // upper bound
								shcoeff[i][kk][jj],                                         // objective
								SCIP_VARTYPE_CONTINUOUS,                                       // variable type
								true,                                                       // initial
								false,                                                      // forget the rest ...
								0, 0, 0, 0, 0) );
						SCIP_CALL( SCIPaddVar(subscip, yvar) );
						y_var[i][kk][jj] = yvar;
						SCIPinfoMessage(subscip, NULL, "New Y variable is Y_%d_%d_%d_%d_%d", i, t, j, t+kk, j * (pow(_K, t+kk-1)/pow(_K, t-1))+jj);
					} // jj
				}
			}

			/*Add variables z_m */
			SCIP_VAR* pz_var;
			SCIP_Real zcoef=0;
			for (int q=0; q<pow(_K,_T-t); q++)
			{
				zcoef += gamma[j*lenght+q][t-2];
			}
			SCIP_CALL( SCIPcreateVar(subscip,
					&pz_var,                                                      // returns new index
					var_name,                                                   // name
					0.0,                                                        // lower bound
					1,                                                          // upper bound
					zcoef,                                                      // objective
					SCIP_VARTYPE_CONTINUOUS,                                        // variable type
					true,                                                       // initial
					false,                                                      // forget the rest ...
					0, 0, 0, 0, 0) );
			SCIP_CALL( SCIPaddVar(subscip, pz_var) );


			/* Add te constraint 6b */
			char con_name[255];
			int NN;
			if (_Nper_to_spoil > _T-t)
				NN = _T-t;
			else
				NN = _Nper_to_spoil;
			vector< SCIP_CONS* > wt_con(pow(_K,NN));
			for (int n=0; n<pow(_K,NN); n++)
			{
				SCIP_CONS* ww_con = NULL;
				SCIPsnprintf(con_name, 255, "wt%d", n);
				SCIP_CALL( SCIPcreateConsLinear(subscip, & ww_con, con_name, 0, NULL, NULL,
						-SCIPinfinity(subscip),    /* lhs */
						_Wastbeta,                    /* rhs */
						true,                   /* initial */
						false,                  /* separate */
						true,                   /* enforce */
						true,                   /* check */
						true,                   /* propagate */
						false,                  /* local */
						true,                   /* modifiable */
						false,                  /* dynamic */
						false,                  /* removable */
						false) );               /* stickingatnode*/
				wt_con[n]=ww_con;
				for (int ii=0; ii<_I; ii++)
				{
					SCIP_CALL( SCIPaddCoefLinear(subscip, wt_con[n], x_var[ii], _vialsize[ii]) );
					SCIPinfoMessage(subscip, NULL, "X at this constraint X_%d", ii);
					for(int kk=0; kk < NN+1 ; kk++)
					{
						int leng = _NScen/pow(_K,t-1);
						int _jjmat = _jmat[t+kk-1][leng*j+n* _NScen/pow(_K,t+NN-1)];
						//assert(y_var[ii][t+kk-2][_jjmat] != NULL);
						SCIP_CALL( SCIPaddCoefLinear(subscip, wt_con[n], y_var[ii][kk][_jjmat-j*pow(_K,kk)], -1) );
						//SCIPinfoMessage(subscip, NULL, "Y at this constraint Y_%d_%d_%d", ii, t+kk, _jjmat);
					}
				}
				assert (pz_var != NULL );
				assert (wt_con[n] != NULL);
				SCIP_CALL( SCIPaddCoefLinear(subscip, wt_con[n], pz_var, -1000) );
				SCIPinfoMessage(subscip, NULL, "z at this constraint Z" );
				SCIP_CALL( SCIPaddCons(subscip, wt_con[n]) );

			}

			//Add constraint 6c
			//char con_name[255];
			vector< vector < SCIP_CONS* > > sps_con(_I);
			for (int ii=0; ii<_I; ii++)
			{
				sps_con[ii].resize(pow(_K,NN));
				for (int n=0; n<pow(_K,NN); n++)
				{
					SCIP_CONS* s_con = NULL;
					SCIPsnprintf(con_name, 255, "sps%d_%d", ii , n);
					SCIP_CALL( SCIPcreateConsLinear(subscip, & s_con, con_name, 0, NULL, NULL,
							0,
							SCIPinfinity(subscip),
							true,
							false,
							true,
							true,
							true,
							false,
							true,
							false,
							false,
							false) );
					sps_con[ii][n] = s_con;
					SCIP_CALL( SCIPaddCoefLinear(subscip, sps_con[ii][n], x_var[ii], _vialsize[ii]) );
					for(int kk=0; kk < NN+1 ; kk++)
					{
						int leng = _NScen/pow(_K,t-1);
						int _jjmat = _jmat[t+kk-1][leng*j+n* _NScen/pow(_K,t+NN-1)];
						//assert (y_var[ii][t+kk-2][_jjmat] != NULL);
						SCIP_CALL( SCIPaddCoefLinear(subscip, sps_con[ii][n], y_var[ii][kk][_jjmat-j*pow(_K,kk)], -1) );
					}
					SCIP_CALL( SCIPaddCons(subscip,sps_con[ii][n]) );
				}
			}

			SCIP_CALL( SCIPsolve(subscip) );
			SCIPinfoMessage(subscip, NULL, "solution status", SCIPgetStatus (subscip) );
			SCIP_CALL( SCIPprintStatistics(subscip, NULL) );

			SCIP_CALL( SCIPprintBestSol(subscip, NULL, FALSE) );


			//Get objective value
			objval [t-2][j] = SCIPgetPrimalbound(subscip);
			SCIPinfoMessage(subscip, NULL, "obj value", objval[t-2][j]);
			SCIP_SOL* sol;
			sol = SCIPgetBestSol (subscip);
			//SCIPinfoMessage(subscip, NULL, "SOL value", SCIPgetSolVal(subscip, sol, x_var[1]));
			for (int i=0; i<_I ; i++)
			{
				int floorvalx = floor (SCIPgetSolVal(subscip, sol, x_var[i]));
				int ceilvalx = ceil(SCIPgetSolVal(subscip, sol, x_var[i]));
				//if (sol == NULL)
				//cout<<"error sol"<<endl;
				if (SCIPgetSolVal(subscip, sol, x_var[i])- floorvalx < ceilvalx - SCIPgetSolVal(subscip, sol, x_var[i]))
					x_newcolumn[i][t-2][j] = floorvalx;
				else
					x_newcolumn[i][t-2][j] = ceilvalx;

				for (int kk=0; kk < NN+1 ; kk++)
				{
					for(int jj=0; jj < pow(_K,kk) ; jj++)
					{
						//int leng = pow(_K,t+kk-1)/pow(_K,t-1);
						int floorval = floor (SCIPgetSolVal(subscip, sol, y_var[i][kk][jj]));
						int ceilval = ceil(SCIPgetSolVal(subscip, sol, y_var[i][kk][jj]));
						//if (sol == NULL)
						//cout<<"error sol"<<endl;
						if (SCIPgetSolVal(subscip, sol, y_var[i][kk][jj])- floorval < ceilval - SCIPgetSolVal(subscip, sol, y_var[i][kk][jj]))
							y_newcolumn[i][t-2][j][kk][jj] = floorval;
						else
							y_newcolumn[i][t-2][j][kk][jj] = ceilval;
					}
				}
			}//for i
			//if (sol == NULL)
			//cout<<"error sol"<<endl;
			if (SCIPgetSolVal(subscip, sol, pz_var) > 0.5)
				z_newcolumn[t-2][j]= 1;
			else
				z_newcolumn[t-2][j]= 0;

			SCIP_CALL( SCIPfree(&subscip) );

			BMScheckEmptyMemory();

		}//j

	}//t

	return SCIP_OKAY;
}
