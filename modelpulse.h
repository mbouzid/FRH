#pragma once
#include <ilopl/iloopl.h>


typedef IloArray<IloNumVarArray>  NumVarMatrix;
typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloIntArray> IntMatrix;

class ModelPulse
{

	static const char* s_params;

	private:

		IloModel _model;
		IloOplRunConfiguration _dat;
		
		
		NumVarMatrix _x;
		NumVarMatrix _u;

		IloNumVarArray _alpha;
		IloNumVarArray _omega;

		IntMatrix _xvalues;

		static IloOplRunConfiguration loadRC(IloEnv& env, const char* datfile);

		void initVars(IloEnv & env);
		void initObj(IloEnv& env);

		void initConstraints(IloEnv& env);


	public:

		ModelPulse
		(
			IloEnv& env,
			const IloOplRunConfiguration& dat,
			const IntMatrix& values
		):
			_model(env),
			_dat(dat),
			_x(env),
			_u(env),
			_alpha(env),
			_omega(env),
			_xvalues(values)

		{
			initVars(env);
			initObj(env);
			initConstraints(env);
		}

		static ModelPulse* load(IloEnv& env, const char* datfile);

		void prepareRelaxedModel(IloEnv& env, IloModel & model, const IloInt& a, const IloInt& b);


		void convertToBool(IloEnv & env, IloExtractableArray & conversions, const IloInt& from, const IloInt& to);

		void fix(IloModel & model, const IloInt& from, const IloInt& to);

		void addMIPStart(IloEnv & env, IloCplex & cplx, const IloInt& from, const IloInt& to);

		void relaxAndFix(IloEnv & env, const IloInt & sigma, const IloInt & d);

		void get(IloCplex& cplx, const IloInt& from, const IloInt& to);
		
		void printVars(IloCplex& cplx);


};

void printConflict(IloEnv & env, IloCplex& cplx, const IloModel& model);

