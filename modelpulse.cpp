#include "modelpulse.h"
#include <algorithm>
#include <vector>

const char* ModelPulse::s_params = "params.mod";

IloOplRunConfiguration ModelPulse::loadRC(IloEnv& env, const char* datfile)
{
	IloOplErrorHandler errHdlr(env, std::cerr);
	IloOplSettings settings(env, errHdlr);
	IloOplModelSource src(env, s_params);

	IloOplModelDefinition def(src, settings);
	IloOplDataSource dat(env, datfile);

	IloOplDataElements elts(def, dat);
	IloOplRunConfiguration rc(def, elts);
	rc.getOplModel().generate();

	return rc;
}

void ModelPulse::initVars(IloEnv & env)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());

	_x = NumVarMatrix(env, n + 1);
	_u = NumVarMatrix(env, n + 1);
	_alpha = IloNumVarArray(env, n + 1);
	_omega = IloNumVarArray(env, n + 1);


	// x_it
	for (IloInt i(0); i <= n; ++i)
	{
		_x[i] = IloNumVarArray(env, T + 1);

		_u[i] = IloNumVarArray(env, n + 1);



		for (IloInt t(0); t <= T; ++t)
		{
			std::ostringstream oss;
			oss << "x#" << i << "#" << t;
			std::string repr(oss.str());

			_x[i][t] = IloNumVar(env, 0, 1, ILOFLOAT, repr.c_str());
			
		}

		for (IloInt j(0); j <= n; ++j)
		{
			std::ostringstream oss;
			oss << "u#" << i << "#" << j;
			std::string repr(oss.str());
			_u[i][j] = IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str());
		}

	}

	// add alpha variables
	for (IloInt j(0); j <= n; ++j)
	{
		std::ostringstream oss;
		oss << "alpha#" << j;
		std::string repr(oss.str());
		_alpha[j] = IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str());

		std::ostringstream oss1;
		oss1 << "omega#" << j;
		std::string repr1(oss1.str());
		_omega[j] = IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr1.c_str());
	}


}

void ModelPulse::initObj(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloNumMap f(opl.getElement("f").asNumMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	IloNumMap c(opl.getElement("c").asNumMap());
	IloIntMap e(opl.getElement("e").asIntMap());



	IloExpr obj(env);
	for (IloInt j(1); j <= n; ++j)
	{
		IloExpr profit(env);
		for (IloInt t(r.get(j)); t <= db.get(j) - p.get(j) + 1; ++t)
		{
			IloNum energyCostProcessing(0.0);
			for (IloInt tp(0); tp < p.get(j); ++tp)
			{
				energyCostProcessing += c.getSub(j).get(t + tp);

			}

			IloExpr	energyCostSetup(env);
			for (IloInt i(0); i <= n; ++i)
			{
				if (i != j)
				{
					IloNum cost(0.0);
					for (IloInt tp(1); tp <= s.getSub(i).get(j); ++tp)
					{
						if (t - tp >= r.get(j))
						{
							cost += c.getSub(j).get(t - tp);
						}

					}

					energyCostSetup += _u[i][j] * cost;
				}
			}


			profit += _x[j][t] * (f.getSub(j).get(t + p.get(j) - 1) - energyCostProcessing - energyCostSetup );

		}

		obj += profit;

	}


	_model.add(IloMaximize(env, obj));

}

void ModelPulse::initConstraints(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap d(opl.getElement("d").asIntMap());



	// capacity	 

	for (IloInt t(0); t <= T; ++t)
	{
		IloExpr capacity(env);
		for (IloInt j(1); j <= n; ++j)
		{
			capacity += _x[j][t];
		}
		_model.add(capacity <= 1);
	}

	// capacity 1
	for (IloInt j(1); j <= n; ++j)
	{
		IloExpr sum(env);
		for (IloInt t(r.get(j)); t <= db.get(j) - p.get(j) + 1; ++t)
		{
			sum += _x[j][t];
		}
		_model.add(sum <= 1);
	}


	// def1	& def2
	for (IloInt j(1); j <= n; ++j)
	{
		IloExpr sum1(env);
		for (IloInt t(0); t < r.get(j); ++t)
		{
			sum1 += _x[j][t];
		}
		_model.add(sum1 == 0);

		IloExpr sum2(env);
		for (IloInt t((db.get(j) - p.get(j) + 1) + 1); t <= T; ++t)
		{
			sum1 += _x[j][t];
		}
		_model.add(sum1 == 0);

	}

	_model.add(_x[0][0] == 1);

	IloExpr sumZ(env);

	for (IloInt t(1); t <= T; ++t)
	{
		sumZ += _x[0][t];
	}

	_model.add(sumZ == 0);



}

ModelPulse* ModelPulse::load(IloEnv& env, const char* datfile)
{

	IloOplRunConfiguration rc(loadRC(env, datfile));


	IloInt T(rc.getOplModel().getElement("T").asInt());
	IloInt n(rc.getOplModel().getElement("n").asInt());

	IntMatrix values(env, n + 1);
	for (IloInt i(0); i <= n; ++i)
	{
		values[i] = IloIntArray(env, T + 1);
		for (IloInt t(0); t <= T; ++t)
		{
			values[i][t] = 0;
		}
	}




	return new ModelPulse(env, rc,values);


}

void ModelPulse::prepareRelaxedModel(IloEnv & env, IloModel & model, const IloInt & a, const IloInt & b)
{
	IloOplModel opl(_dat.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap db(opl.getElement("db").asIntMap());


	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt j(1); j <= n; ++j)
		{
			if (i != j)
			{
				IloExpr sumA(env);
				IloInt ub1(/*std::min(b - p.get(i) + 1,*/ db.get(i) - p.get(i) + 1/*)*/);
				for (IloInt t(r.get(i)); t <= ub1; ++t)
				{
					sumA += (t)*_x[i][t];
				}
				sumA += ((s.getSub(i).get(j) + p.get(i)) * _u[i][j]) - (db.get(i)) * (1 - _u[i][j])
					- (db.get(i)) * (1 - _omega[j]) - db.get(i) * (1 - _omega[i]) + db.get(i) * (_omega[i] - _omega[j])
					
					;

				if (i == 0)
					sumA += _x[0][0];

				IloInt lb(r.get(j) + s.getSub(i).get(j) + 1);

				
				IloInt ub(/*std::min(*/db.get(j) - p.get(j) + 1/*, b - p.get(j) + 1)*/);
				
				IloExpr sumB(env);
				for (IloInt t(lb); t <= ub; ++t)
				{
					sumB += t * _x[j][t];
				}

				IloConstraint prec(sumA <= sumB);
				std::ostringstream oss;
				oss << "prec#" << i << "#" << j;
				std::string name(oss.str());

				prec.setName(name.c_str());

				model.add(prec);


			}
		}
	}


	for (IloInt i(0); i <= n; ++i)
	{
		IloExpr sum1(env);
		for (IloInt j(1); j <= n; ++j)
		{
			if (i != j)
			{
				sum1 += _u[i][j];
			}
		}

		IloExpr sum2(env);
		for (IloInt t(r.get(i)); t <= db.get(i) - p.get(i) + 1; ++t)
		{
			sum2 += _x[i][t];
		}

		model.add(sum1 <= _omega[i]);
	}


	// pred
	for (IloInt i(1); i <= n; ++i)
	{
		IloExpr sum1(env);
		for (IloInt j(0); j <= n; ++j)
		{
			if (i != j)
			{
					sum1 += _u[j][i];
			}
		}


		model.add(sum1 == _omega[i]);
	}

	// approx

	IloExpr sum(env);
	for (IloInt j(1); j <= n; ++j)
	{
		IloInt setup(10);
		IloExpr overlap(env);
		if (b + 1 <= T)
		{
			for (IloInt t(std::max(IloInt(1), (b + 1) - p.get(j) + 1)); t <= b + 1; ++t)
			{
				overlap += ((t + p.get(j)) - b) * (_x[j][t]);
			}
		}

		sum += overlap +  p.get(j) * _alpha[j] + setup * _alpha[j];
	}

	IloConstraint approx(sum <= T - b);
	approx.setName("approx");
	model.add(approx);

	std::cout << "T-b=" << T - b << std::endl;
	

	// define alpha
	for (IloInt j(1); j <= n; ++j)
	{
		IloExpr subsum1(env);

		for (IloInt t(b + 1); t <= T; ++t)
		{
			subsum1 += _x[j][t];
		}

		model.add(subsum1 == _alpha[j]);

	}



	// define omega
	for (IloInt j(0); j <= n; ++j)
	{
		IloExpr subsum2(env);

		for (IloInt t(0); t <= b; ++t)
		{
			subsum2 += _x[j][t];
		}

		model.add(_omega[j] >= subsum2);
	}

	// link alpha omega
	for (IloInt j(1); j <= n; ++j)
	{
		model.add(_alpha[j] <= 1 - _omega[j]);
		model.add(_omega[j] <= 1 - _alpha[j]);
	}

	model.add(_omega[0] == 1);

}

void ModelPulse::convertToBool(IloEnv & env, IloExtractableArray& conversions, const IloInt& from, const IloInt& to)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());

	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(from); t <= to; ++t)
		{
		
			conversions.add(IloConversion(env, _x[i][t], ILOBOOL));
		}
	}

}

void ModelPulse::fix(IloModel & model, const IloInt& from, const IloInt& to)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());

	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(from); t <= to; ++t)
		{
			if (_xvalues[i][t] == IloInt(1))
				std::cout << "x#" << i << "#" << t << std::endl;
			model.add( _x[i][t] == _xvalues[i][t]);
		}
	}
}

void ModelPulse::addMIPStart(IloEnv & env, IloCplex & cplx, const IloInt& from, const IloInt& to)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());

	IloNumVarArray vars(env);
	IloNumArray vals(env);

	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(from); t <= to; ++t)
		{
			if (_xvalues[i][t] == IloInt(1))
				std::cout << "x#" << i << "#" << t << std::endl;
			vars.add(_x[i][t]);
			vals.add(_xvalues[i][t]);
		}
	}


	cplx.addMIPStart(vars, vals);

}

void ModelPulse::relaxAndFix(IloEnv & env, const IloInt& sigma, const IloInt& delta)
{
	IloInt a(0), b(sigma - 1);
	IloOplModel opl(_dat.getOplModel());

	IloInt T(opl.getElement("T").asInt());
	std::cout << "T=" << T << std::endl;
	IloInt k(0);


	IloInt bprev(b);
	bool infeasible(false);

	while (b < T)
	{
		std::cout << "a=" << a << ",b=" << b << std::endl;
		IloModel subProblem(env);

		subProblem.add(_model.getImpl());

		if (k > 0)
		{
			IloInt to(std::max(IloInt(0), a -1 ));
			std::cout << "fix from " << 0 << " to " << to << std::endl;
			fix(subProblem, 0, to);
		}

		IloExtractableArray conversions(env);
		convertToBool(env, conversions,0, b);
		/*if (k == 0)
		{
			std::cout << "cast to boolean from " << a << " to " << b << std::endl;
			convertToBool(env, conversions,a, b);
		}
		else
		{
			std::cout << "cast to boolean from " << bprev+1 << " to " << b << std::endl;
			convertToBool(env, conversions,bprev+1, b);
		}*/

		subProblem.add(conversions);


		prepareRelaxedModel(env, subProblem, a, b);


		IloCplex cplx(env);
		cplx.extract(subProblem);
		
		/*std::ostringstream oss;
		oss << "subModel_" << a << "_" << b << "_pulse" << ".LP";
		std::string fname(oss.str());
		cplx.exportModel(fname.c_str()); */
		

		cplx.setOut(env.getNullStream());
		cplx.setParam(IloCplex::Param::Threads, 4);
		cplx.setParam(IloCplex::Param::TimeLimit, 10);

		if (cplx.solve())
		{


			std::cerr << cplx.getStatus() << std::endl;
			//printVars(cplx);

			IloInt to(b - delta -1);
			std::cout << "get from " << a << " to " << to << std::endl;
			get(cplx, a, to);

			IloInt n(opl.getElement("n").asInt());
			for (IloInt i(0); i <= n; ++i)
			{
				IloInt val(cplx.getValue(_omega[i]));
				if (val==1)
					std::cout << "omega#" << i  << "=" << val << std::endl;
			}

			for (IloInt i(0); i <= n; ++i)
			{
				for (IloInt j(1); j <= n; ++j)
				{
					if (cplx.isExtracted(_u[i][j]))
					{
						IloInt val(cplx.getValue(_u[i][j]));
						if (val==1)
							std::cout << "u#" << i << "#" << j << "=" << val << std::endl;
					}
				}
			}

			std::cout << "partial obj=" << cplx.getObjValue() << std::endl;

		}
		else
		{  
			std::cerr << cplx.getStatus() << std::endl;


			infeasible = true;
			
			//printConflict(env, cplx, subProblem);

			break;
		}

		a = b - delta;
		bprev = b;
		b = b + sigma - delta ;
		++k;

		if (b > T)
		{
			b = T;
		}

		conversions.endElements();

	}

	if (not infeasible)
	{
		IloInt n(opl.getElement("n").asInt());
		for (IloInt i(0); i <= n; ++i)
		{
			for (IloInt t(0); t <= T; ++t)
			{
				std::cout << _xvalues[i][t];
			}
			std::cout << std::endl;
		}
	}

}

void ModelPulse::get(IloCplex& cplx, const IloInt& from, const IloInt& to)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());

	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(from); t <= to; ++t)
		{
			IloNumVar x(_x[i][t]);
			
			_xvalues[i][t] = IloRound(cplx.getValue(_x[i][t]));
		}
	}
}

void ModelPulse::printVars(IloCplex& cplx)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());

	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(0); t <= T; ++t)
		{
			std::cout << cplx.getValue(_x[i][t]);
		}
		std::cout << std::endl;
	}
}

void printConflict(IloEnv & env, IloCplex& cplx, const IloModel& model)
{
	IloConstraintArray infeas(env);
	IloNumArray preferences(env);

	IloModel::Iterator it(model);

	while (it.ok())
	{
		if ((*it).isConstraint())
		{
			infeas.add((*it).asConstraint());
		}

		++it;
	}


	for (IloInt i = 0; i < infeas.getSize(); i++)
	{
		preferences.add(1.0);  // user may wish to assign unique preferences
	}

	if (cplx.refineConflict(infeas, preferences)) {
		IloCplex::ConflictStatusArray conflict = cplx.getConflict(infeas);
		env.getImpl()->useDetailedDisplay(IloTrue);
		std::cout << "Conflict :" << std::endl;
		for (IloInt i = 0; i < infeas.getSize(); i++) {
			if (conflict[i] == IloCplex::ConflictMember)
				std::cout << "Proved  : " << infeas[i] << std::endl;
			if (conflict[i] == IloCplex::ConflictPossibleMember)
				std::cout << "Possible: " << infeas[i] << std::endl;
		}
	}
	else
	{
		std::cout << "Conflict could not be refined" << std::endl;
	}
	std::cout << std::endl;

}
