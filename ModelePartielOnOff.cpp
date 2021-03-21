#include "ModelePartielOnOff.h"

#include<algorithm>
#include <numeric>
#include "utils.h"

const char* ModelePartielOnOff::s_params = "params.mod";

IloOplRunConfiguration ModelePartielOnOff::loadRC(IloEnv& env, const char* datfile)
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

void ModelePartielOnOff::initVars(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());

	
	// non proc
	for (IloInt i : _nonProcessed) 
	{
		_x.emplace(i, std::map<IloInt, IloNumVar>()); 
		_y.emplace(i, std::map<IloInt, IloNumVar>());
		_u.emplace(i, std::map<IloInt, IloNumVar>());

		std::ostringstream ossC;
		ossC << "C#" << i;
		std::string reprC(ossC.str());
		_C.emplace(i,IloIntVar(env, reprC.c_str()));

		std::ostringstream ossR;
		ossR << "Tardiness#" << i;
		std::string reprR(ossR.str());
		_T.emplace(i,IloIntVar(env, reprR.c_str()));

		std::ostringstream oss2;
		oss2 << "a#" << i;
		std::string repr2(oss2.str());
		_A.emplace(i,IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr2.c_str()));

		for (IloInt t(_a); t <= _b; ++t)
		{
			std::ostringstream oss;
			oss << "x#" << i << "#" << t;
			std::string repr(oss.str());
			_x[i].emplace(t,IloNumVar(env, 0, 1, ILOBOOL, repr.c_str()));
			
			std::ostringstream oss1;
			oss1 << "y#" << i << "#" << t;
			std::string repr1(oss1.str());
			_y[i].emplace(t, IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr1.c_str()));
		}

		for (IloInt j : _nonProcessed)
		{
			std::ostringstream oss;
			oss << "u#" << i << "#" << j;
			std::string repr(oss.str());
			_u[i].emplace(j,IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str()));
		}

	}

	Orders E(_nonProcessed);
	E.push_back(_precOrder);
	// add alpha and omega variables
	for (IloInt i : E)
	{
		std::ostringstream oss3;
		oss3 << "alpha#" << i;
		std::string repr3(oss3.str());
		_alpha.emplace(i, IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr3.c_str()));

		std::ostringstream oss4;
		oss4 << "omega#" << i;
		std::string repr4(oss4.str());
		_omega.emplace(i, IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr4.c_str()));


	}


	// prec Order
	_x.emplace(_precOrder, std::map<IloInt, IloNumVar>());
	for (IloInt t(_a); t<= _b ; ++t)
	{
		std::ostringstream ossPrec;
		ossPrec << "x#" << _precOrder << "#" << t;
		std::string reprPrec(ossPrec.str());

		_x[_precOrder].emplace(t, IloNumVar(env, 0, 1, ILOBOOL, reprPrec.c_str()));
	}
	for (IloInt t : _tt)
	{
		std::ostringstream ossPrec;
		ossPrec << "x#" << _precOrder << "#" << t;
		std::string reprPrec(ossPrec.str());

		if (_x.at(_precOrder).find(t) == _x.at(_precOrder).end())
		{
			_x[_precOrder].emplace(t, IloNumVar(env, 1, 1, ILOBOOL, reprPrec.c_str()));
		}
	}

	
	// a 
	std::ostringstream oss2;
	oss2 << "a#" << _precOrder;
	std::string repr2(oss2.str());
	_A.emplace(_precOrder, IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr2.c_str()));

	// T & C
	std::ostringstream ossC;
	ossC << "C#" << _precOrder;
	std::string reprC(ossC.str());
	_C.emplace(_precOrder, IloIntVar(env, reprC.c_str()));

	std::ostringstream ossR;
	ossR << "Tardiness#" << _precOrder;
	std::string reprR(ossR.str());
	_T.emplace(_precOrder, IloIntVar(env, reprR.c_str()));

	// u[_precOrder][j]
	for (IloInt j : _nonProcessed)
	{
		std::ostringstream oss;
		oss << "u#" << _precOrder << "#" << j;
		std::string repr(oss.str());
		_u[_precOrder].emplace(j, IloNumVar(env, IloFalse, IloTrue, ILOBOOL, repr.c_str()));
	}



}

void ModelePartielOnOff::initObj(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloNumMap w(opl.getElement("w").asNumMap());
	IloIntMap e(opl.getElement("e").asIntMap());
	IloNumMap c(opl.getElement("c").asNumMap());

	IloExpr obj(env);
	for (IloInt j : _nonProcessed)
	{
		IloExpr energyCost(env);
		for (IloInt t(std::max(_a,r.get(j))); t <= std::min(_b,db.get(j)); ++t)
		{
			energyCost += (_x[j][t] + _y[j][t]) * c.getSub(j).get(t);
		}

		obj += _A[j] * e.get(j) - w.get(j) * _T[j] - energyCost;
		//(0.8)*_alpha[j]*e.get(j) + _omega[j] * e.get(j) - w.get(j) * _T[j] - energyCost;
	}


	_model.add(IloMaximize(env, obj));


}

void ModelePartielOnOff::initConstraints(IloEnv& env)
{
	IloOplModel opl(_dat.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap d(opl.getElement("d").asIntMap());


	std::vector<IloInt> E(_nonProcessed);
	E.push_back(_precOrder);

	// C0
	for (IloInt t(_a); t <= _b; ++t)
	{
		IloExpr sum(env);
		for (IloInt j : _nonProcessed)
		{
			sum += (_x[j][t] + _y[j][t]);
		}
		_model.add(sum <= 1);

	}

	// Completion times
	for (IloInt i : _nonProcessed)
	{
		_model.add(_C[i] <= db.get(i) * _omega[i]);
		for (IloInt t(std::max(_a,r.get(i))); t <= std::min(_b-1,db.get(i)); ++t)
		{
			_model.add(_C[i] >= (t + 1) * (_x[i][t] - _x[i][t + 1]));
		}
		_model.add(_C[i] >= IloInt(0));
	}

	// Retard
	for (IloInt i : _nonProcessed)
	{
		_model.add(_T[i] >= IloInt(0));

		_model.add(_T[i] >= _C[i] - d.get(i) * _omega[i]);
		_model.add(_T[i] >= IloInt(0));
		_model.add(_T[i] <= (db.get(i) - d.get(i)) * _omega[i]);
	}


	// succ
	for (IloInt i : E)
	{
		IloExpr sum(env);
		for (IloInt j : _nonProcessed)
		{
			if (j != i)
			{
				sum += _u[i][j];
			}
		}

		_model.add(sum <= _omega[i]);
	}


	// pred
	for (IloInt i : _nonProcessed)
	{
		IloExpr sum(env);
		for (IloInt j:E)
		{
			if (j != i)
			{
				sum += _u[j][i];
			}
		}

		_model.add(sum == _omega[i]);
	}

	_model.add(_A[_precOrder] == 1);

	// C3
	for (IloInt j : _nonProcessed)
	{
		IloInt ub(std::min(_b,db.get(j) - p.get(j)));
		IloInt lb(std::max(_a,r.get(j)));
		//IloInt ub(db.get(j));
		for (IloInt t(lb); t <= ub; ++t)
		{

			IloExpr term1(env);

			for (IloInt tp(std::max(_a,r.get(j))); tp <= t - p.get(j); ++tp)
			{
				term1 += _x[j][tp];

			}

			IloExpr term2(env);

			for (IloInt tp(std::min(_b,t + p.get(j))); tp <= std::min(_b,db.get(j)); ++tp)
			{
				term2 += _x[j][tp];

			}

			_model.add(term1 + term2 <= p.get(j) * (1 - _x[j][t]) );

		}


	}

	//	C4
	for (IloInt j : _nonProcessed)
	{
		IloExpr sum(env);
		for (IloInt t(std::max(_a,r.get(j))); t <= std::min(_b,db.get(j)); ++t)
		{
			sum += _x[j][t];
		}

		_model.add(sum == p.get(j) * _A[j]);
	}

	for (IloInt i : E)
	{
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{
				IloInt lb(std::max(_a+1, r.get(j)+1));
				IloInt ub(std::min(_b,db.get(j) - p.get(j)));
				for (IloInt t(lb); t <= ub; ++t)
				{

					IloExpr sumY(env);
					for (IloInt tp(std::min(_a,r.get(j))); tp <= t - 1; ++tp)
					{
						if (tp >= _a)
						{
							sumY += _y[j][tp];
						}
					}


					_model.add(sumY >=
						s.getSub(i).get(j) * (_u[i][j] + _x[j][t] - 1)
						

					);


				}
			}
		}
	}

	// C6
	for (IloInt j : _nonProcessed)
	{
		for (IloInt t(std::max(_a, r.get(j))); t <= std::min(_b-1,db.get(j)); ++t)
		{
			_model.add(_y[j][t] - _y[j][t + 1] - _x[j][t + 1] <= IloInt(0));
		}
	}

	for (IloInt j : _nonProcessed)
	{
		IloExpr sum1(env);
		IloInt lb(std::max(_a, r.get(j)));
		IloInt ub(std::min(_b,db.get(j) - p.get(j)));
		for (IloInt t(lb); t <= ub; ++t)
		{
			sum1 += _y[j][t];
		}

		IloExpr sum2(env);
		for (IloInt i : E)
		{
			if (i != j)
			{
				sum2 += (_u[i][j]) * s.getSub(i).get(j);
			}
		}

		_model.add(sum1 <= sum2);
	}


	 // C8
	for (IloInt i : E)
	{
		for (IloInt j : _nonProcessed)
		{
			if (i != j)
			{

				IloInt lb(std::max(_a, r.get(j)));
				IloInt ub(std::min(_b, db.get(j)));

				for (IloInt t(lb); t <= ub; ++t)
				{
					IloExpr sum(env);
					IloInt lb1(std::min(_b,std::max(_a, r.get(i) + 1)));
					
					for (IloInt tp(lb1); tp <= t - 1; ++tp)
					{  
						sum += _x[i][tp];
						
					}
					_model.add(sum >= p.get(i) * (_u[i][j] + _y[j][t] - 1)

					);

				}
			}
		}
	}

	IloConstraintArray C11(env, n), C11bis(env, n), C12(env, n), C12bis(env, n);
	for (IloInt j : _nonProcessed)
	{
		IloExpr sum1(env), sum2(env);
		for (IloInt t(std::max(_a,IloInt(0))); t < std::min(_b,r.get(j)); ++t)
		{
			sum1 += _x[j][t];
			sum2 += _y[j][t];
		}
		_model.add(sum1 == 0);
		_model.add(sum2 == 0);

		IloExpr sum3(env), sum4(env);
		for (IloInt t(std::min(_b,db.get(j) + 1)); t <= std::min(_b,T); ++t)
		{
			if (t >= _a)
			{
				sum3 += _x[j][t];
			}
		}
		for (IloInt t(std::min(_b,db.get(j)-p.get(j)+1)); t <= std::min(_b,T); ++t)
		{
			if (t >= _a)
			{
				sum4 += _y[j][t];
			}
		}
		_model.add(sum3 == 0);
		_model.add(sum4 == 0);
	}

	for (IloInt t : _tt)
	{
		std::cout << "x" << _precOrder << t << std::endl;
		_model.add(_x[_precOrder][t] == 1);
	}


	/*IloExpr sumZz(env);
	for (IloInt t(_a); t <= _b; ++t)
	{
		sumZz += _x[_precOrder][t];
	}
	IloConstraint nonZ(sumZz == 0);
	nonZ.setName("nonZ");
	opl.getModel().add(nonZ); */



	IloExpr sumApprox(env);
	for (IloInt j : _nonProcessed)
	{

		IloInt setup(utils::computeSetup(_setup, j, E, s));

		sumApprox += (p.get(j) + setup) * _alpha[j];
	}

	IloConstraint approx(sumApprox <= T - _b);
	approx.setName("approx");
	_model.add(approx);

	std::cout << "T-b=" << T - _b << std::endl;


	// define omega
	for (IloInt j : _nonProcessed)
	{
		IloExpr subsum2(env);

		for (IloInt t(_a); t <= _b; ++t)
		{
			subsum2 += _x[j][t];
		}

		_model.add(_omega[j] >= p.get(j)*subsum2);
	}

	// link alpha omega
	for (IloInt j : _nonProcessed)
	{
		_model.add(_alpha[j] == 1 - _omega[j]);
		//_model.add(_omega[j] <= 1 - _alpha[j]);
	}

	_model.add(_alpha[_precOrder] == 0);
	_model.add(_omega[_precOrder] == 1);


	for (IloInt j : _nonProcessed)
	{
		if (_b >= db.get(j))
		{
			_model.add(_alpha[j] == 0);
		}
	}

	
}

ModelePartielOnOff* ModelePartielOnOff::load(IloEnv& env, IloOplRunConfiguration& rc, SETUP setup, const IloInt& a, const IloInt& b, const Orders& nonProcessed, const IloInt& precOrder, const std::vector<IloInt>& tt)
{
	return new ModelePartielOnOff(env, rc, setup, a, b, nonProcessed, precOrder, tt);
}

void ModelePartielOnOff::relaxAndFix(IloEnv& env, const char* datfile, const IloInt& sigma, const IloInt& delta, SETUP setup)
{
	IloOplRunConfiguration rc(loadRC(env, datfile));

	IloInt T(rc.getOplModel().getElement("T").asInt());
	IloInt n(rc.getOplModel().getElement("n").asInt());


	IloInt b(sigma - 1);
	IloInt a(0);
	IloInt k(0);

	Orders nonProc;

	IloInt precOrder(0);
	std::vector<IloInt> tt;
	tt.push_back(0);

	IntMatrix xvals(env, n + 1);
	for (IloInt i(0); i <= n; ++i)
	{
		xvals[i] = IloIntArray(env, T + 1);
		for (IloInt t(0); t <= T; ++t)
		{
			xvals[i][t] = 0;
		}

		if (i != 0)
		{
			nonProc.push_back(i);
		}
	}




	double timeBudg(3600);

	while (b < T)
	{

		relaxAndFixLoop(rc, setup, k, a, b, delta, precOrder, nonProc, tt, xvals);

		a = b - delta;
		b = b + sigma - delta;
		++k;

		if (b > T)
		{
			b = T;
		}
	}


	std::cout << "last loop" << std::endl;
	relaxAndFixLoop(rc, setup, k, a, b, delta, precOrder, nonProc, tt, xvals);



	Orders seq;
	std::map<IloInt, IloInt> u;
	// get sequence dirty
	IloInt previous(0);
	for (IloInt t(0); t <= T; ++t)
	{
		for (IloInt i(1); i <= n; ++i)
		{
			if (xvals[i][t])
			{
				u.emplace(previous, i);
				previous = i;
				seq.push_back(i);
			}
		}
	}

	IloOplModel opl(rc.getOplModel());
	IloIntMap p(opl.getElement("p").asIntMap());
	IloIntMap s(opl.getElement("s").asIntMap());
	for (IloInt i(0); i <= n; ++i)
	{
		for (IloInt t(0); t <= T; ++t)
		{
			if (xvals[i][t] == 1)
			{
				if (u.find(i) != u.end())
				{
					IloInt prev(u.at(i));
					//std::cout << "ST[" << i << "]==" << t - s.getSub(prev).get(i) -1 << ";" << std::endl;
					std::cout << "C[" << i << "]==" << t + p.get(i) - 1 << ";" << std::endl;
				}
			}
		}

	}



	for (const std::pair<IloInt, IloInt>& p : u)
	{
		std::cout << "u[" << p.first << "][" << p.second << "]==" << 1 << ";" << std::endl;
	}
	std::cout << std::endl;

	//std::cout << "obj=" << computeObj(rc, xvals, seq) << std::endl;

}

void ModelePartielOnOff::relaxAndFixLoop(IloOplRunConfiguration& rc, SETUP setup, const IloInt& k, const IloInt& a, const IloInt& b, const IloInt& delta, IloInt& precOrder, Orders& nonProc, std::vector< IloInt> & tt, IntMatrix& vals)
{
	std::cout << "non processed orders: " << std::endl;
	for (IloInt i : nonProc)
	{
		std::cout << i << " ";
	}
	std::cout << std::endl;

	std::cout << "precOrder=" << precOrder << std::endl;

	IloEnv env1;
	std::cout << "a=" << a << "b=" << b << std::endl;
	ModelePartielOnOff* subProblem = load(env1, rc, setup, a, b, nonProc, precOrder, tt);


	if (k > 0)
	{
		IloInt to(std::max(IloInt(0), a - 1));

		size_t len(tt.size());
		IloInt from(*tt.begin());
		std::cout << "fix from " << from << " to " << to << std::endl;
		subProblem->fix(env1, vals, from, to);
	}


	IloCplex cplx(env1);
	cplx.extract(subProblem->_model.getImpl());

	cplx.setOut(env1.getNullStream());
	cplx.setParam(IloCplex::Param::Threads, 4);
	cplx.setParam(IloCplex::Param::TimeLimit, 15);

	/*std::ostringstream oss;
	oss << "subModel_" << a << "_" << b << "_pulse" << ".LP";
	std::string fname(oss.str());
	cplx.exportModel(fname.c_str());*/

	if (cplx.solve())
	{

		std::cout << cplx.getStatus() << std::endl;

		IloInt to(b - delta - 1);
		std::cout << "get from " << a << " to " << to << std::endl;
		subProblem->get(cplx, vals, a, to, nonProc, precOrder, tt);

		subProblem->printCplx(std::cout, cplx);

		for (IloInt i : nonProc)
		{
			std::cout << "omega#" << i << "=" << cplx.getValue(subProblem->_omega[i]) << std::endl;
		}

		for (IloInt i : nonProc)
		{
			std::cout << "alpha#" << i << "=" << cplx.getValue(subProblem->_alpha[i]) << std::endl;
		}

		std::cout << "partial obj=" << cplx.getObjValue() << std::endl;
		//std::cout << "tt=" << tt << std::endl;

	}
	else
	{
		std::cerr << cplx.getStatus() << std::endl;


		//printConflict(env, cplx, subProblem);

		return;
	}


	delete subProblem;
	cplx.end();
	env1.end();
}

void ModelePartielOnOff::fix(IloEnv& env, const IntMatrix& vals, const IloInt& from, const IloInt& to)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());


	for (IloInt t : _tt)
	{
		_model.add(_x[_precOrder][t] == vals[_precOrder][t]);
	}

}

void ModelePartielOnOff::get(IloCplex& cplx, IntMatrix& vals, const IloInt& from, const IloInt& to, std::vector<IloInt>& orders, IloInt& precOrder, std::vector<IloInt>& tt)
{
	IloOplModel opl(_dat.getOplModel());
	IloInt n(opl.getElement("n").asInt());
	IloIntMap p(opl.getElement("p").asIntMap());

	Orders E(orders);

	IloInt prevOrder(_precOrder);

	//E.push_back(0);

	std::map<IloInt, std::vector<IloInt> > O;
	for (IloInt t(from); t <= to; ++t)
	{
		for (IloInt i : E)
		{
			IloNumVar x(_x[i][t]);
			
			vals[i][t] = IloRound(cplx.getValue(_x[i][t]));


			if (vals[i][t] == 1 and i != 0)
			{
				if (O.find(i) == O.end())
				{
					O.emplace(i, std::vector<IloInt>());
				}
				
			}


			if (O.find(i) != O.end())
			{
				IloInt sum(std::accumulate(O.at(i).begin(), O.at(i).end(), 0));

				if (sum >= p.get(i))
				{
					orders.erase(std::find(orders.begin(), orders.end(), i));
					O[i].push_back(t);
					precOrder = i;
					tt.push_back(t);
				}
			}

			if (vals[i][t] and i == 0)
			{
				if (O.find(i) == O.end())
				{
					O.emplace(i, std::vector<IloInt>());
				}
				O[i].push_back(0);
				precOrder = i;
				
			}
		}
	}

	if (precOrder == 0)
	{
		vals[0][0] = 1;
	}



}

double ModelePartielOnOff::computeObj(IloOplRunConfiguration& rc, const IntMatrix& xvals, const IntMatrix & yvals, const IntArray& A, const IntArray& Tardiness, const Orders& sequence)
{

	IloOplModel opl(rc.getOplModel());

	IloInt n(opl.getElement("n").asInt());
	IloInt T(opl.getElement("T").asInt());
	IloIntMap db(opl.getElement("db").asIntMap());
	IloIntMap r(opl.getElement("r").asIntMap());
	IloNumMap w(opl.getElement("w").asNumMap());
	IloIntMap e(opl.getElement("e").asIntMap());
	IloNumMap c(opl.getElement("c").asNumMap());

	IloNum obj(0.0);
	for (IloInt j(1); j <= n; ++j)
	{
		IloNum energyCost(0.0);
		for (IloInt t(r.get(j)); t <= db.get(j); ++t)
		{
			energyCost += (xvals[j][t] + yvals[j][t]) * c.getSub(j).get(t);
		}

		obj += A[j]*e.get(j) - w.get(j) * Tardiness[j] - energyCost;
	}

	return obj;

}

void ModelePartielOnOff::printCplx(std::ostream& os, IloCplex & cplx)
{
	for (IloInt i : _nonProcessed)
	{
		for (IloInt t(_a); t <= _b; ++t)
		{
			IloNumVar x(_x[i][t]);
			IloInt val(cplx.getValue(x));
			std::cout << val;
		}
		std::cout << std::endl;
	}

}
