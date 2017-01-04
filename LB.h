#ifndef LB_H
#define LB_H
#include"CGraph.h"
#include <ilcplex/ilocplex.h>

double ilpSolver(CGraph *G,vector<demand> & req,int ornum){
	IloEnv env;
	IloModel model(env);
	IloCplex EEsolver(model);

	int num = req.size();	
	IloArray<IloIntVarArray> x(env, num); 
	for(int d = 0; d < num; d++)
		x[d] = IloIntVarArray(env, G->m, 0, 1); 	

	//优化目标
	IloNumVar z(env, 0, IloInfinity);
	model.add(IloMinimize(env, z));

	// 对每个点进行流量守恒约束  
	for(int d = 0; d < num; d++){
		for(int i = 0; i < G->n; i++){    // n为顶点数
			IloExpr constraint(env);
			for(unsigned int k = 0; k < G->adjL[i].size(); k++) // 出度边
				constraint += x[d][G->adjL[i][k]->id];
			for(unsigned int k = 0; k < G->adjRL[i].size(); k++) // 入度边
				constraint -= x[d][G->adjRL[i][k]->id];

			if(i == req[d].org)
				model.add(constraint == 1);
			else if(i == req[d].des)
				model.add(constraint == -1);
			else
				model.add(constraint == 0);
		}
	}

	for(int i = 0; i < G->m; i++){
		IloExpr load(env);
		for(int d = 0; d <  num; d++)
			load += req[d].flow*x[d][i];
		model.add( G->Link[i]->capacity > load);  
		model.add( z >= load);	
	}	

	EEsolver.setOut(env.getNullStream());
	double obj = INF;
	G->mlu = INF;
	if(EEsolver.solve()){
	    obj = EEsolver.getObjValue(); 
		G->mlu = EEsolver.getObjValue(); //mlu
		for(int i = 0;i < G->m;i++){  
			double loadc = 0;
			for(int d=0;d < num;d++)
				loadc += EEsolver.getValue(x[d][i])*req[d].flow;
			G->Link[i]->latency = linearCal(loadc,G->Link[i]->capacity); //拟合
		}
		
		double del = 0;
		for(int d = 0; d < ornum; d++){  
			for(int i = 0;i < G->m;i++)
				del +=  req[d].flow*EEsolver.getValue(x[d][i])*G->Link[i]->latency;
		}
		G->delay = del;
	}
	else{
		cout << "ilpsolver unfeasible"<<endl;
	}
	for(int i = 0; i < req.size(); i++)
		x[i].end();
	x.end();
	env.end();
	return obj;
}

#endif