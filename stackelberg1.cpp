#include"LB.h"
#include"evolutionbit.h"
#include"nash.h"

int main(){
	srand((unsigned)time(NULL));
	int Time = 3;
	int CASEnum= 30;	

	vector<double>CONSIDER;
	int CON_VALUE = 25;
	double c = 0;
	for(int i=0;i <= CON_VALUE;i++){
		CONSIDER.push_back(c);
		c += 0.2;
	}
	int LOOP  = 100;

	for(int i =0;i<Time;i++){

		int conN = CONSIDER.size();
		vector<double> smlu(conN,0) ;
		vector<double> sdelay(conN,0);

		double mlu = 0;
		double mludelay = 0;

		double delaymlu = 0;
		double delay = 0;

		double nashmlu = 0;
		double nashdelay = 0;

		vector<int> successCase (conN, 0) ;

		int sucCaseLB = 0,sucCaseOR = 0,sucCaseNash = 0;

		for(int casenum = 0; casenum < CASEnum; casenum++){

			genGraph(12,64,"inputFile//graph.txt");
			genGraphOR(12,5,10,"inputFile//graphOR.txt");
			CGraph *G = new CGraph("inputFile//graph.txt");
			CGraph *GOR = new CGraph("inputFile//graphOR.txt");

			vector<demand> eqOR;
			vector<demand> eqTE;
			vector<demand> eqbase;

			//eqbase.clear();//background流 
			for(int i = 0; i < BGNUM; i++){
				int s = rand()%G->n, t;
				do{
					t = rand()%G->n;
				}while( s == t || G->canNotReach(s,t));
				eqbase.push_back(demand(s, t, rand()%(MAXDEMAND)+1));
			}

			////Overlay  产生demand
			//eqOR.clear(); 
			for(int i =0 ;i<GOR->m;i++)
				if(G->canNotReach(GOR->Link[i]->tail, GOR->Link[i]->head))
					continue;
				else
					eqOR.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,rand()%(MAXDEMAND)+1));

			//eqTE.clear(); 
			for(unsigned int i=0;i<eqOR.size();i++)
				eqTE.push_back(eqOR[i]);
			int ornum = eqTE.size();
			for(unsigned int i=0;i<eqbase.size();i++)
				eqTE.push_back(eqbase[i]);

			FILE *cur = fopen("outputFile//mludelay.csv", "a");
			fprintf(cur,"\n\n\n%d\n",casenum);
			
			double lbdic = 0,ordic=0;
			G->clearOcc();
			LBdictor(G,eqTE,ornum);

			if( G->mlu + 1e-5 <= INF ){
				sucCaseLB++;
				mlu += G->mlu;
				mludelay += G->delay;
				lbdic = G->mlu;
			}
			else
				break;

			fprintf(cur,",LB,,%lf,%lf\n",G->mlu,G->delay);

			G->clearOcc();
			ORdictor(G,eqbase,eqOR);

			if( G->delay + 1e-5 <= INF ){
				sucCaseOR++;
				delay += G->delay;
				delaymlu += G->mlu;
				ordic = G->delay;
			}
			else
				break;
			fprintf(cur,",OR,,%lf,%lf\n",G->mlu,G->delay);

			G->clearOcc();
			if(!G->GAinit(eqTE)){
				cout << "*****GA init failed***** "<<endl;
				break;
			}

			//// nash		
			FILE *nash = fopen("outputFile//nash.csv", "a");
			int nacase = 0;
			double loopnashmlu=0,loopnashdelay=0;
			fprintf(nash,"\n\n nash \n");
			for(int i =0;i<LOOP;i++){
				G->clearOcc();
				GOR->clearOcc();
				double curmlu = NashLB(G,GOR,eqTE);
				if(curmlu + 1e-5 >= INF){
					fprintf(nash,"NashEE unfeasible\n");
					break;
				}
				double curdelay = NashOR(GOR,eqOR);
				if( curdelay + 1e-5 >= INF){
					fprintf(nash,"NashOR unfeasible\n");
					break;
				}
				eqTE.clear();
				for(int i=0;i<GOR->m;i++){
					if(GOR->Link[i]->use>0)
						eqTE.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,GOR->Link[i]->use));
				}
				for(unsigned int i=0;i<eqbase.size();i++)
					eqTE.push_back(eqbase[i]);

				nacase++;
				loopnashmlu += curmlu;
				loopnashdelay += curdelay;
				fprintf(nash,"%lf,%lf\n",curmlu,curdelay);
				cout<<curmlu<<"\t"<<curdelay<<endl;
			}
			fclose(nash);

			if(nacase){
				nashmlu += (loopnashmlu/nacase);
				nashdelay += (loopnashdelay/nacase);
				sucCaseNash++;
			}	
			fprintf(cur,",Nash,,%lf,%lf\n\n",loopnashmlu/nacase,loopnashdelay/nacase);
			fclose(cur);

			for(unsigned int con = 0;con < CONSIDER.size();con++){

				int n = 150;//种群个体数目
				int m = eqTE.size();
				evoluPopubit popubit(n,m,G,GOR,&eqTE,&eqOR,lbdic,ordic,CONSIDER[con]);
				(popubit).evolution();
				cout<<"S\t"<<popubit.hero.mlu<<"\t"<<popubit.hero.delay <<endl;

				if( (popubit.hero.mlu +1e-5) >= INF ||  (popubit.hero.delay  + 1e-5) >= INF ){
					break;
				}

				cur = fopen("outputFile//mludelay.csv", "a");

				fprintf(cur,",S,%lf,%lf,%lf\n",CONSIDER[con],popubit.hero.mlu,popubit.hero.delay);
				fclose(cur);

				successCase[con] += 1;
				smlu[con] += popubit.hero.mlu;
				sdelay[con] += popubit.hero.delay;

			} // end of CONSIDER for

			delete G;
			delete GOR;

		} // end of CASENum for

		FILE *res = fopen("outputFile//result.csv", "a");		
		fprintf(res,"%d case average\n",CASEnum);
		fprintf(res,",,CONSIDER,successCase,Energy Efficiency,throughput\n");
		for(unsigned int con = 0;con < CONSIDER.size();con++){
			fprintf(res, ",S,%lf,%d,%lf,%lf\n",CONSIDER[con],successCase[con],smlu[con]/successCase[con],sdelay[con]/successCase[con]); 
		}
		fprintf(res,",,,successCase,Energy Efficiency,throughput\n");
		fprintf(res, "\n\n,EE,,%d,%lf,%lf\n",sucCaseLB,mlu/sucCaseLB,mludelay/sucCaseLB);
		fprintf(res, ",OR,,%d,%lf,%lf\n",sucCaseOR,delay/sucCaseOR,delaymlu/sucCaseOR);
		fprintf(res, ",Nash,,%d,%lf,%lf\n",sucCaseNash,nashmlu/sucCaseNash,nashdelay/sucCaseNash);
		fclose(res);
	}
	system("pause");
	return 0;
}