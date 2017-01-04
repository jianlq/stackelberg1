#ifndef EVOLUTION
#define EVOLUTION
#include "Common.h"
#include"CGraph.h"

extern int ORNUM;


class evoluDiv{
	private:
		static const int MUT = 5;
		static const int CUL = 30;
		static const int HOR = 20;
		CGraph *G;
		CGraph *GOR;
		vector<demand> *dem;
	public:
		vector<double> x;
		double ability;
		double delay;
		evoluDiv() {;}
		evoluDiv(int m, CGraph *g, CGraph *gor, vector<demand> *d){
			x.resize(g->m);
			G = g;
			GOR = gor;
			dem = d;
			randNature();	
		}
		evoluDiv(vector<double> &tx, CGraph *g,  CGraph *gor,vector<demand> *d){
			x.clear();
			G = g;
			GOR = gor;
			dem = d;
			for(int i = 0; i < tx.size(); i++)
				x.push_back(tx[i]);
		}

		evoluDiv mate(evoluDiv other){ 
			vector<double> nx;
			for(int i = 0; i < x.size(); i++){
				double r = 1.0 * rand() / (RAND_MAX);
				nx.push_back(r * x[i] + (1 - r) * other.x[i]);
			}
			return evoluDiv(nx, G, GOR,dem);
		}

		void CalORDelay(){
			for(int i = 0; i < G->m; i++){
				G->Link[i]->latency = linearCal( G->Link[i]->use,G->Link[i]->capacity);
			}
			for(int i = 0; i < GOR->m; i++)
				GOR->Link[i]->latency = G->dijkstra(GOR->Link[i]->tail,GOR->Link[i]->head,0);
		}

		////ability
		void calAbility(){  
			ability = 0;
			delay = 0;
			for(int i = 0; i < G->m; i++)
				G->Link[i]->dist = x[i];

			G->clearOcc();
			double can = 0;
			for(int i = 0; i < (*dem).size(); i++){
				double dis = G->dijkstraLB(i,(*dem)[i].org, (*dem)[i].des, (*dem)[i].flow);
				can = max(can, dis);
			}

			double util1 = 0;
			for(int i = 0; i < G->m; i++)
				util1 = max(util1, (double)G->Link[i]->use);
			
			double del1 = 0;
			for(int d = 0; d < ORNUM; d++){
				for(unsigned int ij = 0; ij < G->reqPathID[d].size(); ij++)
					del1 +=  (*dem)[d].flow*linearCal(G->Link[G->reqPathID[d][ij]]->use,G->Link[G->reqPathID[d][ij]]->capacity);
			}

			GOR->clearOcc();
			CalORDelay();
			double del2 = 0;
			for(int i = 0; i < ORNUM; i++)
				del2 += (*dem)[i].flow * GOR->dijkstra((*dem)[i].org,(*dem)[i].des,(*dem)[i].flow);
			
			/*if(del1 > INF || del2 >INF){
				double del = 0;
				for(int i = 0; i < ORNUM; i++){
					del += (*dem)[i].flow * GOR->dijkstra((*dem)[i].org,(*dem)[i].des,(*dem)[i].flow);
					cout <<GOR->canNotReach((*dem)[i].org,(*dem)[i].des)<< "  && "<<del<<"    ";
				}
				cout << endl<<"INF"<<endl;
				for(int i=0;i<G->m;i++)
					cout << G->Link[i]->latency<<" G   "<<endl;
				cout <<endl;
				for(int i=0;i<GOR->m;i++)
					cout << GOR->Link[i]->latency<< "  GOR  "<<G->canNotReach(GOR->Link[i]->tail,GOR->Link[i]->head)<<"   "<<G->dijkstra(GOR->Link[i]->tail,GOR->Link[i]->head,0)<<endl;
				cout << del1 << " delay "<<del2<<endl;
				cout <<"util "<< util1 <<endl;
				exit(1);
			}*/

			if(del1 >=INF && del2>=INF){
				ability = INF;
				return;
			}

			if( del1 < del2){
				ability = util1;
				delay = del1;
				if(can + SMALL >= INF)
					ability = INF;
			}
			else{
				//printf("overlay  adjust\n");
				delay = del2;
				vector<demand> req;
				for(int i = 0;i <GOR->m;i++){
					if(GOR->Link[i]->use > 0)
						req.push_back(demand(GOR->Link[i]->tail,GOR->Link[i]->head,GOR->Link[i]->use));
				}
				for( int i = ORNUM; i<(*dem).size(); i++)
					req.push_back((*dem)[i]);

				G->clearOcc();
				for(int i = 0; i < req.size(); i++){
					double dis = G->dijkstraLB(i,req[i].org, req[i].des, req[i].flow);
					can = max(can, dis);
				}

				double util2 = 0;
				for(int i = 0; i < G->m; i++)
					util2 = max(util2, (double)G->Link[i]->use);
				
				ability = util2;
				if(can + SMALL >= INF)
					ability = INF;
			}		
		}

		void randNature(){
			for(int i = 0; i < x.size(); i++){
				x[i] = 1.0 * HOR * rand() / RAND_MAX;
				if(rand()%2)
					x[i] = -x[i];
				x[i] += G->Link[i]->dist;
				x[i] = max(0.0, x[i]);
				x[i] = min(1.0*MAXWEIGHT, x[i]);
			}
		}

		void mutation(){  //变异
			for(int i = 0; i < x.size(); i++){
				x[i] += -(MAXWEIGHT * MUT/100.0) + 2 * (MAXWEIGHT * MUT/100.0) * rand() / RAND_MAX;
				x[i] = max(0.0, x[i]);
				x[i] = min(1.0 * MAXWEIGHT, x[i]);
			}
		}
		void culture(evoluDiv hero){
			for(int i = 0; i < x.size(); i++){
				double r = (CUL/100.0) * rand() / RAND_MAX;
				x[i] = (1 - r) * x[i] + r * hero.x[i];
			}
		}
};

bool evoluCmp(evoluDiv a, evoluDiv b){
	return a.ability < b.ability;
}

class evoluPopu{
	private:
		static const int YEAR = 100;
		static const int YEARCUL = 50;
		vector<evoluDiv> popu;
		CGraph *G;
		CGraph *GOR;
		vector<demand> *dem;
		double pm;
	public:
		evoluDiv hero;
		evoluPopu(int n, int m, CGraph *g, CGraph *gor,vector<demand> *d ){
			pm = 0.25;
			popu.clear();
			G = g;
			GOR = gor;
			dem = d;
			for(int i = 0; i < n; i++){
				evoluDiv divi(m, G, GOR,dem);
				popu.push_back(divi);
			}
			hero = evoluDiv(m, G, GOR,dem);
		}
		evoluDiv evolution(){
			for(int i = 0; i < hero.x.size(); i++)
				hero.x[i] = G->Link[i]->dist;
			
			hero.calAbility();
			for(int i = 0; i < popu.size(); i++){
				popu[i].calAbility();
			}
			
			sort(popu.begin(), popu.end(), evoluCmp);
			

			for(int curYear = 1; curYear <= YEAR; curYear++){
				int n = popu.size(), getMore = 0;
				vector<evoluDiv> sons;
				for(int i = 0; i+1 < n; i+=2){
					sons.push_back(popu[i].mate(popu[i+1]));
					sons.push_back(popu[i+1].mate(popu[i]));
					sons.push_back(popu[i].mate(popu[i+1]));
				}
				int m = sons.size();
				for(int i = 0; i < m; i++){
					double p = rand()%100*0.01;
					if( p < pm ) 
						sons[i].mutation();
					if(curYear > YEARCUL)
						sons[i].culture(hero);
					sons[i].calAbility();
				}
				sort(sons.begin(), sons.end(), evoluCmp);
				popu.clear();
				for(int i = 0; i < n; i++){
					popu.push_back(sons[i]);
					if(sons[i].ability < hero.ability){
						hero = sons[i];
						getMore = 1;
					}
				}
				if(getMore){
					;
					//printf("Year %d: find hero \n", curYear);
					//printf("%f\n", hero.ability);
				}
			}
			return hero;
		}
};

#endif
