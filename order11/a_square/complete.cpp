#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
int n,words;
vector<array<int,11>> rows;
vector<vector<uint64_t>> adj;
long long nodes=0,found=0,edges=0;
vector<int> clique;
vector<vector<int>> solutions;
inline bool connected(int a,int b){return (adj[a][b/64]>>(b%64))&1;}
// Greedy independent-set coloring gives a rigorous upper bound on clique
// size. Each color class is explicitly checked to be independent by the
// adjacency operations in color_order.
void color_order(const vector<int>& candidates,vector<int>& order,vector<int>& bounds){
 vector<int> remaining=candidates;
 int color=0;
 while(!remaining.empty()){
  ++color;vector<int> colorclass,next;
  for(int v:remaining){
   bool allowed=true;
   for(int u:colorclass)if(connected(v,u)){allowed=false;break;}
   if(allowed){colorclass.push_back(v);order.push_back(v);bounds.push_back(color);}
   else next.push_back(v);
  }
  remaining.swap(next);
 }
}
void search(const vector<int>& candidates){
 ++nodes;
 if(clique.size()==9){++found;solutions.push_back(clique);return;}
 if(clique.size()+candidates.size()<9)return;
 vector<int> order,bounds;color_order(candidates,order,bounds);
 for(int pos=(int)order.size()-1;pos>=0;--pos){
  if(clique.size()+bounds[pos]<9)return;
  int v=order[pos];vector<int> next;
  for(int k=0;k<pos;++k)if(connected(v,order[k]))next.push_back(order[k]);
  clique.push_back(v);search(next);clique.pop_back();
 }
}
int main(int argc,char**argv){
 if(argc!=3){cerr<<"Usage: complete rows.txt output.json\n";return 2;}
 ifstream input(argv[1]);input>>n;rows.resize(n);
 for(auto &r:rows)for(int &a:r)input>>a;
 words=(n+63)/64;adj.assign(n,vector<uint64_t>(words));
 for(int a=0;a<n;++a)for(int b=a+1;b<n;++b){
  unsigned used=0;bool ok=true;
  for(int j=0;j<11;++j){int d=(rows[a][j]-rows[b][j]+11)%11;if(used>>d&1){ok=false;break;}used|=1u<<d;}
  if(ok){adj[a][b/64]|=1ull<<(b%64);adj[b][a/64]|=1ull<<(a%64);++edges;}
 }
 vector<int> all(n);iota(all.begin(),all.end(),0);
 cout<<"vertices "<<n<<", edges "<<edges<<endl;
 search(all);
 cout<<"search nodes "<<nodes<<", cliques of size 9 "<<found<<endl;
 ofstream out(argv[2]);out<<"{\n  \"vertices\": "<<n<<",\n  \"edges\": "<<edges<<",\n  \"search_nodes\": "<<nodes<<",\n  \"completion_cliques\": [";
 for(size_t i=0;i<solutions.size();++i){if(i)out<<",";out<<"\n    [";for(size_t j=0;j<solutions[i].size();++j){if(j)out<<",";out<<solutions[i][j];}out<<"]";}
 out<<"\n  ]\n}\n";
 for(const auto& sol:solutions){cout<<"clique:";for(int v:sol)cout<<' '<<v;cout<<endl;}
}
