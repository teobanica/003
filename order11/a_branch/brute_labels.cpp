// Independent Cartesian-product count. No phase solver or incremental
// all-different pruning: check all 11 differences at each full assignment.
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <utility>
int main(){
 int cases;std::cin>>cases;
 for(int cas=0;cas<cases;++cas){
  std::string name;int vertices,left,fixed;long long expected;
  std::cin>>name>>vertices>>left>>fixed>>expected;
  std::vector<int> value(vertices,-1),free;
  for(int i=0;i<fixed;++i){int v,a;std::cin>>v>>a;value[v]=a;}
  std::vector<std::pair<int,int>> edges;
  for(int e=0;e<11;++e){int a,b;std::cin>>a>>b;edges.push_back({a,b});}
  for(int v=0;v<vertices;++v)if(value[v]<0)free.push_back(v);
  long long count=0,assignments=0;
  std::function<void(int)> run=[&](int depth){
   if(depth==(int)free.size()){
    ++assignments;unsigned used=0;
    for(auto [i,j]:edges){int d=(value[j]-value[i]+11)%11;if(used&(1u<<d))return;used|=1u<<d;}
    if(used==2047)++count;return;
   }
   for(int a=0;a<11;++a){value[free[depth]]=a;run(depth+1);}
  };
  run(0);
  std::cout<<name<<" graph "<<cas+1<<": "<<count<<" normalized labels, "<<assignments<<" assignments; expected "<<expected<<std::endl;
  if(count!=expected)return 1;
 }
}
