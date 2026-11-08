#include <cstdio>
#include <vector>
#include <algorithm>
using namespace std;
int main() {
    const int orders[]={19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103};
    for(int p:orders){
        vector<int> Q(p),inv(p);
        for(int i=1;i<p;i++){Q[(i*i)%p]=1;for(int j=1;j<p;j++)if(i*j%p==1)inv[i]=j;}
        vector<int> a{0};
        for(int i=2;i<p;i++)if(Q[i])a.push_back(i);
        vector<int> f{0,a[1]};long nodes=0,solutions=0;
        auto visit=[&](auto&& visit)->void{
            nodes++;
            int k=int(f.size());
            if(k==int(a.size())){if(f!=a){puts("nonlinear");}solutions++;return;}
            for(int t=0;t<p;t++){
                bool good=true;
                for(int j=0;j<k;j++){
                    int u=(p+t-f[j])*inv[a[k]-a[j]]%p;
                    if(u==0||(Q[u]&&u!=1)){good=false;break;}
                }
                if(good){f.push_back(t);visit(visit);f.pop_back();}
            }
        };
        visit(visit);printf("%d %ld %ld\n",p,nodes,solutions);fflush(stdout);
    }
}
