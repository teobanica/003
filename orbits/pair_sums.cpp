// Independent enumeration using all unordered pair sums, including doubles.
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
int main(int argc, char** argv) {
    const std::vector<int> supported{5,7,11,13,17,19,23,29,31};
    const int p = argc == 2 ? std::atoi(argv[1]) : 0;
    if (std::find(supported.begin(), supported.end(), p) == supported.end()) {
        std::fprintf(stderr, "Choose an odd prime between 5 and 31.\n");
        return 2;
    }
    const int m=(p-1)/2;
    std::vector<int> squares, word{1};
    for(int i=1;i<p;++i)squares.push_back(i*i%p);
    std::sort(squares.begin(),squares.end());
    squares.erase(std::unique(squares.begin(),squares.end()),squares.end());
    std::array<std::uint32_t,16> sums{};
    sums[0]=1u<<(2%p);
    auto extend=[&](auto&& extend, std::uint32_t used){
        int j=word.size();
        if(j==m){
            for(int a:word)std::printf("%d ",a);
            std::puts("");
            return;
        }
        for(int x:squares){
            if(used&(1u<<x))continue;
            auto fresh=sums;
            bool ok=true;
            for(int i=0;i<=j && ok;++i){
                auto bit=1u<<((x+(i==j?x:word[i]))%p);
                auto& target=fresh[(i+j)%m];
                ok=!(target&bit);
                target|=bit;
            }
            if(!ok)continue;
            std::swap(fresh,sums);
            word.push_back(x); extend(extend,used|(1u<<x)); word.pop_back();
            std::swap(fresh,sums);
        }
    };
    extend(extend,2);
}
