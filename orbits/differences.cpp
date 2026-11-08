// Enumerate square orderings with a_0=1, using directed differences.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char** argv) {
    const std::vector<int> supported{5, 7, 11, 13, 17, 19, 23, 29, 31};
    const int p = argc == 2 ? std::atoi(argv[1]) : 0;
    if (std::find(supported.begin(), supported.end(), p) == supported.end()) {
        std::fprintf(stderr, "Choose an odd prime between 5 and 31.\n");
        return 2;
    }
    const int m = (p - 1) / 2;
    std::vector<int> squares;
    for (int j = 1; j <= m; ++j) squares.push_back(j * j % p);
    std::sort(squares.begin(), squares.end());
    std::vector<int> word{1};
    using Masks = std::array<std::uint32_t, 16>;
    auto visit = [&](auto&& visit, std::uint32_t used, const Masks& gaps) -> void {
        int j = static_cast<int>(word.size());
        if (j == m) {
            for (int value : word) std::printf("%d ", value);
            std::puts("");
            return;
        }
        for (int x : squares) {
            if ((used >> x) & 1u) continue;
            Masks next = gaps;
            bool unique = true;
            for (int i = 0; i < j && unique; ++i) {
                const int displacement = j - i;
                const int difference = (p + x - word[i]) % p;
                const auto forward = 1u << difference;
                const auto backward = 1u << (p - difference);
                if ((next[displacement] & forward) || (next[m-displacement] & backward)) {
                    unique = false;
                    break;
                }
                next[displacement] |= forward;
                next[m-displacement] |= backward;
            }
            if (unique) {
                word.push_back(x);
                visit(visit, used | (1u << x), next);
                word.pop_back();
            }
        }
    };
    visit(visit, 1u << 1, Masks{});
}
