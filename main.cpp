#include <vector>
#include <cstdint>
#include <complex>
#include <bitset>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>

constexpr std::uint64_t MAX_OPERANDS = 3;
constexpr std::uint8_t SENTINEL = -1;

struct Entry {
    std::uint32_t gateIndex = 0;
    std::array<std::uint8_t, MAX_OPERANDS> operands = {SENTINEL, SENTINEL, SENTINEL};
    std::uint8_t i = 0;
    std::uint8_t j = 0;
    std::complex<double> c;
};

using Circuit = std::vector<Entry>;
using Ket = std::bitset<256>;
using State = std::unordered_map<Ket, std::complex<double>>;
using BiggestAmplitudes = std::vector<std::pair<Ket, std::complex<double>>>;

struct StackEntry {
    Circuit::const_iterator iterator;
    std::complex<double> c;
};

using Stack = std::vector<StackEntry>;

inline bool matchesEntry(Entry entry, Ket ket) {
    assert(!ket[SENTINEL]);

    auto ketJ = 0;
    assert(!ket.test(entry.operands[SENTINEL]));
    for (auto opIndex = 0; opIndex < MAX_OPERANDS; ++opIndex) {
        ketJ |= (ket.test(entry.operands[opIndex]) << opIndex);
    }

    return entry.j == ketJ;
}

void down(Stack& stack, Ket& ket, Circuit const& circuit) {
    while (true) {
        auto it = stack.empty() ? circuit.cbegin() :
            std::find_if(stack.back().iterator, circuit.cend(), [&](auto const& x) {
                assert(stack.back().iterator != circuit.cend());
                return x.gateIndex != stack.back().iterator->gateIndex && matchesEntry(x, ket);
            });
        
        if (it == circuit.cend()) {
            return;
        }

        stack.push_back({ .iterator = it, .c = stack.empty() ? it->c : (stack.back().c * it->c) });

        auto i = it->i;
        for (auto op: it->operands) {
            assert(op != SENTINEL || i == 0);
            
            ket[op] = (i & 1U);
            i >>= 1U;
        }
    }
}

void upAndAdvance(Stack& stack, Ket& ket, Circuit const& circuit) {
    if (stack.empty()) {
        return;
    }

    assert(stack.back().iterator != circuit.end());

    auto gateIndex = stack.back().iterator->gateIndex;
    auto currentJ = stack.back().iterator->j;
    auto newIterator = std::next(stack.back().iterator);
    while (newIterator != circuit.end() && newIterator->gateIndex == gateIndex && newIterator->j != currentJ) {
        ++newIterator;
    }

    if (newIterator != circuit.end() && newIterator->gateIndex == gateIndex) {
        stack.back().iterator = newIterator;
        if (stack.size() >= 2) {
            stack.back().c = std::next(stack.rbegin())->c * newIterator->c;
        } else {
            stack.back().c = newIterator->c;
        }

        auto i = stack.back().iterator->i;
        for (auto op: stack.back().iterator->operands) {
            assert(op != SENTINEL || i == 0);

            ket[op] = (i & 1U);
            i >>= 1;
        }
    } else {
        for (auto op: stack.back().iterator->operands) {
            assert(op != SENTINEL || currentJ == 0);

            ket[op] = (currentJ & 1U);
            currentJ >>= 1;
        }
        stack.pop_back();
        upAndAdvance(stack, ket, circuit);
    }
}

BiggestAmplitudes getMostLikelyAmplitudes(Circuit circuit) {
    auto numberOfGates = circuit.back().gateIndex + 1;
    std::uint8_t qubitCount = 0;
    for (auto const& entry: circuit) {
        for (auto op: entry.operands) {
            if (op != SENTINEL) {
                qubitCount = std::max(qubitCount, static_cast<std::uint8_t>(op + 1));
            }
        }
    }

    State currentPartialState;

    std::uint64_t nBiggestAmplitudes = 10;
    BiggestAmplitudes biggestAmplitudes;
    biggestAmplitudes.reserve(nBiggestAmplitudes);

    Stack stack;
    stack.reserve(numberOfGates);

    Ket ket;

    std::uint64_t rightBits = 10;
    std::uint64_t maxBitsInBucket = (qubitCount <= rightBits) ? 0 : (qubitCount - rightBits);
    std::uint64_t maxBucket = (maxBitsInBucket < 64) ? ((1 << maxBitsInBucket) - 1) : ~0ULL;

    auto getBucket = [rightBits](Ket ket) -> std::uint64_t {
        Ket mask;
        mask.set();
        mask >>= 256 - 64;
        return ((ket >> rightBits) & mask).to_ullong();
    };

    for (std::uint64_t bucket = 0; bucket <= maxBucket; ++bucket) {
        assert(currentPartialState.empty());

        do {
            down(stack, ket, circuit);
            
            assert(stack.size() <= numberOfGates);

            if (stack.size() == numberOfGates && getBucket(ket) == bucket) {
                auto [it, inserted] = currentPartialState.insert({ket, stack.back().c});
                if (!inserted) {
                    it->second += stack.back().c;
                    if (std::norm(it->second) < 0.000001) {
                        currentPartialState.erase(it);
                    }
                }
            }

            upAndAdvance(stack, ket, circuit);
        } while (!stack.empty());

        BiggestAmplitudes partialBiggestAmplitudes(nBiggestAmplitudes, std::make_pair(Ket{}, 0.));
        auto top = std::partial_sort_copy(currentPartialState.begin(), currentPartialState.end(), partialBiggestAmplitudes.begin(), partialBiggestAmplitudes.end(),
            [](auto left, auto right) {
                return std::norm(left.second) > std::norm(right.second);
            });
        
        BiggestAmplitudes outputBiggestAmplitudes;
        std::merge(biggestAmplitudes.begin(), biggestAmplitudes.end(), partialBiggestAmplitudes.begin(), top, std::back_inserter(outputBiggestAmplitudes),
            [](auto left, auto right) {
                return std::norm(left.second) > std::norm(right.second);
            });

        biggestAmplitudes.swap(outputBiggestAmplitudes);

        currentPartialState.clear();
    }

    return biggestAmplitudes;
}

class CircuitBuilder {
public:
    CircuitBuilder& x(std::uint8_t operand) {
        push({operand}, {{0, 1, 1.}, {1, 0, 1.}});
        return *this;
    }

    CircuitBuilder& h(std::uint8_t operand) {
        auto q = sqrt(0.5);
        push({operand}, {{0, 0, q}, {0, 1, q}, {1, 0, q}, {1, 1, -q}});
        return *this;
    }

    CircuitBuilder& cnot(std::uint8_t operand1, std::uint8_t operand2) {
        push({operand1, operand2}, {{0, 0, 1.}, {1, 1, 1.}, {2, 3, 1.}, {3, 2, 1.}});
        return *this;
    }

    Circuit get() {
        return circuit;
    }

private:
    void push(std::initializer_list<std::uint8_t> operands, std::initializer_list<std::tuple<std::uint8_t, std::uint8_t, std::complex<double>>> m) {
        assert(operands.size() <= MAX_OPERANDS);
        std::array<std::uint8_t, MAX_OPERANDS> operandsArray = {SENTINEL, SENTINEL, SENTINEL};
        std::uint64_t i = operands.size() - 1;
        for (auto op: operands) {
            operandsArray[i] = op;
            --i;
        }

        for (auto [i, j, c]: m) {
            circuit.push_back({
                .gateIndex = gateIndex,
                .operands = operandsArray,
                .i = i,
                .j = j,
                .c = c
            });
        }

        ++gateIndex;
    }

    std::uint32_t gateIndex = 0;
    Circuit circuit;
};

int main() {
    auto circuit = CircuitBuilder().h(0).h(1).h(2).h(3).h(4).h(5).get();

    auto result = getMostLikelyAmplitudes(circuit);

    for (auto [k, v]: result) {
        std::cout << k << "   ->    " << v << std::endl;
    }

    return 0;
}