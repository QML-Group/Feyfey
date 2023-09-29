#include <vector>
#include <cstdint>
#include <complex>
#include <bitset>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>

constexpr std::uint64_t MAX_OPERANDS = 3;
constexpr std::uint8_t SENTINEL = -1;
constexpr std::uint64_t N_BIGGEST_AMPLITUDES = 10;

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
using BiggestAmplitudes = std::array<std::pair<Ket, std::complex<double>>, N_BIGGEST_AMPLITUDES>;

struct StackEntry {
    Circuit::const_iterator iterator;
    std::complex<double> c;
};

using Stack = std::vector<StackEntry>;

inline bool matchesEntry(Entry entry, Ket ket) {
    assert(!ket[SENTINEL]);

    auto ketJ = 0;
    assert(!ket.test(SENTINEL));
    for (auto opIndex = 0; opIndex < MAX_OPERANDS; ++opIndex) {
        ketJ |= (ket.test(entry.operands[opIndex]) << opIndex);
    }

    return entry.j == ketJ;
}

void down(Stack& stack, Ket& ket, Circuit const& circuit) {
    while (true) {
        auto it = stack.empty() ? std::find_if(circuit.begin(), circuit.cend(), [ket](auto const& x) {
                return matchesEntry(x, ket);
                }) :
            std::find_if(stack.back().iterator, circuit.cend(), [&](auto const& x) {
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

std::pair<BiggestAmplitudes, std::size_t> getMostLikelyAmplitudes(Circuit circuit) {
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

    BiggestAmplitudes biggestAmplitudes;

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

    std::size_t numberBiggestAmplitudes = 0;

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

        BiggestAmplitudes partialBiggestAmplitudes;
        auto top = std::partial_sort_copy(currentPartialState.begin(), currentPartialState.end(), partialBiggestAmplitudes.begin(), partialBiggestAmplitudes.end(),
            [](auto left, auto right) {
                return std::norm(left.second) > std::norm(right.second);
            });
        
        numberBiggestAmplitudes = std::min(N_BIGGEST_AMPLITUDES, std::distance(partialBiggestAmplitudes.begin(), top) + numberBiggestAmplitudes);
        
        std::array<std::pair<Ket, std::complex<double>>, 2 * N_BIGGEST_AMPLITUDES> outputBiggestAmplitudes;
        std::merge(biggestAmplitudes.begin(), biggestAmplitudes.end(), partialBiggestAmplitudes.begin(), top, outputBiggestAmplitudes.begin(),
            [](auto left, auto right) {
                return std::norm(left.second) > std::norm(right.second);
            });

        std::copy(outputBiggestAmplitudes.begin(), std::next(outputBiggestAmplitudes.begin(), N_BIGGEST_AMPLITUDES), biggestAmplitudes.begin());

        currentPartialState.clear();
    }

    return { biggestAmplitudes, numberBiggestAmplitudes };
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

    CircuitBuilder& random(std::uint8_t operand) {
        // This ignores Haar measure - just for benchmarking.

        double alpha = dis(re);
        double beta = dis(re);
        double delta = dis(re);
        double gamma = dis(re);

        auto s = std::sin(gamma / 2.);
        auto c = std::cos(gamma / 2.);

        push({operand}, {
            {0, 0,   std::polar(1., alpha - beta / 2 - delta / 2) * c},
            {0, 1, - std::polar(1., alpha - beta / 2 + delta / 2) * s},
            {1, 0,   std::polar(1., alpha + beta / 2 - delta / 2) * s},
            {1, 1,   std::polar(1., alpha + beta / 2 + delta / 2) * c}
        });
        return *this;
    }

    static Circuit createRandom(std::uint8_t num_qubits, std::uint64_t num_gates, double probaOfCnot) {
        assert(0 <= probaOfCnot && probaOfCnot <= 1);

        CircuitBuilder builder;
        std::default_random_engine re;
        std::bernoulli_distribution bernoulli(probaOfCnot);
        std::uniform_int_distribution<std::uint8_t> operands(0, num_qubits);

        while (num_gates > 0) {
            auto operand1 = operands(re);
            auto operand2 = operands(re);

            if (bernoulli(re)) {
                builder.cnot(operand1, operand2);
            } else {
                builder.random(operand1);
            }

            --num_gates;
        }

        return builder.get();
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

    std::uniform_real_distribution<double> dis{0, 1000};
    std::default_random_engine re; // Always the same seed here.

    std::uint32_t gateIndex = 0;
    Circuit circuit;
};

void check(Circuit circuit, Ket ket, std::complex<double> expectedAmplitude) {
    auto [result, nAmplitudes] = getMostLikelyAmplitudes(circuit);

    auto max = std::next(result.begin(), nAmplitudes);

    auto it = std::find_if(result.begin(), max, [ket](auto x) {
        return x.first == ket;
    });

    if (it == max) {
        throw std::runtime_error("[TEST FAILURE] Did not found amplitude for " + ket.to_string());
    }

    if (std::norm(it->second - expectedAmplitude) > 0.0000001) {
        throw std::runtime_error("[TEST FAILURE] Wrong amplitude for " + ket.to_string() + ": got " + std::to_string(it->second.real()) + " + i" + std::to_string(it->second.imag()) + " but expected "
            + std::to_string(expectedAmplitude.real()) + " + i" + std::to_string(expectedAmplitude.imag()) );
    }
}

void checkNot(Circuit circuit, Ket ket) {
    auto [result, nAmplitudes] = getMostLikelyAmplitudes(circuit);

    auto max = std::next(result.begin(), nAmplitudes);

    auto it = std::find_if(result.begin(), max, [ket](auto x) {
        return x.first == ket;
    });

    if (it != max) {
        throw std::runtime_error("[TEST FAILURE] Expected no amplitude but got one for " + ket.to_string());
    }
}

void test() {
    using C = CircuitBuilder;

    check(C().h(0).get(), Ket("01"), sqrt(0.5));
    checkNot(C().h(0).get(), Ket("10"));
    check(C().h(0).h(1).get(), Ket("11"), 0.5);
    check(C().h(0).h(1).get(), Ket("01"), 0.5);
    check(C().h(0).h(0).get(), Ket("0"), 1.);
    checkNot(C().h(0).h(0).get(), Ket("1"));
    check(C().x(4).h(4).get(), Ket("00000"), sqrt(0.5));
    check(C().x(4).h(4).get(), Ket("10000"), -sqrt(0.5));
    check(C().h(0).cnot(0, 1).get(), Ket("11"), sqrt(0.5));
    check(C().x(0).get(), Ket("1"), 1.);
    checkNot(C().x(0).get(), Ket("0"));
    check(C().x(2).get(), Ket("100"), 1.);
    checkNot(C().x(2).get(), Ket("101"));
    check(C().x(3).cnot(3, 5).get(), Ket("101000"), 1.);
    check(C().x(2).cnot(3, 5).get(), Ket("000100"), 1.);
}

int main() {
    test();

    std::uint8_t num_qubits = 20;
    std::uint64_t num_gates = 40;
    double probaOfCnot = 0.5;

    auto circuit = CircuitBuilder::createRandom(num_qubits, num_gates, probaOfCnot);

    auto [result, count] = getMostLikelyAmplitudes(circuit);

    for (auto it = result.begin(); it != std::next(result.begin(), count); ++it) {
        std::cout << it->first << "   ->    " << it->second << std::endl;
    }

    return 0;
}