#include <complex>
#include <array>
#include <vector>
#include <bitset>
#include <climits>
#include <cassert>
#include <optional>
#include <unordered_map>
#include <algorithm>
#include <cmath>

template <std::uint64_t N>
using Matrix = std::array<std::array<std::complex<double>, N>, N>;

using Operand = std::uint8_t;
using MatrixIndex = std::uint8_t;
using GateIndex = std::uint32_t;

struct GateWithOperands {
    Matrix<4> matrix;
    Operand operand1 = 0;
    Operand operand2 = 0;
};

using Circuit = std::vector<GateWithOperands>;

struct Incompatibility {
    GateIndex gateIndex = 0;
    MatrixIndex i = 0;
    MatrixIndex j = 0;
};

static constexpr std::uint64_t MAX_NUMBER_OF_QUBITS = std::numeric_limits<Operand>::max();
static constexpr std::uint64_t MAX_NUMBER_OF_GATES = 1000;

using Ket = std::bitset<MAX_NUMBER_OF_QUBITS>;

using SetOfGates = std::bitset<MAX_NUMBER_OF_GATES>;
using SetOfMatrixElements = std::bitset<4 * 4 * MAX_NUMBER_OF_GATES>;

static constexpr double ATOL = 1e-7;

template <typename T>
bool isNotNull(T t) {
    return std::norm(t) > ATOL;
}

template <typename T1, typename T2>
bool fequal(T1 t1, T2 t2) {
    return std::norm(t1 - t2) < ATOL;
}

SetOfMatrixElements getElements(Circuit const& circuit, Ket ket) {
    std::bitset<MAX_NUMBER_OF_QUBITS> allowedZero;
    allowedZero.set();

    std::bitset<MAX_NUMBER_OF_QUBITS> allowedOne;

    std::bitset<MAX_NUMBER_OF_QUBITS> allOperands;

    SetOfMatrixElements result;

    auto set = [&circuit, &result](GateIndex gateIndex, MatrixIndex i, MatrixIndex j) {
        bool shouldSet = std::norm(circuit[gateIndex].matrix[i][j]) > ATOL;
        result.set(16 * gateIndex + i * 4 + j, shouldSet);
        return shouldSet;
    };

    for (std::uint64_t gateIndex = 0; gateIndex < circuit.size(); ++gateIndex) {
        auto const& g = circuit[gateIndex];
        std::bitset<4> nonNullRows;

        allOperands.set(g.operand1);
        allOperands.set(g.operand2);

        if (allowedZero.test(g.operand1) && allowedZero.test(g.operand2)) {
            if (set(gateIndex, 0b00, 0b00)) {
                nonNullRows.set(0b00);
            }

            if (set(gateIndex, 0b01, 0b00)) {
                nonNullRows.set(0b01);
            }
            
            if (set(gateIndex, 0b10, 0b00)) {
                nonNullRows.set(0b10);
            }
            
            if (set(gateIndex, 0b11, 0b00)) {
                nonNullRows.set(0b11);
            }
        }

        if (allowedZero.test(g.operand1) && allowedOne.test(g.operand2)) {
            if (set(gateIndex, 0b00, 0b01)) {
                nonNullRows.set(0b00);
            }

            if (set(gateIndex, 0b01, 0b01)) {
                nonNullRows.set(0b01);
            }
            
            if (set(gateIndex, 0b10, 0b01)) {
                nonNullRows.set(0b10);
            }
            
            if (set(gateIndex, 0b11, 0b01)) {
                nonNullRows.set(0b11);
            }
        }
        
        if (allowedOne.test(g.operand1) && allowedZero.test(g.operand2)) {
            if (set(gateIndex, 0b00, 0b10)) {
                nonNullRows.set(0b00);
            }

            if (set(gateIndex, 0b01, 0b10)) {
                nonNullRows.set(0b01);
            }
            
            if (set(gateIndex, 0b10, 0b10)) {
                nonNullRows.set(0b10);
            }
            
            if (set(gateIndex, 0b11, 0b10)) {
                nonNullRows.set(0b11);
            }
        }
        
        if (allowedOne.test(g.operand1) && allowedOne.test(g.operand2)) {
            if (set(gateIndex, 0b00, 0b11)) {
                nonNullRows.set(0b00);
            }

            if (set(gateIndex, 0b01, 0b11)) {
                nonNullRows.set(0b01);
            }
            
            if (set(gateIndex, 0b10, 0b11)) {
                nonNullRows.set(0b10);
            }
            
            if (set(gateIndex, 0b11, 0b11)) {
                nonNullRows.set(0b11);
            }
        }

        allowedZero.set(g.operand1, nonNullRows.test(0b00) || nonNullRows.test(0b01));
        allowedZero.set(g.operand2, nonNullRows.test(0b00) || nonNullRows.test(0b10));
        allowedOne.set(g.operand1, nonNullRows.test(0b10) || nonNullRows.test(0b11));
        allowedOne.set(g.operand2, nonNullRows.test(0b01) || nonNullRows.test(0b11));
    }

    if ((ket & (~allOperands)).any()) {
        return SetOfMatrixElements();
    }

    allowedZero = ~ket;
    allowedOne = ket;

    auto reset = [&result](GateIndex gateIndex, MatrixIndex i, MatrixIndex j) {
        result.reset(16 * gateIndex + 4 * i + j);
    };

    auto resetRow = [&reset, &result](GateIndex gateIndex, MatrixIndex rowIndex) {
        reset(gateIndex, rowIndex, 0b00);
        reset(gateIndex, rowIndex, 0b01);
        reset(gateIndex, rowIndex, 0b10);
        reset(gateIndex, rowIndex, 0b11);
    };
    
    auto colIsNonNull = [&circuit, &result](GateIndex gateIndex, MatrixIndex colIndex) {
        return isNotNull(circuit[gateIndex].matrix[0b00][colIndex] * static_cast<double>(result.test(16 * gateIndex + 4 * 0b00 + colIndex))) ||
            isNotNull(circuit[gateIndex].matrix[0b01][colIndex] * static_cast<double>(result.test(16 * gateIndex + 4 * 0b01 + colIndex))) ||
            isNotNull(circuit[gateIndex].matrix[0b10][colIndex] * static_cast<double>(result.test(16 * gateIndex + 4 * 0b10 + colIndex))) ||
            isNotNull(circuit[gateIndex].matrix[0b11][colIndex] * static_cast<double>(result.test(16 * gateIndex + 4 * 0b11 + colIndex)));
    };

    for (std::uint64_t k = 0; k < circuit.size(); ++k) {
        std::uint64_t gateIndex = circuit.size() - k - 1;

        auto const& g = circuit[gateIndex];
        std::bitset<4> nonNullCols;

        if (!allowedZero.test(g.operand1)) {
            resetRow(gateIndex, 0b00);
            resetRow(gateIndex, 0b01);
        }
        
        if (!allowedZero.test(g.operand2)) {
            resetRow(gateIndex, 0b00);
            resetRow(gateIndex, 0b10);
        }

        if (!allowedOne.test(g.operand1)) {
            resetRow(gateIndex, 0b10);
            resetRow(gateIndex, 0b11);
        }

        if (!allowedOne.test(g.operand2)) {
            resetRow(gateIndex, 0b01);
            resetRow(gateIndex, 0b11);
        }

        allowedZero.reset(g.operand1);
        allowedZero.reset(g.operand2);
        allowedOne.reset(g.operand1);
        allowedOne.reset(g.operand2);

        if (colIsNonNull(gateIndex, 0b00)) {
            allowedZero.set(g.operand1, true);
            allowedZero.set(g.operand2, true);
        }

        if (colIsNonNull(gateIndex, 0b01)) {
            allowedZero.set(g.operand1, true);
            allowedOne.set(g.operand2, true);
        }

        if (colIsNonNull(gateIndex, 0b10)) {
            allowedOne.set(g.operand1, true);
            allowedZero.set(g.operand2, true);
        }

        if (colIsNonNull(gateIndex, 0b11)) {
            allowedOne.set(g.operand1, true);
            allowedOne.set(g.operand2, true);
        }

        if (!allowedZero.test(g.operand1) && !allowedOne.test(g.operand1)) {
            return SetOfMatrixElements();
        }

        if (!allowedZero.test(g.operand2) && !allowedOne.test(g.operand2)) {
            return SetOfMatrixElements();
        }
    }

    return result;
}

using Incompatibilities = std::vector<SetOfMatrixElements>;

struct Partition {
    SetOfGates left;
    SetOfMatrixElements leftKeys;
    SetOfGates right;
    SetOfMatrixElements rightKeys;
    Incompatibilities incompatibilities;
};

Partition split(Circuit const& circuit, SetOfMatrixElements elements, SetOfGates left) {
    assert((left >> circuit.size()).none());

    SetOfGates right;
    auto half = left.count() / 2;
    std::size_t pos = 0;
    while (left.count() > half) {
        if (left.test(pos)) {
            left.reset(pos);
            right.set(pos);
        }
        ++pos;
    }

    std::vector<SetOfMatrixElements> incompatibilities;

    SetOfMatrixElements leftKeys;
    SetOfMatrixElements rightKeys;

    std::array<std::optional<GateIndex>, MAX_NUMBER_OF_QUBITS> qubitToLastGate;
    for (std::uint64_t gateIndex = 0; gateIndex < circuit.size(); ++gateIndex) {
        bool isLeft = left.test(gateIndex);
        bool isRight = right.test(gateIndex);

        if (!isLeft && !isRight) {
            continue;
        }

        auto op1 = circuit[gateIndex].operand1;
        auto op2 = circuit[gateIndex].operand2;

        auto lastGate1 = qubitToLastGate[op1];
        assert(!lastGate1 || left.test(*lastGate1) || right.test(*lastGate1));

        auto lastGate1IsRight = lastGate1 && right.test(*lastGate1);
        auto lastGate1IsLeft = lastGate1 && left.test(*lastGate1);

        auto lastGate2 = qubitToLastGate[op2];
        assert(!lastGate2 || left.test(*lastGate2) || right.test(*lastGate2));

        auto lastGate2IsRight = lastGate2 && right.test(*lastGate2);
        auto lastGate2IsLeft = lastGate2 && left.test(*lastGate2);

        SetOfMatrixElements maskG = SetOfMatrixElements(0b1111111111111111) << (16 * gateIndex);

        if ((isLeft && lastGate1IsRight) || (isRight && lastGate1IsLeft)) {
            auto lastG1index = *lastGate1;
            auto lastG = circuit[lastG1index];

            SetOfMatrixElements maskLast1 = SetOfMatrixElements(0b1111111111111111) << (16 * lastG1index);

            if (op1 == lastG.operand1) {
                SetOfMatrixElements inc;
                inc |= SetOfMatrixElements(0b1111111100000000) << (16 * lastG1index);
                inc |= SetOfMatrixElements(0b0011001100110011) << (16 * gateIndex);
                inc &= elements;

                if ((inc & maskG).any() && (inc & maskLast1).any()) {
                    leftKeys |= isLeft ? inc & maskG : inc & maskLast1;
                    rightKeys |= isRight ? inc & maskG : inc & maskLast1;
                    incompatibilities.push_back(inc);
                }

                SetOfMatrixElements inc2;
                inc2 |= SetOfMatrixElements(0b0000000011111111) << (16 * lastG1index);
                inc2 |= SetOfMatrixElements(0b1100110011001100) << (16 * gateIndex);
                inc2 &= elements;

                if ((inc2 & maskG).any() && (inc2 & maskLast1).any()) {
                    leftKeys |= isLeft ? inc2 & maskG : inc2 & maskLast1;
                    rightKeys |= isRight ? inc2 & maskG : inc2 & maskLast1;
                    incompatibilities.push_back(inc2);
                }
            } else {
                assert(op1 == lastG.operand2);

                SetOfMatrixElements inc;
                inc |= SetOfMatrixElements(0b1111000011110000) << (16 * lastG1index);
                inc |= SetOfMatrixElements(0b0011001100110011) << (16 * gateIndex);
                inc &= elements;

                if ((inc & maskG).any() && (inc & maskLast1).any()) {
                    leftKeys |= isLeft ? inc & maskG : inc & maskLast1;
                    rightKeys |= isRight ? inc & maskG : inc & maskLast1;
                    incompatibilities.push_back(inc);
                }

                SetOfMatrixElements inc2;
                inc2 |= SetOfMatrixElements(0b0000111100001111) << (16 * lastG1index);
                inc2 |= SetOfMatrixElements(0b1100110011001100) << (16 * gateIndex);
                inc2 &= elements;

                if ((inc2 & maskG).any() && (inc2 & maskLast1).any()) {
                    leftKeys |= isLeft ? inc2 & maskG : inc2 & maskLast1;
                    rightKeys |= isRight ? inc2 & maskG : inc2 & maskLast1;
                    incompatibilities.push_back(inc2);
                }
            }
        }

        if ((isLeft && lastGate2IsRight) || (isRight && lastGate2IsLeft)) {
            auto lastG2index = *lastGate2;
            auto lastG = circuit[lastG2index];

            SetOfMatrixElements maskLast2 = SetOfMatrixElements(0b1111111111111111) << (16 * lastG2index);

            if (op2 == lastG.operand1) {
                SetOfMatrixElements inc;
                inc |= SetOfMatrixElements(0b1111111100000000) << (16 * lastG2index);
                inc |= SetOfMatrixElements(0b0101010101010101) << (16 * gateIndex);
                inc &= elements;

                if ((inc & maskG).any() && (inc & maskLast2).any()) {
                    leftKeys |= isLeft ? inc & maskG : inc & maskLast2;
                    rightKeys |= isRight ? inc & maskG : inc & maskLast2;
                    incompatibilities.push_back(inc);
                }

                SetOfMatrixElements inc2;
                inc2 |= SetOfMatrixElements(0b0000000011111111) << (16 * lastG2index);
                inc2 |= SetOfMatrixElements(0b1010101010101010) << (16 * gateIndex);
                inc2 &= elements;

                if ((inc2 & maskG).any() && (inc2 & maskLast2).any()) {
                    leftKeys |= isLeft ? inc2 & maskG : inc2 & maskLast2;
                    rightKeys |= isRight ? inc2 & maskG : inc2 & maskLast2;
                    incompatibilities.push_back(inc2);
                }
            } else {
                assert(op2 == lastG.operand2);

                SetOfMatrixElements inc;
                inc |= SetOfMatrixElements(0b1111000011110000) << (16 * lastG2index);
                inc |= SetOfMatrixElements(0b0101010101010101) << (16 * gateIndex);
                inc &= elements;

                if ((inc & maskG).any() && (inc & maskLast2).any()) {
                    leftKeys |= isLeft ? inc & maskG : inc & maskLast2;
                    rightKeys |= isRight ? inc & maskG : inc & maskLast2;
                    incompatibilities.push_back(inc);
                }

                SetOfMatrixElements inc2;
                inc2 |= SetOfMatrixElements(0b0000111100001111) << (16 * lastG2index);
                inc2 |= SetOfMatrixElements(0b1010101010101010) << (16 * gateIndex);
                inc2 &= elements;

                if ((inc2 & maskG).any() && (inc2 & maskLast2).any()) {
                    leftKeys |= isLeft ? inc2 & maskG : inc2 & maskLast2;
                    rightKeys |= isRight ? inc2 & maskG : inc2 & maskLast2;
                    incompatibilities.push_back(inc2);
                }
            }
        }


        qubitToLastGate[op1] = gateIndex;
        qubitToLastGate[op2] = gateIndex;
    }

    return Partition{ .left = left, .leftKeys = leftKeys & elements, .right = right, .rightKeys = rightKeys & elements, .incompatibilities = incompatibilities };
}

std::unordered_map<SetOfMatrixElements, std::complex<double>> computeAmplitudeImpl(Circuit const& circuit, SetOfMatrixElements elements, SetOfGates gates, SetOfMatrixElements recursionKeys) {
    std::unordered_map<SetOfMatrixElements, std::complex<double>> result;

    if (gates.count() == 1) {
        std::uint64_t gateIndex = 0;
        while (!gates.test(gateIndex)) {
            ++gateIndex;
        }

        for (std::uint64_t elementIndex = 16 * gateIndex; elementIndex < 16 * (gateIndex + 1); ++elementIndex) {
            if (!elements.test(elementIndex)) {
                continue;
            }

            auto key = SetOfMatrixElements();
            key.set(elementIndex);
            key &= recursionKeys;
            result[key] += circuit[gateIndex].matrix[(elementIndex % 16) / 4][elementIndex % 4];
        }

        return result;
    }

    auto partition = split(circuit, elements, gates);

    auto left = computeAmplitudeImpl(circuit, elements, partition.left, partition.leftKeys | recursionKeys);
    auto right = computeAmplitudeImpl(circuit, elements, partition.right, partition.rightKeys | recursionKeys);

    for (auto [leftKey, leftValue]: left) {
        for (auto [rightKey, rightValue]: right) {
            auto leftAndRight = leftKey | rightKey;

            auto correspondingIncompat = std::find_if(partition.incompatibilities.begin(), partition.incompatibilities.end(),
                [leftAndRight](auto incompat) {
                    return (incompat & leftAndRight).count() == 2;
                });
            
            if (correspondingIncompat != partition.incompatibilities.end()) {
                continue;
            }

            result[leftAndRight & recursionKeys] += leftValue * rightValue;
        }
    }

    return result;
}

std::complex<double> computeAmplitude(Circuit const& circuit, Ket ket) {
    if (circuit.empty()) {
        return ket == Ket() ? 1. : 0.;
    }

    auto elements = getElements(circuit, ket);

    if (elements.none()) {
        return 0.;
    }

    SetOfGates allGates = ~(~SetOfGates() << circuit.size());

    auto result = computeAmplitudeImpl(circuit, elements, allGates, SetOfMatrixElements());

    assert(result.count(SetOfMatrixElements()) == 1);

    return result[SetOfMatrixElements()];
}


Matrix<4> kron(Matrix<2> a, Matrix<2> b) {
    return {{
        {a[0][0] * b[0][0], a[0][0] * b[0][1], a[0][1] * b[0][0], a[0][1] * b[0][1]},
        {a[0][0] * b[1][0], a[0][0] * b[1][1], a[0][1] * b[1][0], a[0][1] * b[1][1]},
        {a[1][0] * b[0][0], a[1][0] * b[0][1], a[1][1] * b[0][0], a[1][1] * b[0][1]},
        {a[1][0] * b[1][0], a[1][0] * b[1][1], a[1][1] * b[1][0], a[1][1] * b[1][1]},
    }};
}

Matrix<2> I{{
    {1, 0},
    {0, 1}
}};

Matrix<2> X{{
    {0, 1},
    {1, 0}
}};

Matrix<2> H{{
    {sqrt(.5), sqrt(.5)},
    {sqrt(.5), -sqrt(.5)}
}};

Matrix<4> CNOT{{
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 0, 1},
    {0, 0, 1, 0}
}};

auto IX = GateWithOperands{ .matrix = kron(I, X), .operand1 = 1, .operand2 = 0 };
auto IH = GateWithOperands{ .matrix = kron(I, H), .operand1 = 1, .operand2 = 0 };
auto CNOT01 = GateWithOperands{ .matrix = {{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}}, .operand1 = 0, .operand2 = 1 };

void test() {
    using C = Circuit;

    auto k00 = Ket("00");
    auto k01 = Ket("01");
    auto k10 = Ket("10");
    auto k11 = Ket("11");

    C c1 = {IX, CNOT01};

    SetOfMatrixElements expectedElements1("0100000000000000" "0000000000010000");
    assert(getElements(c1, k11) == expectedElements1);
    assert(getElements(c1, k00) == SetOfMatrixElements());
    assert(getElements(c1, k01) == SetOfMatrixElements());
    assert(getElements(c1, k10) == SetOfMatrixElements());

    C c2 = {IH, CNOT01};

    SetOfMatrixElements expectedElements2("0000000000000001" "0000000000000001");
    assert(getElements(c2, k00) == expectedElements2);
    SetOfMatrixElements expectedElements3("0100000000000000" "0000000000010000");
    assert(getElements(c2, k11) == expectedElements3);
    assert(getElements(c2, k01) == SetOfMatrixElements());
    assert(getElements(c2, k10) == SetOfMatrixElements());

    assert(getElements({IH, IH, IH}, Ket(0b00)) == SetOfMatrixElements("0000000000000011" "0000000000110011" "0000000000010001"));
    assert(getElements({IH, IH, IH}, Ket(0b100)) == SetOfMatrixElements());

    auto p1 = split(c1, expectedElements1, SetOfGates("11"));
    assert(p1.left == SetOfGates("10"));
    assert(p1.right == SetOfGates("01"));
    assert(p1.incompatibilities == Incompatibilities());

    auto p2 = split({IH, IH, IH}, getElements({IH, IH, IH}, Ket("01")), SetOfGates("111"));
    assert(p2.incompatibilities == Incompatibilities({SetOfMatrixElements(0b1000000000000001100000000000000000000), SetOfMatrixElements(0b10000000000000000000110000000000000000)}));

    assert(fequal(computeAmplitude({IH, IH, IH}, Ket(0b00)), sqrt(.5)));
    assert(fequal(computeAmplitude({IH, IH, IH}, Ket(0b01)), sqrt(.5)));
    assert(fequal(computeAmplitude({IH, IH, IH}, Ket(0b10)), 0.));
    assert(fequal(computeAmplitude({IH, IH, IH}, Ket(0b11)), 0.));
    assert(fequal(computeAmplitude({IH, CNOT01}, Ket(0b00)), sqrt(.5)));
    assert(fequal(computeAmplitude({IH, CNOT01}, Ket(0b01)), 0.));
    assert(fequal(computeAmplitude({IH, CNOT01}, Ket(0b10)), 0.));
    assert(fequal(computeAmplitude({IH, CNOT01}, Ket(0b11)), sqrt(.5)));
    
    assert(fequal(computeAmplitude({IH, CNOT01, {CNOT, 0, 2}}, Ket(0b000)), sqrt(.5)));
    assert(fequal(computeAmplitude({IH, CNOT01, {CNOT, 0, 2}}, Ket(0b111)), sqrt(.5)));

    // assert(fequal(computeAmplitude({IH, CNOT01, {CNOT, 0, 2}, {}}, Ket(0b111)), sqrt(.5)));
}


int main() {
    test();

    return 0;
}