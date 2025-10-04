"""Script to test the correctness of implicit equavalent SM4 encryption for a fixed key."""
import sys
sys.path.append("..")
import sage.all
from sage.all import GF, vector
from config.constants import CONST
from utils.functionalequations import find_fixed_vars
from utils.utilities import (
    int2vector, vector2int,
    substitute_variables,
    BooleanPolynomialRing
)

def bitvectors_to_gf2vector(bitvectors, ws):
    parts = []
    for bit in bitvectors:
        vec_part = list(int2vector(bit, ws))
        parts.append(vec_part)

    return sage.all.vector(sage.all.GF(2), parts)

def gf2vector_to_bitvectors(v, ws):
    return vector2int(v[:ws]), vector2int(v[ws:])

def equivalent_SM4_encrypt(implicit_unencoded_round_function):
    ws = CONST.WORD_SIZE
    rounds = 64
    assert len(implicit_unencoded_round_function) == rounds
    bpr_sbox = BooleanPolynomialRing(
        names=['x%d' % i for i in range(128)] + ['y%d' % i for i in range(128)]
    )

    ordered_replacement = []
    assert len(bpr_sbox.gens()) == 8 * ws

    output_vars = bpr_sbox.gens()[4 * ws: 8 * ws]
    ordered_replacement = [None] * (4 * ws) + list(output_vars)

    def eval_round_function(v, round_index):
        ordered_replacement_copy = ordered_replacement[:]
        for i in range(4 * ws):
            ordered_replacement_copy[i] = bpr_sbox(v[i])

        list_outputs = []
        systems_in_round_i = [implicit_unencoded_round_function[round_index]]
        assert len(systems_in_round_i) == 1

        for index_irf, implicit_rf in enumerate(systems_in_round_i):
            system = [
                substitute_variables(bpr_sbox, ordered_replacement_copy, f)
                for f in implicit_rf
            ]

            if not all(f.degree() <= 1 for f in system):  # assuming QUASILINEAR_RP is True
                raise ValueError(
                    f"implicit round function {index_irf} is not quasilinear "
                    f"(has degrees {[f.degree() for f in system]} after fixing the input variables)"
                )

            try:
                fixed_vars, new_equations = find_fixed_vars(
                    system,
                    only_linear=True,
                    initial_r_mode="gauss",
                    repeat_with_r_mode=None,
                    initial_fixed_vars=None,
                    bpr=bpr_sbox, check=False,
                    verbose=False,
                    debug=False,
                    filename=None
                )
            except ValueError as e:
                raise ValueError(f"implicit round function {index_irf} has no solution, found error {e}")
                assert str(e).startswith("found 0 == 1")
                continue

            found_non_cta = any(v not in [0, 1] for v in fixed_vars.values())
            if found_non_cta or len(new_equations) >= 1:
                raise ValueError(f"implicit round function {index_irf} has no unique solution")
                continue

            assert len(new_equations) == 0, f"{fixed_vars}\n{list(new_equations)}"

            sol = [fixed_vars.get(v, None) for v in output_vars]

            if None in sol:
                raise ValueError(f"implicit round function {index_irf} has no unique solution")

            list_outputs.append(tuple(sol))

        assert len(list_outputs) == 1
        v0 = list_outputs[0]

        return sage.all.vector(sage.all.GF(2), v0)

    def eval_implicit_implementation(v):
        assert len(v) == 4 * ws, f"输入向量长度必须为 {4 * ws}"
        v = vector(GF(2), v)  # 确保为 GF(2) 向量
        for i in range(rounds):
            print(f"正在处理轮次 {i}...")
            v = eval_round_function(v, i)

        return v

    return eval_implicit_implementation



if __name__ == '__main__':
    input_file = "implicit_unencoded_layer.sobj"
    implicit_unencoded_round_function = sage.all.load(input_file, compress=True)

    plaintext = 0x0123456789abcdeffedcba9876543210
    # plaintext = 0x00000000000000000000000000000001
    binary_str = format(plaintext, "0128b")

    gf2_vec = vector(GF(2), [int(bit) for bit in binary_str])
    eval_wb = equivalent_SM4_encrypt(implicit_unencoded_round_function)
    ciphertext = eval_wb(gf2_vec)
    print("加密结果:", ciphertext)
    ciphertext = ''.join(str(bit) for bit in ciphertext)
    hex_str = ""
    for i in range(0, len(ciphertext), 8):  # 每8位（1字节）处理
        byte = ciphertext[i:i + 8]
        hex_byte = format(int(byte, 2), "02x")  # 固定2位十六进制
        hex_str += hex_byte + " "

    print(hex_str.strip())



