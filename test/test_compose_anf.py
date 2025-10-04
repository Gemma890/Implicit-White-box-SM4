import sage.all
from functools import partial

from sympy.abc import alpha

from config.constants import CONST
from core.sm4 import rotate_left
from utils.utilities import (matrix2anf, compose_affine, int2vector, compose_anf_fast, BooleanPolynomialRing)
from core.get_implicit_S_hat_layer import get_4contenated_sbox_implicit_quatratic_anf
from utils.modularaddition import get_implicit_modadd_anf
from core.get_implicit_S_hat_layer import get_4contenated_sbox_implicit_quatratic_anf, get_sbox_implicit_quatratic_matrix,matrix_to_anf_corrected
from sage.rings.polynomial.pbori.pbori import BooleanPolynomialVector
# bpr = sage.all.GF(2)
#
# identity_matrix = partial(sage.all.identity_matrix, bpr)
# zero_matrix = partial(sage.all.zero_matrix, bpr)
# ws = 256
# zz = identity_matrix(ws)
# # zz = zero_matrix(ws, ws)
# anf1 = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_BOX, 4)
# bpr_pmodadd = anf1[0].parent()
# bpr_pmodadd = BooleanPolynomialRing(names=bpr_pmodadd.variable_names(), order="deglex")
# anf1 = [bpr_pmodadd(str(f)) for f in anf1]
#
# anf2 = matrix2anf(zz, bool_poly_ring=bpr_pmodadd)
#
# anf = compose_anf_fast(anf1, anf2)
#
# for i, eq in enumerate(anf):
#         print(f"方程 {i}: {eq} = 0")
#
#
#
def get_implicit_unencoded_affine_layers():
    bpr = sage.all.GF(2)
    ws = 3
    alpha = 0
    beta = 1
    identity_matrix = partial(sage.all.identity_matrix, bpr)
    zero_matrix = partial(sage.all.zero_matrix, bpr)

    ra = sage.all.block_matrix(bpr, 2, 2, [
        [zero_matrix(ws - alpha, alpha), identity_matrix(ws - alpha)],
        [identity_matrix(alpha), zero_matrix(alpha, ws - alpha)]])
    print(ra)
    lb = sage.all.block_matrix(bpr, 2, 2, [
        [zero_matrix(beta, ws - beta), identity_matrix(beta)],
        [identity_matrix(ws - beta), zero_matrix(ws - beta, beta)]])
    assert ra.is_square() and lb.is_square()

    ii = identity_matrix(ws)
    zz = zero_matrix(ws, ws)


    rotateright_identity_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ra, zz],
        [zz, ii]])
    print(rotateright_identity_matrix)

    implicit_pmodadd = get_implicit_modadd_anf(ws, permuted=True, only_x_names=False)
    bpr_pmodadd = implicit_pmodadd[0].parent()
    bpr_pmodadd = BooleanPolynomialRing(names=bpr_pmodadd.variable_names(), order="deglex")
    implicit_pmodadd = [bpr_pmodadd(str(f)) for f in implicit_pmodadd]
    for i, eq in enumerate(implicit_pmodadd):
        print(f"方程 {i}: {eq} = 0")

    implicit_round_functions = []
    matrix = sage.all.block_matrix(bpr, 2, 2, [
            [rotateright_identity_matrix, zero_matrix(2 * ws, 2 * ws)],
            [zero_matrix(2 * ws, 2 * ws), identity_matrix(2 * ws)]])
    print(matrix)
    # anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd)
    anf = matrix2anf(matrix, bool_poly_ring=bpr_pmodadd)
    for i, eq in enumerate(anf):
            print(f"方程 {i}: {eq} = 0")
    print("last")
    implicit_round_functions.append(compose_anf_fast(implicit_pmodadd, anf))
    # print(implicit_round_functions)
    for i, eq in enumerate(implicit_round_functions[0]):
            print(f"方程 {i}: {eq}")

def get_implicit_unencoded_affine_layers1():
    bpr = sage.all.GF(2)
    S = [7, 6, 0, 4, 2, 5, 1, 3]
    ws = 3
    alpha = 1
    beta = 1
    identity_matrix = partial(sage.all.identity_matrix, bpr)
    zero_matrix = partial(sage.all.zero_matrix, bpr)

    ra = sage.all.block_matrix(bpr, 2, 2, [
        [zero_matrix(ws - alpha, alpha), identity_matrix(ws - alpha)],
        [identity_matrix(alpha), zero_matrix(alpha, ws - alpha)]])
    print(ra)
    lb = sage.all.block_matrix(bpr, 2, 2, [
        [zero_matrix(beta, ws - beta), identity_matrix(beta)],
        [identity_matrix(ws - beta), zero_matrix(ws - beta, beta)]])
    assert ra.is_square() and lb.is_square()

    ii = identity_matrix(ws)
    zz = zero_matrix(ws, ws)



    rotateright_identity_matrix = sage.all.block_matrix(bpr, 2, 2, [
        [ra, zz],
        [zz, ii]])
    print(rotateright_identity_matrix)

    m = get_sbox_implicit_quatratic_matrix(3, S)
    s_anf = matrix_to_anf_corrected(m, 3, False, 0)
    bpr_sanf = s_anf[0].parent()
    bpr_sanf = BooleanPolynomialRing(names=bpr_sanf.variable_names(), order="deglex")
    s_anf = [bpr_sanf(str(f)) for f in s_anf]

    for i, eq in enumerate(s_anf):
        print(f"方程 {i}: {eq} = 0")

    implicit_round_functions = []
    # matrix = sage.all.block_matrix(bpr, 2, 2, [
    #         [rotateright_identity_matrix, zero_matrix(2 * ws, 2 * ws)],
    #         [zero_matrix(2 * ws, 2 * ws), identity_matrix(2 * ws)]])
    # print(matrix)
    anf = matrix2anf(rotateright_identity_matrix, bool_poly_ring=bpr_sanf)
    for i, eq in enumerate(anf):
            print(f"方程 {i}: {eq}")
    print("last")
    implicit_round_functions.append(compose_anf_fast(s_anf, anf))
    # print(implicit_round_functions)
    for i, eq in enumerate(implicit_round_functions[0]):
            print(f"方程 {i}: {eq}")

def test_S_compose_Sinv_anf():
    S = [7, 6, 0, 4, 2, 5, 1, 3]
    S_inv = [2, 6, 4, 7, 3, 5, 1, 0]
    s_matrix = get_sbox_implicit_quatratic_matrix(3, S)
    sinv_matrix = get_sbox_implicit_quatratic_matrix(3, S_inv)

    sanf = matrix_to_anf_corrected(s_matrix, 3, False, 0)
    sinv_anf = matrix_to_anf_corrected(sinv_matrix, 3, False, 0)
    for i, eq in enumerate(sanf):
        print(f"方程 {i}: {eq} = 0")
    # bpr_sanf = sanf[0].parent()
    # bpr_sanf = BooleanPolynomialRing(names=bpr_sanf.variable_names(), order="deglex")
    # sanf = [bpr_sanf(str(f)) for f in sanf]
    #
    # bpr_sinvanf = sinv_anf[0].parent()
    # bpr_sinvanf = BooleanPolynomialRing(names=bpr_sinvanf.variable_names(), order="deglex")
    # sinv_anf = [bpr_sinvanf(str(f)) for f in sinv_anf]
    #
    anf = compose_anf_fast(sinv_anf, sanf)
    #
    # # equ = compose_anf_fast(sanf, sinv_anf)
    for i, eq in enumerate(anf):
        print(f"方程 {i}: {eq} = 0")

# test_S_compose_Sinv_anf()




# sbox_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_BOX, 4)

# bpr_g = sbox_layer_anf[0].parent()
# print("左侧函数的布尔环变量数:", bpr_g.n_variables())
# padded_sbox = BooleanPolynomialVector(sbox_layer_anf)
# num_padding = 256 - len(padded_sbox)
# for _ in range(num_padding):
#     padded_sbox.append(bpr_g(0))
# sbox_layer_anf1 = compose_anf_fast(sbox_layer_anf, padded_sbox)
#
#
#
# for i, eq in enumerate(sbox_layer_anf1):
#     print(f"方程 {i}: {eq} = 0")


# print("lulululululuullllllllllllllllllllllll")
# sbox_inv_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_INV, 4)
# for i, eq in enumerate(sbox_inv_layer_anf):
#     print(f"方程 {i}: {eq} = 0")
# bpr_g = sbox_inv_layer_anf[0].parent()
# padded_sbox_inv = BooleanPolynomialVector(sbox_inv_layer_anf)
# num_padding = 256 - len(padded_sbox_inv)
# for _ in range(num_padding):
#     padded_sbox_inv.append(bpr_g(0))
# sbox_layer_inv1 = compose_anf_fast(sbox_inv_layer_anf, padded_sbox_inv)


# equ = compose_anf_fast(sbox_layer_anf, sbox_inv_layer_anf)
# for i, eq in enumerate(equ):
#     print(f"方程 {i}: {eq} = 0")
# get_implicit_unencoded_affine_layers1()

def compose_immodp():
    implicit_pmodadd1 = get_implicit_modadd_anf(3, permuted=True, only_x_names=False)
    bpr_pmodadd1 = implicit_pmodadd1[0].parent()
    bpr_pmodadd1 = BooleanPolynomialRing(names=bpr_pmodadd1.variable_names(), order="deglex")
    implicit_pmodadd1 = [bpr_pmodadd1(str(f)) for f in implicit_pmodadd1]

    implicit_pmodadd2 = get_implicit_modadd_anf(3, permuted=True, only_x_names=False)
    bpr_pmodadd2 = implicit_pmodadd2[0].parent()
    bpr_pmodadd2 = BooleanPolynomialRing(names=bpr_pmodadd1.variable_names(), order="deglex")
    implicit_pmodadd2 = [bpr_pmodadd1(str(f)) for f in implicit_pmodadd2]

    anf = compose_anf_fast(implicit_pmodadd2, implicit_pmodadd1)
    for i, eq in enumerate(anf):
        print(f"方程 {i}: {eq} = 0")

compose_immodp()