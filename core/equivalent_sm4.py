"""Script to generate the implicit (unencoded) and explicit affine layers of SM4 for a fixed key."""
import sys
sys.path.append("..")
import sage.all
from sage.all import *
from functools import partial
from config.constants import CONST
from sage.rings.polynomial.pbori.pbori import BooleanPolynomialVector
from utils.utilities import (
    int2vector, vector2int, matrix2anf,
    BooleanPolynomialRing, compose_affine, compose_anf_fast
)
from core.get_implicit_S_hat_layer import get_4contenated_sbox_implicit_quatratic_anf



def wd_to_byte(wd, bys):
    bys.extend([(wd >> i) & 0xff for i in range(24, -1, -8)])


def bys_to_wd(bys):
    ret = 0
    for i in range(4):
        bits = 24 - i * 8
        ret |= (bys[i] << bits)
    return ret



def get_round_keys(master_key):

    def s_box(wd):
        """
        Look up in sbox table
        :param wd: 32bit word
        :return: 32bit word  ->int
        """
        ret = []
        for i in range(0, 4):
            byte = (wd >> (24 - i * 8)) & 0xff
            row = byte >> 4
            col = byte & 0x0f
            index = (row * 16 + col)
            ret.append(CONST.S_BOX[index])
        return bys_to_wd(ret)

    def rotate_left(wd, bit):
        """
        :param wd: 32bit block needs to be changed
        :param bit: num of bits to rotate_left
        :return:
        """
        return (wd << bit & 0xffffffff) | (wd >> (32 - bit))

    def linear_transformation(wd):
        """
        Linear exchange L
        :param wd: 32bits word
        """
        return wd ^ rotate_left(wd, 2) ^ rotate_left(wd, 10) ^ rotate_left(wd, 18) ^ rotate_left(wd, 24)

    def Tx(k1, k2, k3, ck):
        """
        Key Expansion Exchange Algorithm
        """
        t = s_box(k1 ^ k2 ^ k3 ^ ck)
        return t ^ rotate_left(t, 13) ^ rotate_left(t, 23)

    MK = [(master_key >> (128 - (i + 1) * 32)) & 0xffffffff for i in range(4)]
    # turn 128bit string into 4 words
    keys = [CONST.FK[i] ^ MK[i] for i in range(4)]
    # generate K0~K3
    RK = []
    for i in range(32):
        t = Tx(keys[i + 1], keys[i + 2], keys[i + 3], CONST.CK[i])
        k = keys[i] ^ t
        keys.append(k)
        RK.append(k)

    return RK

def get_implicit_unencoded_affine_layers(
        master_key,
        return_implicit_round_functions=True  # only needed for debugging
):
    rounds = CONST.ROUNDS
    ws = CONST.WORD_SIZE
    rl = CONST.ROUND_SIZE
    bpr = sage.all.GF(2)
    identity_matrix = partial(sage.all.identity_matrix, bpr)
    zero_matrix = partial(sage.all.zero_matrix, bpr)

    ii = identity_matrix(ws)
    zz = zero_matrix(ws, ws)

    A1_matrix = sage.all.block_matrix(bpr, 4, 4, [
        [zz, ii, ii, ii],
        [ii, ii, ii, ii],
        [ii, zz, ii, zz],
        [ii, zz, zz, ii]
    ])

    # Compute L matrix
    def define_rotateleft(num):
        """Construct a rotate left shift by num bits for ws bits matrix over GF(2)"""
        return sage.all.block_matrix(bpr, 2, 2, [
        [zero_matrix(ws - num, num), identity_matrix(ws - num)],
        [identity_matrix(num), zero_matrix(num, ws - num)]
    ])

    rotateleft_0 = define_rotateleft(0)
    rotateleft_2 = define_rotateleft(2)
    rotateleft_10 = define_rotateleft(10)
    rotateleft_18 = define_rotateleft(18)
    rotateleft_24 = define_rotateleft(24)

    L = rotateleft_0 + rotateleft_2 + rotateleft_10 + rotateleft_18 + rotateleft_24

    A2_matrix = sage.all.block_matrix(bpr, 4, 4, [
        [ii, zz, zz, zz],
        [L, ii, zz, zz],
        [zz, ii, ii, zz],
        [zz, ii, zz, ii]
    ])

    A3_matrix = sage.all.block_matrix(bpr, 4, 4, [
        [ii, zz, ii, ii],
        [ii, zz, ii, zz],
        [ii, zz, zz, ii],
        [ii, ii, zz, zz]
    ])

    reverse_order = sage.all.block_matrix(bpr, 4, 4, [
        [zz, zz, zz, ii],
        [zz, zz, ii, zz],
        [zz, ii, zz, zz],
        [ii, zz, zz, zz]
    ])

    sbox_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_BOX, 4)
    bpr_sbox = BooleanPolynomialRing(names=['x%d' % i for i in range(128)] + ['y%d' % i for i in range(128)])
    sbox_layer_anf = [bpr_sbox(str(f)) for f in sbox_layer_anf]

    sbox_inv_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_INV, 4)
    sbox_inv_layer_anf = [bpr_sbox(str(f)) for f in sbox_inv_layer_anf]

    matrix1 = sage.all.block_matrix(bpr, 2, 2, [
        [A1_matrix, zero_matrix(4 * ws, 4 * ws)],
        [zero_matrix(4 * ws, 4 * ws), identity_matrix(4 * ws)]])

    matrix2 = sage.all.block_matrix(bpr, 2, 2, [
        [A2_matrix, zero_matrix(4 * ws, 4 * ws)],
        [zero_matrix(4 * ws, 4 * ws), identity_matrix(4 * ws)]])

    matrix3 = sage.all.block_matrix(bpr, 2, 2, [
        [A3_matrix, zero_matrix(4 * ws, 4 * ws)],
        [zero_matrix(4 * ws, 4 * ws), identity_matrix(4 * ws)]])

    round_keys = get_round_keys(master_key)

    def int2vector2(x, size):
        bits = [(x >> (size - 1 - i)) & 1 for i in range(size)]
        return tuple(bits)

    def bitvectors_to_gf2vector(x, y, z, w):
        return sage.all.vector(bpr, list(int2vector2(x, ws)) + list(int2vector2(y, ws)) + list(int2vector2(z, ws)) + list(int2vector2(w, ws)))

    implicit_round_functions = []

    for i in range(rounds + 1):
        if i == 0:
            # Layer A1
            cta1 = list(bitvectors_to_gf2vector(round_keys[0], round_keys[0], 0, 0)) + [0 for _ in range(4 * ws)]
            anf1 = matrix2anf(matrix1, bpr_sbox, bin_vector=cta1)

            if not return_implicit_round_functions:
                implicit_round_functions.append(anf1)
            else:
                # Layer S /circ A1
                implicit_round_functions.append(compose_anf_fast(sbox_layer_anf, anf1))

        elif 1 <= i < rounds:
            # Layer A2
            cta2 = list(bitvectors_to_gf2vector(0, 0, round_keys[i - 1], round_keys[i - 1])) + [0 for _ in range(4 * ws)]
            anf2 = matrix2anf(matrix2, bpr_sbox, bin_vector=cta2)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf2)
            else:
                # Layer inverse_S /circ A2
                implicit_round_functions.append(compose_anf_fast(sbox_inv_layer_anf, anf2))

            # constant of Layer A3
            c3 = bitvectors_to_gf2vector(round_keys[i - 1], round_keys[i - 1], round_keys[i - 1], 0)
            # constant of Layer A1
            c1 = bitvectors_to_gf2vector(round_keys[i], round_keys[i], 0, 0)
            # Layer A31
            affine = compose_affine(A1_matrix, c1, A3_matrix, c3)
            matrix31 = sage.all.block_matrix(bpr, 2, 2, [
                [affine[0], zero_matrix(4 * ws, 4 * ws)],
                [zero_matrix(4 * ws, 4 * ws), identity_matrix(4 * ws)]])
            cta31 = list(affine[1]) + [0 for _ in range(4 * ws)]
            anf31 = matrix2anf(matrix31, bpr_sbox, bin_vector=cta31)

            # Layer S /circ A31
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf31)
            else:
                # Layer S /circ A31
                implicit_round_functions.append(compose_anf_fast(sbox_layer_anf, anf31))

        else:
            cta2_last = list(bitvectors_to_gf2vector(0, 0, round_keys[i - 1], round_keys[i - 1])) + [0 for _ in range(4 * ws)]
            anf2_last = matrix2anf(matrix2, bpr_sbox, bin_vector=cta2_last)

            cta3 = bitvectors_to_gf2vector(round_keys[i - 1], round_keys[i - 1], round_keys[i - 1], 0)
            affine = compose_affine(reverse_order, 0, A3_matrix, cta3)
            aux = affine[0] ** (-1)
            matrix_last = sage.all.block_matrix(bpr, 2, 2, [
                [identity_matrix(4 * ws), zero_matrix(4 * ws, 4 * ws)],
                [zero_matrix(4 * ws, 4 * ws), aux]])
            c_last = [0 for _ in range(4 * ws)] + list(aux * affine[1])
            anf_last = matrix2anf(matrix_last, bpr_sbox, bin_vector=c_last)
            anf = compose_anf_fast(anf2_last, anf_last)
            if not return_implicit_round_functions:
                implicit_round_functions.append(anf)
            else:
                implicit_round_functions.append(compose_anf_fast(sbox_inv_layer_anf, anf))


    return implicit_round_functions

if __name__ == '__main__':
    master_key = 0x0123456789abcdeffedcba9876543210
    # output_file = "implicit_unencoded_layer.sobj"
    #
    # implicit_affine_layers = get_implicit_unencoded_affine_layers(master_key, True)
    # for i, affine_layer in enumerate(implicit_affine_layers):
    #     # Wrap in tuple because BooleanPolynomialVector can't be pickled.
    #     implicit_affine_layers[i] = tuple(affine_layer)
    #
    # sage.all.save((implicit_affine_layers), output_file , compress=True)

    output_file = "implicit_affine_layer.sobj"

    implicit_affine_layers = get_implicit_unencoded_affine_layers(master_key, False)
    for i, affine_layer in enumerate(implicit_affine_layers):
        # Wrap in tuple because BooleanPolynomialVector can't be pickled.
        implicit_affine_layers[i] = tuple(affine_layer)

    sage.all.save((implicit_affine_layers), output_file, compress=True)