import sys
sys.path.append("..")
import sage.all
import collections
from sage.all import *
from config.constants import CONST
from functools import partial
from utils.utilities import (
	matrix2anf, BooleanPolynomialRing, compose_anf_fast
)

AffineEncoding = collections.namedtuple('AffineEncoding', ['matrix', 'cta', 'bitsize', 'inverse'])

def get_random_affine_permutations(bitsize, number_of_permutations, TRIVIAL_AE=False, bpr=None):
	vs = sage.all.VectorSpace(sage.all.GF(2), bitsize)

	if bpr is None:
		bpr = sage.all.GF(2)

	def _get_affine_encoding():
		if not TRIVIAL_AE:
			# while loop faster than sage.all.random_matrix(..., algorithm="unimodular")
			while True:
				matrix = sage.all.matrix(bpr, bitsize, entries=[vs.random_element() for _ in range(bitsize)])
				if matrix.is_invertible():
					break
			cta = sage.all.vector(bpr, list(vs.random_element()))
			# matrix = sage.all.random_matrix(sage.all.GF(2), nrows=bitsize, ncols=bitsize, algorithm="unimodular")
			# cta = sage.all.random_vector(sage.all.GF(2), bitsize)
		else:
			matrix = sage.all.identity_matrix(bpr, bitsize)
			cta = sage.all.vector(bpr, [0 for _ in range(bitsize)])

		# inverse_matrix = matrix.inverse()
		# inverse_affine_encoding = AffineEncoding(
			# inverse_matrix,
			# inverse_matrix * cta,
			# bitsize,
			# inverse=None
		# )

		affine_encoding = AffineEncoding(
			matrix,
			cta,
			bitsize,
			inverse=None, # inverse_affine_encoding
		)
		return affine_encoding

	affine_encodings = []
	for _ in range(number_of_permutations):
		affine_encodings.append(_get_affine_encoding())

	return affine_encodings

def get_implicit_affine_round_encodings(TRIVIAL_EE = False, TRIVIAL_AE=False):
	if not TRIVIAL_EE and TRIVIAL_AE:
		raise ValueError("using non-trivial external encoding with trivial affine encodings is not supported")

	rounds = 64
	ws = CONST.WORD_SIZE

	bpr = sage.all.GF(2)
	identity_matrix = partial(sage.all.identity_matrix, bpr)
	zero_matrix = partial(sage.all.zero_matrix, bpr)

	names = ["x" + str(i) for i in range(4 * ws)]
	names += ["y" + str(i) for i in range(4 * ws)]
	bpr_sbox = BooleanPolynomialRing(names=names, order="deglex")

	implicit_round_encodings = [None for _ in range(rounds)]

	affine_encodings = get_random_affine_permutations(4 * ws, rounds, TRIVIAL_AE)

	for i in range(rounds):
		if i == 0:
			if TRIVIAL_EE:
				input_ee_matrix = identity_matrix(4 * ws)
				input_ee_cta = [0 for _ in range(4 * ws)]
			else:
				input_ee = get_random_affine_permutations(4 * ws, 1, TRIVIAL_AE)[0]
				input_ee_matrix = input_ee.matrix
				input_ee_cta = list(input_ee.cta)

			matrix = sage.all.block_matrix(bpr, 2, 2, [
				[input_ee_matrix, zero_matrix(4 * ws, 4 * ws)],
				[zero_matrix(4 * ws, 4 * ws), affine_encodings[i].matrix]])
			cta = input_ee_cta + list(affine_encodings[i].cta)
			implicit_round_encodings[i] = matrix2anf(matrix, bool_poly_ring=bpr_sbox, bin_vector=cta)

		elif 1 <= i < rounds - 1:
			matrix = sage.all.block_matrix(bpr, 2, 2, [
				[affine_encodings[i - 1].matrix, zero_matrix(4 * ws, 4 * ws)],
				[zero_matrix(4 * ws, 4 * ws), affine_encodings[i].matrix]])
			cta = list(affine_encodings[i - 1].cta) + list(affine_encodings[i].cta)
			implicit_round_encodings[i] = matrix2anf(matrix, bool_poly_ring=bpr_sbox, bin_vector=cta)

		else:
			assert i == rounds - 1
			if TRIVIAL_EE:
				output_ee_matrix = identity_matrix(4 * ws)
				output_ee_cta = [0 for _ in range(4 * ws)]
			else:
				output_ee = get_random_affine_permutations(4 * ws, 1, TRIVIAL_AE)[0]
				output_ee_matrix = output_ee.matrix
				output_ee_cta = list(output_ee.cta)

			matrix = sage.all.block_matrix(bpr, 2, 2, [
				[affine_encodings[i - 1].matrix, zero_matrix(4 * ws, 4 * ws)],
				[zero_matrix(4 * ws, 4 * ws), output_ee_matrix]])
			cta = list(affine_encodings[i - 1].cta) + output_ee_cta
			implicit_round_encodings[i] = matrix2anf(matrix, bool_poly_ring=bpr_sbox, bin_vector=cta)

	if TRIVIAL_EE:
		explicit_extin_anf = bpr_sbox.gens()[:4 * ws]
		explicit_extout_anf = bpr_sbox.gens()[:4 * ws]

	else:
		aux_matrix = input_ee.matrix.inverse()
		explicit_extin_anf = matrix2anf(aux_matrix, bool_poly_ring=bpr_sbox, bin_vector=aux_matrix * input_ee.cta)
		explicit_extout_anf = matrix2anf(output_ee.matrix, bool_poly_ring=bpr_sbox, bin_vector=output_ee.cta)

	bpr_x = BooleanPolynomialRing(names=bpr_sbox.variable_names()[:4 * ws], order="deglex")
	explicit_extin_anf = [bpr_x(str(f)) for f in explicit_extin_anf]
	explicit_extout_anf = [bpr_x(str(f)) for f in explicit_extout_anf]

	return implicit_round_encodings, explicit_extin_anf, explicit_extout_anf

def get_graph_automorphisms(TRIVIAL_GA = False):
	ws = CONST.WORD_SIZE
	rounds = 64
	bpr = sage.all.GF(2)
	identity_matrix = partial(sage.all.identity_matrix, bpr)
	zero_matrix = partial(sage.all.zero_matrix, bpr)
	names = ["x" + str(i) for i in range(4 * ws)] + ["y" + str(i) for i in range(4 * ws)]
	bpr_sbox = BooleanPolynomialRing(names=names, order="deglex")

	if TRIVIAL_GA == True:
		identity_matrix = partial(sage.all.identity_matrix, bpr_sbox)

		list_graph_automorphisms = []
		for i in range(rounds):
			list_graph_automorphisms.append(matrix2anf(identity_matrix(8 * ws), bool_poly_ring=bpr_sbox))

		return list_graph_automorphisms

	else:
		from affine_se_sm4 import get_random_affine_self_equivalence_sm4sbox as ase
		graph_automorphisms = [None for _ in range(rounds + 1)]
		for i in range(rounds):
			random_ase = ase(4)
			random_affine = get_random_affine_permutations(3 * ws, 1)
			matrix_a_tmp = sage.all.block_matrix(bpr, 4, 4, [
				[random_ase[0].matrix_A, zero_matrix(8, 8), zero_matrix(8, 8), zero_matrix(8, 8)],
				[zero_matrix(8, 8), random_ase[1].matrix_A, zero_matrix(8, 8), zero_matrix(8, 8)],
				[zero_matrix(8, 8), zero_matrix(8, 8), random_ase[2].matrix_A, zero_matrix(8, 8)],
				[zero_matrix(8, 8), zero_matrix(8, 8), zero_matrix(8, 8), random_ase[3].matrix_A]])
			cta_a_tmp = list(random_ase[0].cta_A) + list(random_ase[1].cta_A) + list(random_ase[2].cta_A) + list(random_ase[3].cta_A)
			A = sage.all.block_matrix(bpr, 2, 2, [
				[matrix_a_tmp, zero_matrix(ws, 3 * ws)],
				[zero_matrix(3 * ws, ws), random_affine[0].matrix]])
			cta_A = cta_a_tmp + list(random_affine[0].cta)

			matrix_b_tmp = sage.all.block_matrix(bpr, 4, 4, [
				[random_ase[0].matrix_B, zero_matrix(8, 8), zero_matrix(8, 8), zero_matrix(8, 8)],
				[zero_matrix(8, 8), random_ase[1].matrix_B, zero_matrix(8, 8), zero_matrix(8, 8)],
				[zero_matrix(8, 8), zero_matrix(8, 8), random_ase[2].matrix_B, zero_matrix(8, 8)],
				[zero_matrix(8, 8), zero_matrix(8, 8), zero_matrix(8, 8), random_ase[3].matrix_B]])
			cta_b_tmp = list(random_ase[0].cta_B) + list(random_ase[1].cta_B) + list(random_ase[2].cta_B) + list(random_ase[3].cta_B)
			matrix_b_tmp_inverse = matrix_b_tmp ** (-1)
			# cta_b_tmp = vector(GF(2), cta_b_tmp)
			cta_b_tmp_inverse = matrix_b_tmp_inverse * vector(GF(2), cta_b_tmp)
			B_inverse = sage.all.block_matrix(bpr, 2, 2, [
				[matrix_b_tmp_inverse, zero_matrix(ws, 3 * ws)],
				[zero_matrix(3 * ws, ws), random_affine[0].matrix]])
			cta_B_inverse = list(cta_b_tmp_inverse) + list(random_affine[0].cta)

			matrix = sage.all.block_matrix(bpr, 2, 2, [
				[A, zero_matrix(4 * ws, 4 * ws)],
				[zero_matrix(4 * ws, 4 * ws), B_inverse]])
			cta = list(cta_A) + list(cta_B_inverse)
			graph_automorphisms[i] = matrix2anf(matrix, bool_poly_ring=bpr_sbox, bin_vector=cta)

		return graph_automorphisms

def get_implicit_encoded_round_funcions(implicit_affine_layers,
                                        sbox_layer_anf, sbox_inv_layer_anf, TRIVIAL_AE=False):
	rounds = 64
	ws = CONST.WORD_SIZE
	assert rounds == len(implicit_affine_layers)
	bpr_sbox = BooleanPolynomialRing(
		names=['x%d' % i for i in range(128)] + ['y%d' % i for i in range(128)]
	)
	sbox_layer_anf = [bpr_sbox(str(f)) for f in sbox_layer_anf]
	sbox_inv_layer_anf = [bpr_sbox(str(f)) for f in sbox_inv_layer_anf]
	implicit_round_encodings, explicit_extin_anf, explicit_extout_anf = get_implicit_affine_round_encodings()
	graph_automorphisms = get_graph_automorphisms()
	left_permutations = get_random_affine_permutations(188, rounds, TRIVIAL_AE, bpr=bpr_sbox)

	implicit_round_functions = []
	list_degs = []

	for i in range(rounds):
		if(i % 2 == 0):
			anf = compose_anf_fast(sbox_layer_anf, graph_automorphisms[i])
			anf = compose_anf_fast(anf, implicit_affine_layers[i])
			anf = compose_anf_fast(anf, implicit_round_encodings[i])
			anf = list(left_permutations[i].matrix * sage.all.vector(bpr_sbox, anf))

			degs = [f.degree() for f in anf]
			assert max(degs) == 2
			list_degs.append(degs)

			implicit_round_functions.append(anf)

		else:
			anf = compose_anf_fast(sbox_inv_layer_anf, graph_automorphisms[i])
			anf = compose_anf_fast(anf, implicit_affine_layers[i])
			anf = compose_anf_fast(anf, implicit_round_encodings[i])
			anf = list(left_permutations[i].matrix * sage.all.vector(bpr_sbox, anf))

			degs = [f.degree() for f in anf]
			assert max(degs) == 2
			list_degs.append(degs)

			implicit_round_functions.append(anf)

	return implicit_round_functions, explicit_extin_anf, explicit_extout_anf