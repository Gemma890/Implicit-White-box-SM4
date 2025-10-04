import sys
sys.path.append("..")
import sage.all
from sage.all import *
from config.constants import CONST
from sage.rings.polynomial.pbori.pbori import BooleanPolynomialVector
from utils.utilities import BooleanPolynomialRing
from core.get_implicit_S_hat_layer import get_4contenated_sbox_implicit_quatratic_anf
from core.implicit_wbsm4_with_affine_encodings import get_implicit_encoded_round_funcions

if __name__ == '__main__':
	input_file = "implicit_affine_layer.sobj"
	output_file = "implicit_encoded_layer.sobj"
	sbox_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_BOX, 4)
	# bpr_sbox = BooleanPolynomialRing(names=['x%d' % i for i in range(128)] + ['y%d' % i for i in range(128)])
	# sbox_layer_anf = [bpr_sbox(str(f)) for f in sbox_layer_anf]

	sbox_inv_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_INV, 4)
	# sbox_inv_layer_anf = [bpr_sbox(str(f)) for f in sbox_inv_layer_anf]
	implicit_affine_layers = sage.all.load(input_file, compress=True)
	implicit_encoded_round_functions, explicit_extin_anf, explicit_extout_anf = (
		get_implicit_encoded_round_funcions(
			implicit_affine_layers, sbox_layer_anf, sbox_inv_layer_anf
		)
	)

	sage.all.save((implicit_encoded_round_functions, explicit_extin_anf, explicit_extout_anf),
	              output_file, compress=True)