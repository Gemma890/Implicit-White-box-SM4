from sage.all import *
from config.constants import CONST
from core.get_implicit_S_hat_layer import get_4contenated_sbox_implicit_quatratic_anf

sbox_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_BOX, 4)
bpr_sbox = sbox_layer_anf[0].parent()
bpr_sbox = BooleanPolynomialRing(names=bpr_sbox.variable_names(), order="deglex")
sbox_layer_anf = [bpr_sbox(str(f)) for f in sbox_layer_anf]

sbox_inv_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_INV, 4)
bpr_sbox_inv = sbox_inv_layer_anf[0].parent()
bpr_sbox_inv = BooleanPolynomialRing(names=bpr_sbox_inv.variable_names(), order="deglex")
sbox_inv_layer_anf = [bpr_sbox_inv(str(f)) for f in sbox_inv_layer_anf]

print(bpr_sbox)
print("lalala")
print(bpr_sbox_inv)