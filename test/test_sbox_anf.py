from config.constants import CONST
from core.get_implicit_S_hat_layer import (
        get_4contenated_sbox_implicit_quatratic_anf,
        get_sbox_implicit_quatratic_matrix, matrix_to_anf_corrected,
        )



# test_equ = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_BOX, 4)
# for i, eq in enumerate(test_equ):
#         print(f"方程 {i}: {eq} = 0")
S = [3, 0, 2, 1]
Sinv = [2, 6, 4, 7, 3, 5, 1, 0]
S = [7, 6, 0, 4, 2, 5, 1, 3]
mat = get_sbox_implicit_quatratic_matrix(8, CONST.S_BOX)
test_equ = matrix_to_anf_corrected(mat, 8)
for eq in test_equ:
    print(eq)
