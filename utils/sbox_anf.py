"""Auxiliary functions for S_hat_layer."""
from sage.all import *
from utilities import concatenate_anf


S = [7, 6, 0, 4, 2, 5, 1, 3]
S_BOX = [
            0xD6, 0x90, 0xE9, 0xFE, 0xCC, 0xE1, 0x3D, 0xB7, 0x16, 0xB6, 0x14, 0xC2, 0x28, 0xFB, 0x2C, 0x05,
            0x2B, 0x67, 0x9A, 0x76, 0x2A, 0xBE, 0x04, 0xC3, 0xAA, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99,
            0x9C, 0x42, 0x50, 0xF4, 0x91, 0xEF, 0x98, 0x7A, 0x33, 0x54, 0x0B, 0x43, 0xED, 0xCF, 0xAC, 0x62,
            0xE4, 0xB3, 0x1C, 0xA9, 0xC9, 0x08, 0xE8, 0x95, 0x80, 0xDF, 0x94, 0xFA, 0x75, 0x8F, 0x3F, 0xA6,
            0x47, 0x07, 0xA7, 0xFC, 0xF3, 0x73, 0x17, 0xBA, 0x83, 0x59, 0x3C, 0x19, 0xE6, 0x85, 0x4F, 0xA8,
            0x68, 0x6B, 0x81, 0xB2, 0x71, 0x64, 0xDA, 0x8B, 0xF8, 0xEB, 0x0F, 0x4B, 0x70, 0x56, 0x9D, 0x35,
            0x1E, 0x24, 0x0E, 0x5E, 0x63, 0x58, 0xD1, 0xA2, 0x25, 0x22, 0x7C, 0x3B, 0x01, 0x21, 0x78, 0x87,
            0xD4, 0x00, 0x46, 0x57, 0x9F, 0xD3, 0x27, 0x52, 0x4C, 0x36, 0x02, 0xE7, 0xA0, 0xC4, 0xC8, 0x9E,
            0xEA, 0xBF, 0x8A, 0xD2, 0x40, 0xC7, 0x38, 0xB5, 0xA3, 0xF7, 0xF2, 0xCE, 0xF9, 0x61, 0x15, 0xA1,
            0xE0, 0xAE, 0x5D, 0xA4, 0x9B, 0x34, 0x1A, 0x55, 0xAD, 0x93, 0x32, 0x30, 0xF5, 0x8C, 0xB1, 0xE3,
            0x1D, 0xF6, 0xE2, 0x2E, 0x82, 0x66, 0xCA, 0x60, 0xC0, 0x29, 0x23, 0xAB, 0x0D, 0x53, 0x4E, 0x6F,
            0xD5, 0xDB, 0x37, 0x45, 0xDE, 0xFD, 0x8E, 0x2F, 0x03, 0xFF, 0x6A, 0x72, 0x6D, 0x6C, 0x5B, 0x51,
            0x8D, 0x1B, 0xAF, 0x92, 0xBB, 0xDD, 0xBC, 0x7F, 0x11, 0xD9, 0x5C, 0x41, 0x1F, 0x10, 0x5A, 0xD8,
            0x0A, 0xC1, 0x31, 0x88, 0xA5, 0xCD, 0x7B, 0xBD, 0x2D, 0x74, 0xD0, 0x12, 0xB8, 0xE5, 0xB4, 0xB0,
            0x89, 0x69, 0x97, 0x4A, 0x0C, 0x96, 0x77, 0x7E, 0x65, 0xB9, 0xF1, 0x09, 0xC5, 0x6E, 0xC6, 0x84,
            0x18, 0xF0, 0x7D, 0xEC, 0x3A, 0xDC, 0x4D, 0x20, 0x79, 0xEE, 0x5F, 0x3E, 0xD7, 0xCB, 0x39, 0x48
        ]

S_INV = [
            0x71, 0x6C, 0x7A, 0xB8, 0x16, 0x0F, 0x1E, 0x41, 0x35, 0xEB, 0xD0, 0x2A, 0xE4, 0xAC, 0x62, 0x5A,
            0xCD, 0xC8, 0xDB, 0x1A, 0x0A, 0x8E, 0x08, 0x46, 0xF0, 0x4B, 0x96, 0xC1, 0x32, 0xA0, 0x60, 0xCC,
            0xF7, 0x6D, 0x69, 0xAA, 0x61, 0x68, 0x1B, 0x76, 0x0C, 0xA9, 0x14, 0x10, 0x0E, 0xD8, 0xA3, 0xB7,
            0x9B, 0xD2, 0x9A, 0x28, 0x95, 0x5F, 0x79, 0xB2, 0x86, 0xFE, 0xF4, 0x6B, 0x4A, 0x06, 0xFB, 0x3E,
            0x84, 0xCB, 0x21, 0x2B, 0x19, 0xB3, 0x72, 0x40, 0xFF, 0x1C, 0xE3, 0x5B, 0x78, 0xF6, 0xAE, 0x4E,
            0x22, 0xBF, 0x77, 0xAD, 0x29, 0x97, 0x5D, 0x73, 0x65, 0x49, 0xCE, 0xBE, 0xCA, 0x92, 0x63, 0xFA,
            0xA7, 0x8D, 0x2F, 0x64, 0x55, 0xE8, 0xA5, 0x11, 0x50, 0xE1, 0xBA, 0x51, 0xBD, 0xBC, 0xED, 0xAF,
            0x5C, 0x54, 0xBB, 0x45, 0xD9, 0x3C, 0x13, 0xE6, 0x6E, 0xF8, 0x27, 0xD6, 0x6A, 0xF2, 0xE7, 0xC7,
            0x38, 0x52, 0xA4, 0x48, 0xEF, 0x4D, 0x1D, 0x6F, 0xD3, 0xE0, 0x82, 0x57, 0x9D, 0xC0, 0xB6, 0x3D,
            0x01, 0x24, 0xC3, 0x99, 0x3A, 0x37, 0xE5, 0xE2, 0x26, 0x1F, 0x12, 0x94, 0x20, 0x5E, 0x7F, 0x74,
            0x7C, 0x8F, 0x67, 0x88, 0x93, 0xD4, 0x3F, 0x42, 0x4F, 0x33, 0x18, 0xAB, 0x2E, 0x98, 0x91, 0xC2,
            0xDF, 0x9E, 0x53, 0x31, 0xDE, 0x87, 0x09, 0x07, 0xDC, 0xE9, 0x47, 0xC4, 0xC6, 0xD7, 0x15, 0x81,
            0xA8, 0xD1, 0x0B, 0x17, 0x7D, 0xEC, 0xEE, 0x85, 0x7E, 0x34, 0xA6, 0xFD, 0x04, 0xD5, 0x8B, 0x2D,
            0xDA, 0x66, 0x83, 0x75, 0x70, 0xB0, 0x00, 0xFC, 0xCF, 0xC9, 0x56, 0xB1, 0xF5, 0xC5, 0xB4, 0x39,
            0x90, 0x05, 0xA2, 0x9F, 0x30, 0xDD, 0x4C, 0x7B, 0x36, 0x02, 0x80, 0x59, 0xF3, 0x2C, 0xF9, 0x25,
            0xF1, 0xEA, 0x8A, 0x44, 0x23, 0x9C, 0xA1, 0x89, 0x58, 0x8C, 0x3B, 0x0D, 0x43, 0xB5, 0x03, 0xB9,
        ]

FK = [0xA3B1BAC6, 0x56AA3350, 0x677D9197, 0xB27022DC]

CK = [
        0x00070E15, 0x1C232A31, 0x383F464D, 0x545B6269,
        0x70777E85, 0x8C939AA1, 0xA8AFB6BD, 0xC4CBD2D9,
        0xE0E7EEF5, 0xFC030A11, 0x181F262D, 0x343B4249,
        0x50575E65, 0x6C737A81, 0x888F969D, 0xA4ABB2B9,
        0xC0C7CED5, 0xDCE3EAF1, 0xF8FF060D, 0x141B2229,
        0x30373E45, 0x4C535A61, 0x686F767D, 0x848B9299,
        0xA0A7AEB5, 0xBCC3CAD1, 0xD8DFE6ED, 0xF4FB0209,
        0x10171E25, 0x2C333A41, 0x484F565D, 0x646B7279
    ]

def get_sbox_implicit_quatratic_matrix(wordsize, sbox):
    """Return the quadratic anf of the implicit function of the given sbox.

    """
    # 1. get the constant matrix C of the given sbox
    inputs = [[(i >> j) & 1 for j in range(wordsize)] for i in range(2 ** wordsize)]
    outputs = [[(sbox[i] >> j) & 1 for j in range(wordsize)] for i in range(2 ** wordsize)]

    C_cols = []
    for i in range(2 ** wordsize):
        # constant terms and basic terms
        col = [1] + inputs[i] + outputs[i]

        # quatratic cross terms
        cross_terms = [
            inputs[i][m] * outputs[i][n]
            for m in range(wordsize)
            for n in range(wordsize)
        ]

        col += cross_terms
        C_cols.append(col)

    C = Matrix(GF(2), list(zip(*C_cols)))
    c_nrows = C.nrows()
    # C_ref = C.echelon_form()
    # U, T, swaps, rank = C.echelon_form(transformation=True)
    # test data
    # print("\nConstant Matrix C:")
    # print(C)
    # print(C.nrows())

    # 2. get augmented matrix (In|C)
    # simplify matrix C to a echelon one and synchronize these changes to identity matrix on the left.
    def reduced_echelonize_right(I, C):
        M = I.augment(C)
        m = M.nrows()
        n_I = I.ncols()
        n_C = C.ncols()
        total_cols = n_I + n_C
        pivot_cols = []  # Record the columns where the pivots are located (relative to the entire augmented matrix).
        current_row = 0

        # Phase 1: Perform downward elimination to form the row echelon form.
        for col in range(n_I, total_cols):  # Only process the columns of C (i.e., starting from n_I).
            # Find the pivot in the current column.
            pivot = -1
            for r in range(current_row, m):
                if M[r, col] != 0:
                    pivot = r
                    break
            if pivot == -1:
                continue  # This column has no pivot; skip it.

            # Record the columns where the pivots are located.
            pivot_cols.append(col)

            # Swap the rows to the current row position.
            if pivot != current_row:
                M.swap_rows(current_row, pivot)

            # Eliminate the elements below the current column.
            for r in range(current_row + 1, m):
                if M[r, col] != 0:
                    M.add_multiple_of_row(r, current_row, 1)  # addition of GF(2) is oplus
            current_row += 1
            if current_row >= m:
                break  # all rows have been processed

        # Phase 2: Perform upward back substitution to form the reduced form.
        for pivot_row in reversed(range(len(pivot_cols))):
            col = pivot_cols[pivot_row]
            # Eliminate the elements above the pivot column.
            for r in range(pivot_row - 1, -1, -1):
                if M[r, col] != 0:
                    M.add_multiple_of_row(r, pivot_row, 1)

        return M

    # identity_matrix of size (C_rows * C_rows)
    In = identity_matrix(GF(2), c_nrows)
    augmented = reduced_echelonize_right(In, C)

    # Get the result.
    # P = augmented.submatrix(0, 0, c_nrows, c_nrows)
    # C_reduced = augmented.submatrix(0, c_nrows, c_nrows, 2 ** wordsize)

    # I_prime = augmented_matrix.submatrix(0, 0, c_nrows, c_nrows)
    # C_upper = augmented_matrix.submatrix(0, c_nrows, c_nrows, 2 ** wordsize)

    # print("变换后的I矩阵（P）：")
    # print(P)
    # print("\nC的最简行阶梯形：")
    # print(C_reduced)

    r = C.rank()
    AS = augmented.submatrix(r, 0, c_nrows - r, c_nrows)

    return AS

def matrix_to_anf_corrected(P, ws, only_x_names = False, start_index = 0):
    """
        Fixed version ensuring cross terms are generated by GF(2) multiplication:
            - Column structure:
            - Column 0: Constant term
            - Columns 1~ws: x0~x(ws-1)
            - Columns (ws+1)~(2ws): y0~y(ws-1)
            - Columns (2ws+1)~(2ws+ws²): Coefficients for x_i*y_j (lexicographic order)
    """
    ncols = P.ncols()
    expected_ncols = 1 + ws + ws + ws ** 2
    if ncols != expected_ncols:
        raise ValueError(f"Invalid column count. Expected {expected_ncols}, got {ncols}")

    x_vars = [f"x{start_index + i}" for i in range(ws)]
    y_prefix = "x" if only_x_names else "y"
    y_start = start_index + ws if only_x_names else start_index
    y_vars = [f"{y_prefix}{y_start + i}" for i in range(ws)]

    BPR = BooleanPolynomialRing(names=x_vars + y_vars)
    x = BPR.gens()[:ws]
    y = BPR.gens()[ws:ws * 2]

    # Column Index Mapping
    # Boundary indices for column types
    COL_CONSTANT = 0
    COL_X_START = 1
    COL_X_END = COL_X_START + ws
    COL_Y_START = COL_X_END
    COL_Y_END = COL_Y_START + ws
    COL_CROSS_START = COL_Y_END

    # Equation Generation
    equations = []
    for row in P.rows():
        # Handle constant term
        equation = BPR(row[COL_CONSTANT])

        # Linear terms (x_i and y_j)
        for col in range(COL_X_START, COL_X_END):
            if row[col]:
                equation += x[col - COL_X_START]

        for col in range(COL_Y_START, COL_Y_END):
            if row[col]:
                equation += y[col - COL_Y_START]

        # Cross terms (x_i * y_j)
        for k in range(COL_CROSS_START, ncols):
            if row[k]:
                # Calculate cross term indices (i,j)
                cross_idx = k - COL_CROSS_START
                i = cross_idx // ws
                j = cross_idx % ws
                equation += x[i] * y[j]

        equations.append(equation)

    # equations = matrix_to_anf_corrected(AS, wordsize)
    # for i, eq in enumerate(equations):
    #     print(f"方程 {i}: {eq} = 0")

    return equations

def get_4contenated_sbox_implicit_quatratic_anf(sbox, num = 4):
    x_vars = [f"x{i + 32}" for i in range(96)]
    y_vars = [f"y{i + 32}" for i in range(96)]
    BPR = BooleanPolynomialRing(names=x_vars + y_vars)

    sbox_mat = get_sbox_implicit_quatratic_matrix(8, sbox)
    equation_group0 = list(matrix_to_anf_corrected(sbox_mat, 8, False, 0))
    equation_group1 = list(matrix_to_anf_corrected(sbox_mat, 8, False, 8))
    equation_group2 = list(matrix_to_anf_corrected(sbox_mat, 8, False, 16))
    equation_group3 = list(matrix_to_anf_corrected(sbox_mat, 8, False, 24))

    concatenated = concatenate_anf(equation_group1, equation_group0, prefix=None)
    concatenated = concatenate_anf(equation_group2, list(concatenated), prefix=None)
    concatenated = concatenate_anf(equation_group3, list(concatenated), prefix=None)

    # Generate new equations x_i + y_i = 0 for i from 32 to 127
    new_equations = []
    for i in range(32, 128):
        x_i = BPR(f"x{i}")
        y_i = BPR(f"y{i}")
        new_equations.append(x_i + y_i)

    # Merge the new equations into the existing system
    concatenated = concatenate_anf(list(concatenated), new_equations, prefix=None)

    # Print the result
    # for i, eq in enumerate(concatenated):
    #     print(f"方程 {i}: {eq} = 0")

    return concatenated










# test line
# AS = get_sbox_implicit_quatratic_matrix(3, S)
# AS = get_sbox_implicit_quatratic_matrix(8, S_BOX)

# matrix_to_anf_corrected(AS, 3, False, 3)
get_4contenated_sbox_implicit_quatratic_anf(S_BOX, 4)
