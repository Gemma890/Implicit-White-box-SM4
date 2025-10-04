from sage.all import *
from config.constants import CONST
import re
from sage.sat.converters.polybori import CNFEncoder
from utils.modularaddition import get_implicit_modadd_anf
from core.get_implicit_S_hat_layer import get_4contenated_sbox_implicit_quatratic_anf, get_sbox_implicit_quatratic_matrix,matrix_to_anf_corrected
from sage.rings.polynomial.pbori.pbori import BooleanPolynomialVector
from sage.sat.converters.polybori import CNFEncoder
from sage.sat.solvers.cryptominisat import CryptoMiniSat


def split_linear_nonlinear(eqs):
	"""分离线性和非线性方程"""
	linear, nonlinear = [], []
	for eq in eqs:
		if eq.degree() <= 1:
			linear.append(eq)
		else:
			nonlinear.append(eq)
	return linear, nonlinear


def gaussian_elimination(linear_eqs, BPR):
	"""高斯消元处理线性方程"""
	if not linear_eqs:
		return []
	variables = BPR.gens()
	matrix_rows = []
	for eq in linear_eqs:
		row = [eq.monomial_coefficient(var) for var in variables] + [eq.constant_coefficient()]
		matrix_rows.append(row)
	A = matrix(GF(2), matrix_rows)
	A_echelon = A.echelon_form()
	reduced_eqs = []
	for row in A_echelon:
		lhs = sum(c * var for c, var in zip(row[:-1], variables) if c != 0)
		reduced_eqs.append(lhs + row[-1])
	return reduced_eqs


def optimized_compose_anf_equations(eqs1, eqs2, input_prefix="x", output_prefix="y", intermediate_prefix="z"):
	# 提取输入变量
	input_vars = sorted([var for var in eqs1[0].parent().variable_names() if var.startswith(input_prefix)],
	                    key=lambda x: int(x[len(input_prefix):]))
	num_vars = len(input_vars)
	intermediate_vars = [f"{intermediate_prefix}{i}" for i in range(num_vars)]
	output_vars = [f"{output_prefix}{i}" for i in range(num_vars)]

	# 创建多项式环（包含所有变量）
	all_vars = input_vars + intermediate_vars + output_vars
	R = PolynomialRing(GF(2), names=all_vars)

	# 获取中间变量的环对象列表
	intermediate_vars_in_ring = [R.gen(R.variable_names().index(var)) for var in intermediate_vars]

	# 快速变量替换函数（使用环对象）
	def rename_vars(eqs, src_prefix, tgt_prefix):
		var_map = {
			R.gen(R.variable_names().index(f"{src_prefix}{i}")):
				R.gen(R.variable_names().index(f"{tgt_prefix}{i}"))
			for i in range(num_vars)
		}
		return [eq.subs(var_map) for eq in eqs]

	# 转换方程到标准环并重命名
	eqs1_std = [R(str(eq).replace("^", "**")) for eq in eqs1]
	eqs2_std = [R(str(eq).replace("^", "**")) for eq in eqs2]

	eqs1_renamed = rename_vars(eqs1_std, output_prefix, intermediate_prefix)  # y → z
	eqs2_renamed = rename_vars(eqs2_std, input_prefix, intermediate_prefix)  # x → z

	# 分离线性和非线性方程
	linear_eqs = [eq for eq in eqs1_renamed + eqs2_renamed if eq.degree() <= 1]
	nonlinear_eqs = [eq for eq in eqs1_renamed + eqs2_renamed if eq.degree() > 1]

	# 快速高斯消元
	if linear_eqs:
		A = matrix(GF(2), [
			[eq.monomial_coefficient(var) for var in R.gens()] + [eq.constant_coefficient()]
			for eq in linear_eqs
		])
		A_echelon = A.echelon_form()
		linear_eqs = [sum(c * var for c, var in zip(row, R.gens()) if c) + row[-1] for row in A_echelon]

	# 消元中间变量
	I = R.ideal(linear_eqs + nonlinear_eqs)
	try:
		I_elim = I.elimination_ideal(intermediate_vars_in_ring, algorithm='libsingular:eliminate')
	except:
		I_elim = I.elimination_ideal(intermediate_vars_in_ring)

	return list(I_elim.gens())


# 示例调用
if __name__ == "__main__":
	# 示例矩阵生成（假设P1和P2为两个8x23矩阵）
	ws = 8  # 3比特示例
	S = [7, 6, 0, 4, 2, 5, 1, 3]

	# 生成第一组方程（输入x0,x1,x2，输出z0,z1,z2）
	P1 = get_sbox_implicit_quatratic_matrix(ws, CONST.S_BOX)
	eqs1 = matrix_to_anf_corrected(P1, ws=ws, only_x_names=False, start_index=0)
	#
	# # 生成第二组方程（输入z0,z1,z2，输出y0,y1,y2）
	# P2 = matrix(GF(2), 8, 23)  # 用实际矩阵替换
	eqs2 = matrix_to_anf_corrected(P1, ws=ws, only_x_names=False, start_index=0)

	# 复合方程
	composed_eqs = optimized_compose_anf_equations(eqs1, eqs2)

	# 打印结果
	print("复合后的方程（输入x0,x1,x2，输出y0,y1,y2）:")
	for i, eq in enumerate(composed_eqs):
		print(f"方程 {i}: {eq} = 0")