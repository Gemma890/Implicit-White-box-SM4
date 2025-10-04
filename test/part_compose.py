import re
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.ideal import Ideal
from utils.utilities import BooleanPolynomialRing
from sage.all import *
from core.get_implicit_S_hat_layer import (
	get_sbox_implicit_quatratic_matrix, matrix_to_anf_corrected,
    get_4contenated_sbox_implicit_quatratic_anf)
from config.constants import CONST
from sage.rings.ideal import Ideal
def compose_anf_equations_stepwise(eqs1, eqs2, input_prefix="x", output_prefix="y", intermediate_prefix="z"):


	def extract_vars_by_prefix(eqs, prefix):
		vars_set = set()
		for eq in eqs:
			terms = re.split(r"[+*]+", str(eq).replace(" ", ""))
			for term in terms:
				if re.match(rf"^{prefix}\d+$", term):
					vars_set.add(term)
		# 按数字顺序排序变量
		return sorted(vars_set, key=lambda s: int(re.findall(r"\d+", s)[0]))

	# 提取变量信息
	input_vars_eqs1 = extract_vars_by_prefix(eqs1, input_prefix)
	output_vars_eqs1 = extract_vars_by_prefix(eqs1, output_prefix)
	num_outputs_eqs1 = len(output_vars_eqs1)

	input_vars_eqs2 = extract_vars_by_prefix(eqs2, input_prefix)
	remaining_inputs_eqs2 = input_vars_eqs2[num_outputs_eqs1:]  # 未被替换的输入变量

	# 重命名函数（仅替换前num_replace个变量）
	def rename_vars(eqs, old_prefix, new_prefix, num_replace):
		var_map = {f"{old_prefix}{i}": f"{new_prefix}{i}" for i in range(num_replace)}
		renamed = []
		for eq in eqs:
			eq_str = str(eq)
			for old, new in var_map.items():
				eq_str = re.sub(rf"\b{old}\b", new, eq_str)
			renamed.append(eq_str)
		return renamed

	# 连接两个方程组
	eqs1_renamed = rename_vars(eqs1, output_prefix, intermediate_prefix, num_outputs_eqs1)
	eqs2_renamed = rename_vars(eqs2, input_prefix, intermediate_prefix, num_outputs_eqs1)

	# 构建多项式环（包含所有变量）
	all_vars = (
			input_vars_eqs1 +
			remaining_inputs_eqs2 +
			[f"{intermediate_prefix}{i}" for i in range(num_outputs_eqs1)] +
			extract_vars_by_prefix(eqs2, output_prefix)
	)
	R = PolynomialRing(GF(2), all_vars)
	R.inject_variables()

	# 方程转换
	def parse_eqs(eq_strs):
		return [R(eq.replace("^", "**")) for eq in eq_strs]

	combined_eqs = parse_eqs(eqs1_renamed + eqs2_renamed)
	I = R.ideal(combined_eqs)

	# 仅消去中间变量
	intermediate_vars = [R.gen(all_vars.index(f"{intermediate_prefix}{i}")) for i in range(num_outputs_eqs1)]
	I_elim = I.elimination_ideal(intermediate_vars)

	# 转换回布尔多项式环
	output_vars = extract_vars_by_prefix(eqs2, output_prefix)
	BPR = BooleanPolynomialRing(
		names=input_vars_eqs1 + remaining_inputs_eqs2 + output_vars
	)
	final_eqs = []
	for poly in I_elim.gens():
		expr = BPR(str(poly))
		if expr != 0:
			final_eqs.append(expr)

	return final_eqs


def compose_anf_equations_stepwise_optimized(eqs1, eqs2, input_prefix="x", output_prefix="y", intermediate_prefix="z"):
	"""优化内存版本的方程组复合函数"""

	# 使用生成器替代列表的辅助函数
	def extract_vars_generator(eqs, prefix):
		"""生成器方式提取变量"""
		pattern = re.compile(rf"\b{prefix}\d+\b")
		seen = set()
		for eq in eqs:
			for term in re.split(r"[+*]+", str(eq).replace(" ", "")):
				match = pattern.match(term)
				if match and term not in seen:
					seen.add(term)
					yield term
		# 按数字排序返回列表（仅在最终需要时生成）
		return sorted(seen, key=lambda s: int(s[1:]))

	# 使用惰性求值提取变量
	input_vars_eqs1 = list(extract_vars_generator(eqs1, input_prefix))
	output_vars_eqs1 = list(extract_vars_generator(eqs1, output_prefix))
	num_outputs = len(output_vars_eqs1)

	# 动态构建变量映射（避免预先生成大字典）
	def rename_terms(eqs, src_prefix, tgt_prefix, replace_count):
		"""流式处理变量重命名"""
		pattern = re.compile(rf"\b{src_prefix}(\d+)\b")
		return [
			pattern.sub(
				lambda m: f"{tgt_prefix}{m.group(1)}" if int(m.group(1)) < replace_count else m.group(0),
				str(eq)
			) for eq in eqs
		]

	# 分步处理方程组合
	stage1_eqs = rename_terms(eqs1, output_prefix, intermediate_prefix, num_outputs)
	stage2_eqs = rename_terms(eqs2, input_prefix, intermediate_prefix, num_outputs)

	# 动态构建多项式环（延迟变量创建）
	all_vars = (
			input_vars_eqs1 +
			[v for v in extract_vars_generator(eqs2, input_prefix) if int(v[1:]) >= num_outputs] +
			[f"{intermediate_prefix}{i}" for i in range(num_outputs)] +
			list(extract_vars_generator(eqs2, output_prefix))
	)

	# 使用布尔多项式环节省内存
	BPR = BooleanPolynomialRing(names=all_vars)
	BPR.inject_variables(verbose=False)

	# 流式方程解析
	def parse_equation(eq_str):
		"""内存优化的方程解析"""
		return BPR(eq_str.replace("^", "**").replace("*", "&"))

	# 使用生成器处理方程
	combined_eqs = (parse_equation(eq) for eq in stage1_eqs + stage2_eqs)

	# 分块构建理想
	I = Ideal([next(combined_eqs)])
	for eq in combined_eqs:
		I += eq

	# 消去中间变量（逐步处理）
	intermediate_vars = [BPR(f"{intermediate_prefix}{i}") for i in range(num_outputs)]
	I_elim = I.elimination_ideal(intermediate_vars)

	# 直接输出结果（避免二次转换）
	return [eq for eq in I_elim.gens() if eq != 0]

R = PolynomialRing(GF(2), names=['x0','x1','x2','y0','y1','y2','z0','z1','z2'])
x0, x1, x2, y0, y1, y2, z0, z1, z2 = R.gens()

# 第一次S盒方程 (x → y)
eqs_first = [
    x0*y1 + x2 + y0 + y1,               # 方程0
    x0*y0 + x1 + y1 + 1,                # 方程1
    x0 + x1*y0 + x1 + y0 + y1 + y2 + 1, # 方程2
    x0*y2 + x0 + x1*y1,                 # 方程3
    x0*y2 + x1*y2 + x1 + x2 + y0 + 1,   # 方程4
    x0*y2 + x2*y0 + y0 + y2,            # 方程5
    x0 + x2*y1 + x2 + y0 + y2,          # 方程6
    x0*y2 + x0 + x1 + x2*y2 + y1 + 1    # 方程7
]
eqs_inv = [
x0*y1 + x0 + x1 + x2 + y0 + y1 + 1,
x0*y0 + x1 + y1 + 1,
x0 + x1*y0 + x1 + y2,
x0*y2 + x0 + x1*y1 + x2 + y0,
x0 + x1*y2 + x2 + y0 + y2,
x0*y2 + x0 + x2*y0 + x2,
x0*y2 + x2*y1 + x2 + y1 + y2 + 1,
x0*y2 + x0 + x1 + x2*y2 + x2 + y0 + y1 + 1
]
equ_second = [
	x0 + y1 + 1,
	x0 + x1 + y0 + 1,
	x0*y1,
	x0*y0 + x1*y0,
	x0*y0 + x1*y1 + x1,
	x2 + y2
]
sbox_inv_mat = get_sbox_implicit_quatratic_matrix(8, CONST.S_BOX)
sbox_inv_anf = matrix_to_anf_corrected(sbox_inv_mat, 8)
sbox_layer_anf = get_4contenated_sbox_implicit_quatratic_anf(CONST.S_BOX)
bpr_sbox = sbox_layer_anf[0].parent()
bpr_sbox = BooleanPolynomialRing(names=bpr_sbox.variable_names(), order="deglex")
sbox_layer_anf = [bpr_sbox(str(f)) for f in sbox_layer_anf]

final_eqs = compose_anf_equations_stepwise_optimized(sbox_layer_anf[:23], sbox_layer_anf)
# final_eqs = compose_anf_equations_stepwise(sbox_layer_anf[5:9], final_eqs)
for i, eq in enumerate(final_eqs):
	print(f"方程 {i}: {eq} = 0")