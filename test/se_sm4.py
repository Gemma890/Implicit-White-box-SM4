from sage.all import *


def compute_inverse_vector(vec, modulus_str):
	"""计算 GF(2^8) 向量的乘法逆元（处理零向量）"""
	# 步骤 1：构造基域 GF(2) 和多项式环
	F_base = GF(2)
	P = PolynomialRing(F_base, "x")

	# 步骤 2：将字符串多项式转换为 SageMath 多项式对象
	modulus_poly = P(modulus_str)  # 关键修复：显式转换

	# 步骤 3：构造有限域 GF(2^8)
	F = GF(2 ** 8, name="a", modulus=modulus_poly)

	# 向量转元素
	integer = sum(bit << (7 - i) for i, bit in enumerate(vec))
	element = F.fetch_int(integer)

	# 计算逆元（处理零）
	if element == 0:
		inv_element = element
	else:
		inv_element = element ** -1

	# 逆元素转向量
	inv_integer = inv_element.integer_representation()
	return [(inv_integer >> (7 - i)) & 1 for i in range(8)]


# 示例使用 AES 不可约多项式
modulus_str = "x^8 + x^7 + x^6 + x^5 + x^4 + x^2 + 1"  # 保持多项式字符串
# modulus_str = "x^8 + x^4 + x^3 + x + 1"
input_vec = [0, 0, 0, 1, 1, 0, 0, 0]  # 对应 a^7


# 计算逆向量
output_vec = compute_inverse_vector(input_vec, modulus_str)
print("Input vector:", input_vec)
print("Inverse vector:", output_vec)