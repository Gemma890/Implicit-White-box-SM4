import collections
import random
from sage.all import *
from config.constants import CONST
from sage.crypto.sboxes import SMS4, AES
from utils.equivalence import (
	get_all_self_ae, get_number_self_ae, check_self_ae_anf,
	get_all_self_le, are_linear_equivalent_lut, are_affine_equivalent_lut
)

AffineSelfEquivalence = collections.namedtuple('AffineSelfEquivalence', ['matrix', 'cta', 'AffineSelfEquivalencePartner'])

def get_random_ae_sm4(bitsize, number_of_permutations, SM4BOX=False, bpr=None):
	vs = sage.all.VectorSpace(sage.all.GF(2), 8)

	if bpr is None:
		bpr = sage.all.GF(2)

	def get_affine_self_equivalence():
		if not SM4BOX:
			# while loop faster than sage.all.random_matrix(..., algorithm="unimodular")
			while True:
				matrix = sage.all.matrix(bpr, bitsize, entries=[vs.random_element() for _ in range(8)])
				if matrix.is_invertible():
					break
			cta = sage.all.vector(bpr, list(vs.random_element()))
			inverse_matrix = matrix.inverse()
			AffineSelfEquivalenceB = AffineSelfEquivalence(
				inverse_matrix,
				inverse_matrix * cta,
				AffineSelfEquivalencePartner=None
			)
			AffineSelfEquivalenceA = AffineSelfEquivalence(
				matrix,
				cta,
				AffineSelfEquivalencePartner = None
			)
			return AffineSelfEquivalenceA, AffineSelfEquivalenceB

		else:
			M_matrix = sage.all.matrix(
				bpr,
				[
					[1, 1, 1, 0, 0, 1, 0, 1],
					[1, 1, 1, 1, 0, 0, 1, 0],
					[0, 1, 1, 1, 1, 0, 0, 1],
					[1, 0, 1, 1, 1, 1, 0, 0],
					[0, 1, 0, 1, 1, 1, 1, 0],
					[0, 0, 1, 0, 1, 1, 1, 1],
					[1, 0, 0, 1, 0, 1, 1, 1],
					[1, 1, 0, 0, 1, 0, 1, 1]
				]
			)

			M_cta = sage.all.vector(bpr, list([1, 1, 0, 0, 1, 0, 1, 1]))

			M_inverse = M_matrix.inverse()
			M_cta_inverse = M_inverse * M_cta

			H = sage.all.matrix(
				bpr,
				[
					[1, 0, 0, 0, 1, 0, 0, 0],
					[0, 0, 0, 0, 0, 1, 0, 1],
					[0, 1, 0, 0, 1, 1, 0, 0],
					[0, 0, 0, 0, 0, 1, 1, 1],
					[0, 0, 1, 0, 1, 1, 1, 0],
					[0, 0, 0, 0, 1, 1, 1, 0],
					[0, 0, 0, 1, 1, 0, 1, 0],
					[0, 0, 0, 0, 1, 0, 1, 0]
				]
			)

			H_inverse = H.inverse()

			R = PolynomialRing(GF(2), 'x')  # 使用标准Python语法
			x = R.gen()
			m = x ** 8 + x ** 7 + x ** 6 + x ** 5 + x ** 4 + x ** 2 + x  # 确保m是8次多项式

			# 使用m作为模数定义GF(2^8)
			F = GF(2 ** 8, name='x', modulus=m)
			x = F.gen()  # 获取域的生成元x

			def generate_a_matrix():
				# 随机选择一个非零元素a
				a = F.random_element()
				while a == F.zero():
					a = F.random_element()

				# 初始化8x8矩阵
				matrix_a = []

				# 遍历基底 {1, x, x^2, ..., x^7}
				for i in range(8):
					# 计算 a * x^i
					element = a * (x ** i)
					# 将结果转换为系数列表（最低位在左）
					coeffs = element.polynomial().coefficients(sparse=False)
					# 补齐8位系数（高位补零）
					coeffs += [0] * (8 - len(coeffs))
					coeffs = coeffs[:8]
					# 添加到矩阵的列
					matrix_a.append(coeffs)

				# 转换为矩阵（列向量转置为行）
				matrix_a = matrix(GF(2), matrix_a).transpose()
				return matrix_a

			matrix_a = generate_a_matrix()
			matrix_a_inverse = matrix_a.inverse()
			random_i = random.randint(0, 7)
			matrix_as_a = M_inverse * matrix_a * (H ** random_i) * M_matrix
			cta_as_a = M_inverse * matrix_a * (H ** random_i) * M_cta + M_cta_inverse

			AffineSelfEquivalenceA = AffineSelfEquivalence(
				matrix_as_a,
				cta_as_a,
				AffineSelfEquivalencePartner=None
			)

			matrix_as_b = M_matrix * (H ** (8 - random_i)) * matrix_a * M_inverse
			cta_as_b = M_matrix * (H ** (8 - random_i)) * matrix_a * M_cta_inverse + M_cta

			AffineSelfEquivalenceB = AffineSelfEquivalence(
				matrix_as_b,
				cta_as_b,
				AffineSelfEquivalencePartner=None
			)
			return AffineSelfEquivalenceA, AffineSelfEquivalenceB








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
		affine_encodings.append(get_affine_self_equivalence())

	return affine_encodings

if __name__ == '__main__':
# 	output_file = "all_sm4sbox_ae.sobj"
# 	# sage.all.block_matrix(1, 2, get_all_self_le(list(AES), return_matrices=True)[-1])
# 	# print(get_all_self_le(list(AES), return_matrices=True)[-1])
# 	print(get_number_self_ae(list(AES)))
# 	# sm4_self_ae =get_all_self_ae(list(SMS4))
# # 	# sage.all.save(sm4_self_ae, output_file, compress=True)
# 	M_cta = sage.all.vector(GF(2), list([0, 1, 1]))
# 	print(M_cta)


# 导入SageMath库
from sage.all import *

# 定义GF(2^8)域，使用AES的不可约多项式 x^8 + x^4 + x^3 + x + 1



# 生成并打印结果
a, M_a = generate_a_matrix()
print(f"随机选择的元素 a = {a} (十六进制: {a.integer_representation():02x})")
print("\n对应的8x8二元矩阵 M_a：")
print(M_a)