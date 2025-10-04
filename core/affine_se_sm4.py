import sage.all
from sage.all import *
import collections
import random
from config.constants import CONST
from utils.utilities import (
    int2vector, vector2int,
)

AffineSelfEquivalence = collections.namedtuple('AffineSelfEquivalence',
                                               ['matrix_A', 'cta_A', 'matrix_B', 'cta_B']
                                               )

def get_random_affine_self_equivalence_sm4sbox(number_of_permutations):
	vs = sage.all.VectorSpace(sage.all.GF(2), 8)
	bpr = sage.all.GF(2)
	# sm4_sbox(x)=M(M(x)^(-1)) with M() bing a fixed affine mapping
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

	# 定义GF(2)[x]多项式环，并构造AES的不可约多项式m
	R = PolynomialRing(GF(2), 'x')  # 正确语法
	x = R.gen()
	m = x ** 8 + x ** 7 + x ** 6 + x ** 5 + x ** 4 + x ** 2 + 1  # 确保是8次不可约多项式

	# 使用m作为模数定义GF(2^8)（注意用2**8而非2^8）
	F = GF(2 ** 8, name='x', modulus=m)
	x = F.gen()

	def generate_a_matrix():
		# 随机选择非零元素a
		a = F.random_element()
		while a == F.zero():
			a = F.random_element()

		# 初始化8x8矩阵
		matrix_a = []

		# 遍历基底 {1, x, x^2, ..., x^7}
		for i in range(8):
			element = a * (x ** i)
			coeffs = element.polynomial().list()
			coeffs += [0] * (8 - len(coeffs))  # 补齐8位
			coeffs = coeffs[:8]
			matrix_a.append(coeffs)

		# 构造矩阵（列向量转置为行）
		matrix_a = matrix(GF(2), matrix_a).transpose()
		return a, matrix_a

	def _get_random_affine_selfequivalence():
		a, M_a = generate_a_matrix()
		i = random.randint(0, 7)

		matrix_A = M_inverse * M_a * (H ** (i)) * M_matrix
		cta_A = M_inverse * M_a * (H ** (i)) * M_cta + M_cta_inverse
		matrix_B = M_matrix * (H ** (8 - i)) * M_a * M_inverse
		cta_B = M_matrix * (H ** (8 - i)) * M_a * M_cta_inverse + M_cta

		RandomAffineSelfEquivalence = AffineSelfEquivalence(
			matrix_A,
			cta_A,
			matrix_B,
			cta_B,  # inverse_affine_encoding
		)

		return RandomAffineSelfEquivalence

	random_ase = []
	for _ in range(number_of_permutations):
		random_ase.append(_get_random_affine_selfequivalence())

	return random_ase


if __name__ == '__main__':
	RandomAffineSelfEquivalence = get_random_affine_self_equivalence_sm4sbox()
	y = 75
	print(int2vector(CONST.S_BOX[y], 8))
	tmp = int2vector(y, 8)
	print(tmp)
	tmp = RandomAffineSelfEquivalence.matrix_A * tmp + RandomAffineSelfEquivalence.cta_A
	print(tmp)
	z = vector2int(tmp)
	z = int2vector(CONST.S_BOX[z], 8)
	print(z)
	ans = RandomAffineSelfEquivalence.matrix_B * z + RandomAffineSelfEquivalence.cta_B
	print("ans")
	print(ans)
