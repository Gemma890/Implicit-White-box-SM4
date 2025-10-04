from sage.all import GF, PolynomialRing, matrix, vector

# 步骤 1：正确定义有限域
F2 = GF(2)
P = PolynomialRing(F2, 'x')
modulus_poly = P('x^8 + x^7 + x^6 + x^5 + x^4 + x^2 + 1')
F256 = GF(2**8, name='a', modulus=modulus_poly)  # 关键修正：2**8

# 步骤 2：定义仿射变换矩阵和常数项
L = matrix(F2, [
	[1, 1, 1, 0, 0, 1, 0, 1],
	[1, 1, 1, 1, 0, 0, 1, 0],
	[0, 1, 1, 1, 1, 0, 0, 1],
	[1, 0, 1, 1, 1, 1, 0, 0],
	[0, 1, 0, 1, 1, 1, 1, 0],
	[0, 0, 1, 0, 1, 1, 1, 1],
	[1, 0, 0, 1, 0, 1, 1, 1],
	[1, 1, 0, 0, 1, 0, 1, 1]
])
c = vector(F2, [1,1,0,0, 1,0,1,1])

# 步骤 3：计算逆变换
L_inv = L.inverse()
c_inv = L_inv * c

def inverse_affine_transform(y_vec):
    y_element = F256.fetch_int(sum(bit << (7 - i) for i, bit in enumerate(y_vec)))
    inv_element = y_element**-1 if y_element != 0 else 0
    inv_vec = inv_element._vector_() if inv_element != 0 else vector(F2, 8)
    return L_inv * inv_vec + c_inv

# 验证测试
test_vec = vector(F2, [0,0,0,0,0,0,1,1])
print("Output:", inverse_affine_transform(test_vec))