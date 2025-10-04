import numpy as np

def build_shift_matrix(shift, size=32):
    M = np.zeros((size, size), dtype=int)
    for i in range(size):
        j = (i + shift) % size  # 正确计算循环左移
        M[i][j] = 1
    return M

# 构造变换矩阵
I = np.eye(32, dtype=int)
M_shift_2 = build_shift_matrix(2)
M_shift_10 = build_shift_matrix(10)
M_shift_18 = build_shift_matrix(18)
M_shift_24 = build_shift_matrix(24)

M_f = (I + M_shift_2 + M_shift_10 + M_shift_18 + M_shift_24) % 2

# 生成随机32位向量
t = np.random.randint(0, 2, 32)

# 矩阵乘法结果
result_matrix = (M_f @ t) % 2

# 手动计算变换结果
def manual_transform(t):
    rotated_2 = np.roll(t, -2)  # 循环左移2位
    rotated_10 = np.roll(t, -10)
    rotated_18 = np.roll(t, -18)
    rotated_24 = np.roll(t, -24)
    return (t + rotated_2 + rotated_10 + rotated_18 + rotated_24) % 2

result_manual = manual_transform(t)

# 断言并打印结果
try:
    assert np.array_equal(result_matrix, result_manual)
except AssertionError:
    print("断言失败：矩阵变换与手动计算结果不一致！")
else:
    print("通过")  # 断言成功时输出