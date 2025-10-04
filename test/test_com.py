from sage.all import *

# 定义GF(2)上的多项式环（无需特殊项序）
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

# 第二次S盒方程 (y → z)
eqs_second = [
    y0*z1 + y2 + z0 + z1,               # 方程0'
    y0*z0 + y1 + z1 + 1,                # 方程1'
    y0 + y1*z0 + y1 + z0 + z1 + z2 + 1, # 方程2'
    y0*z2 + y0 + y1*z1,                 # 方程3'
    y0*z2 + y1*z2 + y1 + y2 + z0 + 1,   # 方程4'
    y0*z2 + y2*z0 + z0 + z2,            # 方程5'
    y0 + y2*z1 + y2 + z0 + z2,          # 方程6'
    y0*z2 + y0 + y1 + y2*z2 + z1 + 1    # 方程7'
]

# 合并所有方程构造理想
I = R.ideal(eqs_first + eqs_second)

# 直接消除中间变量y0,y1,y2
I_elim = I.elimination_ideal([y0, y1, y2])

# 提取结果方程
composed_equations = list(I_elim.gens())

# 打印结果
print("复合隐函数方程组（GF(2)）:")
for i, eq in enumerate(composed_equations):
    print(f"方程 {i}: {eq} = 0")