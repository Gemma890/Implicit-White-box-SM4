from sage.all import *
import re
from config.constants import CONST
from core.get_implicit_S_hat_layer import (
    get_sbox_implicit_quatratic_matrix, matrix_to_anf_corrected
)

def demo_S():
    # 定义GF(2)上的多项式环（无需特殊项序）
    R = PolynomialRing(GF(2), names=['x0', 'x1', 'x2', 'y0', 'y1', 'y2', 'z0', 'z1', 'z2'])
    x0, x1, x2, y0, y1, y2, z0, z1, z2 = R.gens()

    # 第一次S盒方程 (x → y)
    eqs_first = [
        x0 * y1 + x2 + y0 + y1,  # 方程0
        x0 * y0 + x1 + y1 + 1,  # 方程1
        x0 + x1 * y0 + x1 + y0 + y1 + y2 + 1,  # 方程2
        x0 * y2 + x0 + x1 * y1,  # 方程3
        x0 * y2 + x1 * y2 + x1 + x2 + y0 + 1,  # 方程4
        x0 * y2 + x2 * y0 + y0 + y2,  # 方程5
        x0 + x2 * y1 + x2 + y0 + y2,  # 方程6
        x0 * y2 + x0 + x1 + x2 * y2 + y1 + 1  # 方程7
    ]

    # 第二次S盒方程 (y → z)
    eqs_second = [
        y0 * z1 + y2 + z0 + z1,  # 方程0'
        y0 * z0 + y1 + z1 + 1,  # 方程1'
        y0 + y1 * z0 + y1 + z0 + z1 + z2 + 1,  # 方程2'
        y0 * z2 + y0 + y1 * z1,  # 方程3'
        y0 * z2 + y1 * z2 + y1 + y2 + z0 + 1,  # 方程4'
        y0 * z2 + y2 * z0 + z0 + z2,  # 方程5'
        y0 + y2 * z1 + y2 + z0 + z2,  # 方程6'
        y0 * z2 + y0 + y1 + y2 * z2 + z1 + 1  # 方程7'
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

# AS = get_sbox_implicit_quatratic_matrix(8, CONST.S_BOX)
# anf = matrix_to_anf_corrected(AS, 8, False, 0)
# demo_S()
from sage.all import *
import re

def compose_anf_equations(eqs1, eqs2, input_prefix="x", output_prefix="y", intermediate_prefix="z"):
    """
    复合两组ANF方程（输入/输出均为x和y，自动重命名中间变量为z）
    """
    # 提取变量并确定变量数量
    def extract_vars(eqs):
        vars = set()
        for eq in eqs:
            vars.update(str(eq).replace(" ", "").split("+"))
            vars.update(str(eq).replace(" ", "").split("*"))
        filtered_vars = [v for v in vars if re.match(r"^[a-zA-Z]+\d+$", v)]
        return sorted(filtered_vars, key=lambda s: (s[0], int(re.findall(r"\d+", s)[0])))

    input_vars = [v for v in extract_vars(eqs1) if v.startswith(input_prefix)]
    num_vars = len(input_vars)

    # 重命名函数
    def rename_equations(eqs, old_output, new_output):
        var_map = {f"{old_output}{i}": f"{new_output}{i}" for i in range(num_vars)}
        renamed_eqs = []
        for eq in eqs:
            eq_str = str(eq)
            for old, new in var_map.items():
                eq_str = re.sub(rf"\b{old}\b", new, eq_str)
            renamed_eqs.append(eq_str)
        return renamed_eqs

    # 重命名方程变量
    eqs1_renamed = rename_equations(eqs1, output_prefix, intermediate_prefix)
    eqs2_renamed = rename_equations(eqs2, input_prefix, intermediate_prefix)

    # 创建标准多项式环（非BooleanPolynomialRing）
    all_vars = input_vars + [f"{intermediate_prefix}{i}" for i in range(num_vars)] + [f"{output_prefix}{i}" for i in range(num_vars)]
    R = PolynomialRing(GF(2), names=all_vars)
    R.inject_variables()

    # 将方程转换为标准多项式环
    def parse_eqs(eq_strs):
        eqs = []
        for s in eq_strs:
            s = s.replace("^", "**")
            eq = R(s)
            eqs.append(eq)
        return eqs

    # 合并方程并消元
    combined_eqs = parse_eqs(eqs1_renamed + eqs2_renamed)
    I = R.ideal(combined_eqs)
    I_elim = I.elimination_ideal([R(f"{intermediate_prefix}{i}") for i in range(num_vars)])

    # 转换回布尔多项式形式（可选）
    BPR = BooleanPolynomialRing(names=input_vars + [f"{output_prefix}{i}" for i in range(num_vars)])
    final_eqs = [BPR(str(eq)) for eq in I_elim.gens()]
    return final_eqs

def compose_anf_equations1(eqs1, eqs2, input_prefix="x", output_prefix="y", intermediate_prefix="z"):
    """
    复合两组ANF方程（输入/输出均为x和y，自动重命名中间变量为z）
    """
    # import re
    # from sage.rings.polynomial.pbori import BooleanPolynomialRing
    # from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_generic
    from sage.rings.finite_rings.finite_field_constructor import GF

    # 提取变量并确定变量数量
    def extract_vars(eqs):
        vars = set()
        for eq in eqs:
            vars.update(str(eq).replace(" ", "").split("+"))
            vars.update(str(eq).replace(" ", "").split("*"))
        filtered_vars = [v for v in vars if re.match(r"^[a-zA-Z]+\d+$", v)]
        return sorted(filtered_vars, key=lambda s: (s[0], int(re.findall(r"\d+", s)[0])))

    input_vars = [v for v in extract_vars(eqs1) if v.startswith(input_prefix)]
    num_vars = len(input_vars)

    # 重命名函数
    def rename_equations(eqs, old_output, new_output):
        var_map = {f"{old_output}{i}": f"{new_output}{i}" for i in range(num_vars)}
        renamed_eqs = []
        for eq in eqs:
            eq_str = str(eq)
            for old, new in var_map.items():
                eq_str = re.sub(rf"\b{old}\b", new, eq_str)
            renamed_eqs.append(eq_str)
        return renamed_eqs

    # 重命名方程变量
    eqs1_renamed = rename_equations(eqs1, output_prefix, intermediate_prefix)
    eqs2_renamed = rename_equations(eqs2, input_prefix, intermediate_prefix)

    # 创建标准多项式环（非BooleanPolynomialRing）
    all_vars = input_vars + [f"{intermediate_prefix}{i}" for i in range(num_vars)] + [f"{output_prefix}{i}" for i in range(num_vars)]
    R = PolynomialRing(GF(2), names=all_vars)
    R.inject_variables()

    # 将方程转换为标准多项式环
    def parse_eqs(eq_strs):
        eqs = []
        for s in eq_strs:
            s = s.replace("^", "**")
            eq = R(s)
            eqs.append(eq)
        return eqs

    # 合并方程并消元
    combined_eqs = parse_eqs(eqs1_renamed + eqs2_renamed)
    I = R.ideal(combined_eqs)
    I_elim = I.elimination_ideal([R(f"{intermediate_prefix}{i}") for i in range(num_vars)])

    # 转换回布尔多项式形式并过滤零多项式
    BPR = BooleanPolynomialRing(names=input_vars + [f"{output_prefix}{i}" for i in range(num_vars)])
    final_eqs = []
    for eq in I_elim.gens():
        poly = BPR(str(eq))
        if poly != 0:
            final_eqs.append(poly)
    return final_eqs
# 示例用法
if __name__ == "__main__":
    # 示例矩阵生成（假设P1和P2为两个8x23矩阵）
    ws = 3  # 3比特示例
    S = [7, 6, 0, 4, 2, 5, 1, 3]
    Sinv = [2, 6, 4, 7, 3, 5, 1, 0]

    # 生成第一组方程（输入x0,x1,x2，输出z0,z1,z2）
    P1 = get_sbox_implicit_quatratic_matrix(ws, S)
    eqs1 = matrix_to_anf_corrected(P1, ws=ws, only_x_names=False, start_index=0)
    #
    # # 生成第二组方程（输入z0,z1,z2，输出y0,y1,y2）
    # P2 = matrix(GF(2), 8, 23)  # 用实际矩阵替换
    eqs2 = matrix_to_anf_corrected(P1, ws=ws, only_x_names=False, start_index=0)

    # 复合方程
    composed_eqs = compose_anf_equations(eqs1, eqs2)

    # 打印结果
    print("复合后的方程（输入x0,x1,x2，输出y0,y1,y2）:")
    for i, eq in enumerate(composed_eqs):
        print(f"方程 {i}: {eq} = 0")

