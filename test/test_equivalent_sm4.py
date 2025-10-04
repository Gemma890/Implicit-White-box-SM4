from core.equivalent_sm4 import get_round_keys

# test get_round_key
masterk = 0x0123456789abcdeffedcba9876543210
rk = get_round_keys(masterk)
print(type(rk[0]))
# for index, x in enumerate(rk):  # 解包元组为 index 和 x
#     # 操作整数值 x
#     hex_str = f"{x & 0xFFFFFFFF:08x}"
#     print(f"Index {index}: {hex_str}")

#   enumerate returns a tuple (index, value),
#   whitch means it needs the index to unpack the list


