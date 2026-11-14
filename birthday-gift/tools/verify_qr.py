#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
二维码自检器：把刚生成的位图重新解码一遍。

它不复用编码器的排布逻辑，独立地：
  1. 从格式信息里读出纠错等级和掩码，并验证 BCH 校验；
  2. 重建功能区（定位图形/时序/校正/格式/版本），把数据区按之字形读回来；
  3. 反交织、算 Reed-Solomon 伴随式（必须全为 0）；
  4. 解析头部并把内容还原成原文与输入比对。

任何一步不对都会失败 —— 这说明生成的码真的扫不出来。
"""

from make_qr import (
    blocks_for, bch_format, mask_bit, EXP, gmul,
)


def _div_rem(value, gen):
    while value.bit_length() >= gen.bit_length():
        value ^= gen << (value.bit_length() - gen.bit_length())
    return value


G15 = (1 << 10) | (1 << 8) | (1 << 5) | (1 << 4) | (1 << 2) | (1 << 1) | 1
G18 = (1 << 12) | (1 << 11) | (1 << 10) | (1 << 9) | (1 << 8) | (1 << 5) | (1 << 2) | 1
ECC_NAME = {0: "M", 1: "L", 2: "H", 3: "Q"}


def reserved(version):
    n = version * 4 + 17
    f = [[False] * n for _ in range(n)]

    # 1) 三个定位图形 + 分隔符
    for r0, c0 in ((0, 0), (0, n - 7), (n - 7, 0)):
        for r in range(r0 - 1, r0 + 8):
            for c in range(c0 - 1, c0 + 8):
                if 0 <= r < n and 0 <= c < n:
                    f[r][c] = True

    # 2) 校正图形（与定位图形重叠的不画）
    from make_qr import ALIGN_POS
    for r in ALIGN_POS[version]:
        for c in ALIGN_POS[version]:
            if f[r][c]:
                continue
            for dr in range(-2, 3):
                for dc in range(-2, 3):
                    f[r + dr][c + dc] = True

    # 3) 时序图形
    for i in range(n):
        f[6][i] = True
        f[i][6] = True

    # 4) 格式信息 + 固定黑模块
    for i in range(9):
        f[8][i] = True
        f[i][8] = True
    for i in range(8):
        f[8][n - 1 - i] = True
        f[n - 1 - i][8] = True
    f[n - 8][8] = True

    # 5) 版本信息
    if version >= 7:
        for i in range(18):
            f[i // 3][i % 3 + n - 8 - 3] = True
            f[i % 3 + n - 8 - 3][i // 3] = True

    return f


def read_format(m):
    n = len(m)
    bits = 0
    for i in range(15):
        if i < 6:
            b = m[i][8]
        elif i < 8:
            b = m[i + 1][8]
        else:
            b = m[n - 15 + i][8]
        if b:
            bits |= 1 << i

    raw = bits ^ 0x5412
    if _div_rem(raw, G15) != 0:
        raise ValueError("格式信息 BCH 校验失败")

    data = raw >> 10
    ecc = data >> 3
    mask = data & 7
    if ecc not in ECC_NAME:
        raise ValueError("纠错等级读出异常")
    return ECC_NAME[ecc], mask, bits


def read_codewords(m, version, mask):
    n = len(m)
    f = reserved(version)
    bits = []
    inc = -1
    row = n - 1
    for col0 in range(n - 1, 0, -2):
        col = col0 - 1 if col0 <= 6 else col0
        while True:
            for c in (col, col - 1):
                if not f[row][c]:
                    b = m[row][c]
                    if mask_bit(mask, row, c):
                        b ^= 1
                    bits.append(b)
            row += inc
            if row < 0 or row >= n:
                row -= inc
                inc = -inc
                break

    out = []
    for i in range(0, len(bits) - 7, 8):
        v = 0
        for b in bits[i:i + 8]:
            v = (v << 1) | b
        out.append(v)
    return out


def deinterleave(cws, blocks):
    datas = [[] for _ in blocks]
    eccs = [[] for _ in blocks]
    idx = 0
    for k in range(max(d for d, _ in blocks)):
        for bi, (dcount, _) in enumerate(blocks):
            if k < dcount:
                datas[bi].append(cws[idx])
                idx += 1
    for k in range(max(e for _, e in blocks)):
        for bi, (_, ecount) in enumerate(blocks):
            if k < ecount:
                eccs[bi].append(cws[idx])
                idx += 1
    return datas, eccs


def syndromes(block, ecc_count):
    """伴随式必须全 0，否则数据/纠错码字对不上"""
    out = []
    for i in range(ecc_count):
        acc = 0
        alpha = EXP[i]
        for c in block:
            acc = gmul(acc, alpha) ^ c
        out.append(acc)
    return out


def check(bitmap, text, ecc="M"):
    total = len(bitmap)
    for border in (4, 0, 3, 2, 5, 6, 1, 7, 8):
        inner = total - 2 * border
        if inner < 21 or (inner - 17) % 4:
            continue
        ok, msg = _check_inner(bitmap, border, text, ecc)
        if ok:
            return ok, msg
    return False, "无法识别（尺寸或格式信息不合法）"


def _check_inner(bitmap, border, text, ecc):
    total = len(bitmap)
    n = total - 2 * border
    version = (n - 17) // 4
    m = [[1 if bitmap[r + border][c + border] else 0 for c in range(n)] for r in range(n)]

    try:
        got_ecc, mask, fmt_bits = read_format(m)
    except ValueError as e:
        return False, str(e)

    if got_ecc != ecc:
        return False, "纠错等级读出为 %s，应为 %s" % (got_ecc, ecc)

    blocks = blocks_for(version, got_ecc)
    cws = read_codewords(m, version, mask)
    need = sum(d + e for d, e in blocks)
    if len(cws) < need:
        return False, "码字不足"

    datas, eccs = deinterleave(cws[:need], blocks)

    for i, ((dcount, ecount), d, e) in enumerate(zip(blocks, datas, eccs)):
        syn = syndromes(d + e, ecount)
        if any(syn):
            return False, "第 %d 块 Reed-Solomon 伴随式非零：%s" % (i, syn[:4])

    stream = []
    for d in datas:
        stream += d

    bits = []
    for cw in stream:
        for i in range(7, -1, -1):
            bits.append((cw >> i) & 1)

    def take(k):
        v = 0
        for b in bits[:k]:
            v = (v << 1) | b
        del bits[:k]
        return v

    mode = take(4)
    if mode != 0b0100:
        return False, "模式位不是字节模式（读出 %s）" % bin(mode)

    length = take(8 if version < 10 else 16)
    payload = bytearray()
    for _ in range(length):
        payload.append(take(8))

    try:
        decoded = payload.decode("utf-8")
    except UnicodeDecodeError:
        return False, "内容不是合法 UTF-8"

    if decoded != text:
        return False, "内容不一致：%r != %r" % (decoded, text)

    return True, "版本 %d / 纠错 %s / 掩码 %d，解码内容完全一致" % (version, got_ecc, mask)


if __name__ == "__main__":
    import sys
    from make_qr import make_matrix
    t = sys.argv[1] if len(sys.argv) > 1 else "https://example.com"
    print(check(make_matrix(t), t))
