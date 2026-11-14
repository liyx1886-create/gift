#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
离线二维码生成器（纯标准库，不联网、不装包）

用法：
    python make_qr.py https://你的网址
    python make_qr.py https://你的网址 qr.png --scale 24
    python make_qr.py https://你的网址 --test     # 生成后自我解码校验

生成的是黑白色 PNG，四周留 4 个模块的静区，可以直接打印或贴到礼物上。
"""

import sys
import zlib
import struct

# ---------------------------------------------------------------- GF(256)
EXP = [0] * 512
LOG = [0] * 256
_x = 1
for _i in range(255):
    EXP[_i] = _x
    LOG[_x] = _i
    _x <<= 1
    if _x & 0x100:
        _x ^= 0x11D
for _i in range(255, 512):
    EXP[_i] = EXP[_i - 255]


def gmul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[LOG[a] + LOG[b]]


def poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            if b:
                out[i + j] ^= gmul(a, b)
    return out


def gen_poly(n):
    g = [1]
    for i in range(n):
        g = poly_mul(g, [1, EXP[i]])
    return g


def rs_ecc(data, n):
    """返回 data 的 n 个纠错码字"""
    g = gen_poly(n)
    buf = list(data) + [0] * n
    for i in range(len(data)):
        c = buf[i]
        if c:
            for j, gj in enumerate(g):
                buf[i + j] ^= gmul(gj, c)
    return buf[len(data):]


# ---------------------------------------------------------- 版本 / 分块表
# 每 4 项为一版：[L, M, Q, H]，每项是 (块数, 总码字, 数据码字) 的整数序列
RS_BLOCK_TABLE = (
    (1, 26, 19), (1, 26, 16), (1, 26, 13), (1, 26, 9),
    (1, 44, 34), (1, 44, 28), (1, 44, 22), (1, 44, 16),
    (1, 70, 55), (1, 70, 44), (2, 35, 17), (2, 35, 13),
    (1, 100, 80), (2, 50, 32), (2, 50, 24), (4, 25, 9),
    (1, 134, 108), (2, 67, 43), (2, 33, 15, 2, 34, 16), (2, 33, 11, 2, 34, 12),
    (2, 86, 68), (4, 43, 27), (4, 43, 19), (4, 43, 15),
    (2, 98, 78), (4, 49, 31), (2, 32, 14, 4, 33, 15), (4, 39, 13, 1, 40, 14),
    (2, 121, 97), (2, 60, 38, 2, 61, 39), (4, 40, 18, 2, 41, 19), (4, 40, 14, 2, 41, 15),
    (2, 146, 116), (3, 58, 36, 2, 59, 37), (4, 36, 16, 4, 37, 17), (4, 36, 12, 4, 37, 13),
    (2, 86, 68, 2, 87, 69), (4, 69, 43, 1, 70, 44), (6, 43, 19, 2, 44, 20), (6, 43, 15, 2, 44, 16),
)

ALIGN_POS = {
    1: [], 2: [6, 18], 3: [6, 22], 4: [6, 26], 5: [6, 30],
    6: [6, 34], 7: [6, 22, 38], 8: [6, 24, 42], 9: [6, 26, 46], 10: [6, 28, 50],
}

ECC_OFFSET = {"L": 0, "M": 1, "Q": 2, "H": 3}
ECC_BITS = {"L": 1, "M": 0, "Q": 3, "H": 2}
MAX_VERSION = 10


def blocks_for(version, ecc):
    row = RS_BLOCK_TABLE[(version - 1) * 4 + ECC_OFFSET[ecc]]
    blocks = []
    for i in range(0, len(row), 3):
        cnt, total, data = row[i], row[i + 1], row[i + 2]
        blocks += [(data, total - data)] * cnt
    return blocks


def char_count_bits(version):
    return 8 if version < 10 else 16


def bit_capacity(version, ecc):
    return sum(d for d, _ in blocks_for(version, ecc)) * 8


# ---------------------------------------------------------------- 编码
def encode_bits(text, version, ecc):
    data = text.encode("utf-8")
    bits = []

    def put(val, n):
        for i in range(n - 1, -1, -1):
            bits.append((val >> i) & 1)

    put(0b0100, 4)                      # 字节模式
    put(len(data), char_count_bits(version))
    for b in data:
        put(b, 8)

    cap = bit_capacity(version, ecc)
    if len(bits) > cap:
        raise ValueError("内容太长，放不进这个版本")
    bits += [0] * min(4, cap - len(bits))          # 终止符
    if len(bits) % 8:
        bits += [0] * (8 - len(bits) % 8)
    pad = [0xEC, 0x11]
    k = 0
    while len(bits) < cap:
        put(pad[k % 2], 8)
        k += 1

    out = []
    for i in range(0, len(bits), 8):
        v = 0
        for b in bits[i:i + 8]:
            v = (v << 1) | b
        out.append(v)
    return out


def interleave(codewords, blocks):
    datas, eccs = [], []
    pos = 0
    for dcount, ecount in blocks:
        chunk = codewords[pos:pos + dcount]
        pos += dcount
        datas.append(chunk)
        eccs.append(rs_ecc(chunk, ecount))

    out = []
    for i in range(max(len(d) for d in datas)):
        for d in datas:
            if i < len(d):
                out.append(d[i])
    for i in range(max(len(e) for e in eccs)):
        for e in eccs:
            if i < len(e):
                out.append(e[i])
    return out


def bch_format(ecc, mask):
    data = (ECC_BITS[ecc] << 3) | mask
    d = data << 10
    g = (1 << 10) | (1 << 8) | (1 << 5) | (1 << 4) | (1 << 2) | (1 << 1) | 1
    while d.bit_length() >= g.bit_length():
        d ^= g << (d.bit_length() - g.bit_length())
    return ((data << 10) | d) ^ 0x5412


def bch_version(version):
    d = version << 12
    g = (1 << 12) | (1 << 11) | (1 << 10) | (1 << 9) | (1 << 8) | (1 << 5) | (1 << 2) | 1
    while d.bit_length() >= g.bit_length():
        d ^= g << (d.bit_length() - g.bit_length())
    return (version << 12) | d


def mask_bit(mask, i, j):
    if mask == 0: return (i + j) % 2 == 0
    if mask == 1: return i % 2 == 0
    if mask == 2: return j % 3 == 0
    if mask == 3: return (i + j) % 3 == 0
    if mask == 4: return (i // 2 + j // 3) % 2 == 0
    if mask == 5: return (i * j) % 2 + (i * j) % 3 == 0
    if mask == 6: return ((i * j) % 2 + (i * j) % 3) % 2 == 0
    return ((i * j) % 3 + (i + j) % 2) % 2 == 0


def blank_matrix(version):
    n = version * 4 + 17
    m = [[None] * n for _ in range(n)]

    def finder(r0, c0):
        for r in range(-1, 8):
            for c in range(-1, 8):
                rr, cc = r0 + r, c0 + c
                if 0 <= rr < n and 0 <= cc < n:
                    inner = (0 <= r <= 6 and c in (0, 6)) or (0 <= c <= 6 and r in (0, 6)) or (2 <= r <= 4 and 2 <= c <= 4)
                    m[rr][cc] = 1 if inner else 0

    finder(0, 0)
    finder(n - 7, 0)
    finder(0, n - 7)

    pos = ALIGN_POS[version]
    for r in pos:
        for c in pos:
            if m[r][c] is not None:
                continue
            for dr in range(-2, 3):
                for dc in range(-2, 3):
                    edge = dr in (-2, 2) or dc in (-2, 2)
                    m[r + dr][c + dc] = 1 if (edge or (dr == 0 and dc == 0)) else 0

    for r in range(8, n - 8):
        if m[r][6] is None:
            m[r][6] = 1 if r % 2 == 0 else 0
    for c in range(8, n - 8):
        if m[6][c] is None:
            m[6][c] = 1 if c % 2 == 0 else 0

    # 预留格式信息区（先置 None 之外的占位，最后统一写）
    for i in range(9):
        if m[8][i] is None: m[8][i] = 0
        if m[i][8] is None: m[i][8] = 0
    for i in range(8):
        if m[8][n - 1 - i] is None: m[8][n - 1 - i] = 0
        if m[n - 1 - i][8] is None: m[n - 1 - i][8] = 0
    m[n - 8][8] = 1

    if version >= 7:
        for i in range(18):
            m[i // 3][i % 3 + n - 8 - 3] = 0
            m[i % 3 + n - 8 - 3][i // 3] = 0
    return m


def place_data(m, codewords, mask):
    n = len(m)
    bits = []
    for cw in codewords:
        for i in range(7, -1, -1):
            bits.append((cw >> i) & 1)

    idx = 0
    inc = -1
    row = n - 1
    for col0 in range(n - 1, 0, -2):
        col = col0 - 1 if col0 <= 6 else col0
        while True:
            for c in (col, col - 1):
                if m[row][c] is None:
                    bit = bits[idx] if idx < len(bits) else 0
                    if mask_bit(mask, row, c):
                        bit ^= 1
                    m[row][c] = bit
                    idx += 1
            row += inc
            if row < 0 or row >= n:
                row -= inc
                inc = -inc
                break
    return m


def write_format(m, version, ecc, mask):
    n = len(m)
    bits = bch_format(ecc, mask)
    for i in range(15):
        b = (bits >> i) & 1
        if i < 6:
            m[i][8] = b
        elif i < 8:
            m[i + 1][8] = b
        else:
            m[n - 15 + i][8] = b
    for i in range(15):
        b = (bits >> i) & 1
        if i < 8:
            m[8][n - 1 - i] = b
        elif i < 9:
            m[8][15 - i - 1 + 1] = b
        else:
            m[8][15 - i - 1] = b
    m[n - 8][8] = 1

    if version >= 7:
        vbits = bch_version(version)
        for i in range(18):
            b = (vbits >> i) & 1
            m[i // 3][i % 3 + n - 8 - 3] = b
            m[i % 3 + n - 8 - 3][i // 3] = b


def penalty(m):
    n = len(m)
    score = 0

    for grid in (m, list(zip(*m))):
        for line in grid:
            run = 1
            for i in range(1, n):
                if line[i] == line[i - 1]:
                    run += 1
                else:
                    if run >= 5:
                        score += 3 + (run - 5)
                    run = 1
            if run >= 5:
                score += 3 + (run - 5)

    for r in range(n - 1):
        for c in range(n - 1):
            v = m[r][c]
            if v == m[r][c + 1] == m[r + 1][c] == m[r + 1][c + 1]:
                score += 3

    p1 = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0]
    p2 = [0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1]
    for grid in (m, list(zip(*m))):
        for line in grid:
            line = list(line)
            for i in range(n - 10):
                seg = line[i:i + 11]
                if seg == p1 or seg == p2:
                    score += 40

    dark = sum(sum(r) for r in m)
    ratio = dark * 100 // (n * n)
    score += (abs(ratio - 50) // 5) * 10
    return score


def make_matrix(text, ecc="M", border=4):
    data_len = len(text.encode("utf-8"))
    version = None
    for v in range(1, MAX_VERSION + 1):
        need = 4 + char_count_bits(v) + data_len * 8
        if need <= bit_capacity(v, ecc):
            version = v
            break
    if version is None:
        raise ValueError("内容太长（超过版本 %d 容量）" % MAX_VERSION)

    codewords = interleave(encode_bits(text, version, ecc), blocks_for(version, ecc))

    best, best_score = None, None
    for mask in range(8):
        m = place_data(blank_matrix(version), codewords, mask)
        write_format(m, version, ecc, mask)
        s = penalty(m)
        if best_score is None or s < best_score:
            best, best_score = m, s

    if not border:
        return best
    n = len(best)
    size = n + border * 2
    out = [[0] * size for _ in range(size)]
    for r in range(n):
        for c in range(n):
            out[r + border][c + border] = best[r][c]
    return out


# ---------------------------------------------------------------- PNG
def write_png(path, matrix, scale=16, dark=(17, 17, 26), light=(255, 255, 255)):
    n = len(matrix)
    size = n * scale

    rows = []
    for r in range(n):
        px = bytearray()
        for c in range(n):
            px += bytes(dark if matrix[r][c] else light) * scale
        row = b"\x00" + bytes(px)
        rows.extend([row] * scale)
    raw = b"".join(rows)

    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data +
                struct.pack(">I", zlib.crc32(tag + data) & 0xFFFFFFFF))

    png = (b"\x89PNG\r\n\x1a\n"
           + chunk(b"IHDR", struct.pack(">IIBBBBB", size, size, 8, 2, 0, 0, 0))
           + chunk(b"IDAT", zlib.compress(raw, 9))
           + chunk(b"IEND", b""))
    with open(path, "wb") as f:
        f.write(png)
    return size


def main(argv):
    rest = argv[1:]
    args, flags = [], []
    scale = 16
    i = 0
    while i < len(rest):
        a = rest[i]
        if a == "--scale" and i + 1 < len(rest):
            scale = int(rest[i + 1])
            i += 2
            continue
        if a.startswith("--scale="):
            scale = int(a.split("=", 1)[1])
        elif a.startswith("--"):
            flags.append(a)
        else:
            args.append(a)
        i += 1

    if not args or args[0] in ("-h", "--help"):
        print(__doc__)
        return 0

    text = args[0]
    out = args[1] if len(args) > 1 else "qr.png"

    matrix = make_matrix(text, "M", border=4)
    size = write_png(out, matrix, scale=scale)
    print("已生成 %s  (%dx%d, 模块数 %d)" % (out, size, size, len(matrix)))
    print("内容: %s" % text)

    if "--test" in flags:
        import verify_qr
        ok, msg = verify_qr.check(matrix, text, "M")
        print(("[自检通过] " if ok else "[自检失败] ") + msg)
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
