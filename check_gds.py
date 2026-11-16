import struct, sys

RECNAMES = {0:"HEADER",1:"BGNLIB",2:"LIBNAME",3:"UNITS",4:"ENDLIB",5:"BGNSTR",
            6:"STRNAME",7:"ENDSTR",8:"BOUNDARY",13:"LAYER",14:"DATATYPE",
            16:"XY",17:"ENDEL"}

def real8(b):
    if b == b'\x00'*8: return 0.0
    sign = -1.0 if b[0] & 0x80 else 1.0
    e = (b[0] & 0x7F) - 64
    mant = int.from_bytes(b[1:8], 'big') / 72057594037927936.0  # 16^14：7字节尾数为56位十六进制小数
    return sign * mant * (16.0 ** e)

path = sys.argv[1]
data = open(path, 'rb').read()
pos = 0
counts = {}
layer_counts = {}
xy_all = []
expect = None
nrec = 0
while pos < len(data):
    if pos + 4 > len(data):
        print("TRUNCATED header at", pos); break
    length = struct.unpack('>H', data[pos:pos+2])[0]
    rt, dt = data[pos+2], data[pos+3]
    if length < 4 or pos + length > len(data):
        print("INVALID length", length, "at", pos); break
    payload = data[pos+4:pos+length]
    name = RECNAMES.get(rt, f"RT{rt}")
    counts[name] = counts.get(name, 0) + 1
    if name == "UNITS":
        u1, u2 = real8(payload[0:8]), real8(payload[8:16])
        print(f"UNITS: {u1} user-unit/dbu, {u2} m/dbu")
    if name == "LIBNAME" or name == "STRNAME":
        print(f"{name}: {payload.rstrip().decode('ascii')}")
    if name == "LAYER":
        layer = struct.unpack('>H', payload)[0]
        layer_counts[layer] = layer_counts.get(layer, 0) + 1
    if name == "XY":
        n = len(payload) // 8
        pts = [struct.unpack('>ii', payload[i*8:(i+1)*8]) for i in range(n)]
        xy_all.extend(pts[:4])  # 前4点即矩形四角(第5点闭合)
    if name == "HEADER":
        print("HEADER version", struct.unpack('>H', payload)[0])
    pos += length
    nrec += 1

print(f"\n共 {nrec} 条记录, 覆盖 {pos}/{len(data)} 字节 ({'完整' if pos==len(data) else '不完整!'})")
print("记录统计:", counts)
print("各图层边界数:", layer_counts)
xs = [p[0] for p in xy_all]; ys = [p[1] for p in xy_all]
print(f"坐标范围: x [{min(xs)}, {max(xs)}] nm, y [{min(ys)}, {max(ys)}] nm")
print(f"即 x [{min(xs)/1000}, {max(xs)/1000}] um, y [{min(ys)/1000}, {max(ys)/1000}] um")
