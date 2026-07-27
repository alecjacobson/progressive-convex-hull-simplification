import sys, struct, math
N = int(sys.argv[1]); out = sys.argv[2]
# Deterministic Fibonacci-sphere points (all end up on the convex hull).
pts = []
ga = math.pi * (3.0 - math.sqrt(5.0))
for i in range(N):
    y = 1.0 - 2.0*(i+0.5)/N
    r = math.sqrt(max(0.0,1.0-y*y))
    th = ga*i
    pts.append((r*math.cos(th), y, r*math.sin(th)))
# Dummy triangles (ignored for hull; just need valid indices for the reader).
faces = [(i, (i+1)%N, (i+2)%N) for i in range(0, N-2)]
with open(out,'wb') as f:
    hdr = ("ply\nformat binary_little_endian 1.0\n"
           f"element vertex {N}\nproperty double x\nproperty double y\nproperty double z\n"
           f"element face {len(faces)}\nproperty list uchar int vertex_indices\nend_header\n")
    f.write(hdr.encode())
    for p in pts: f.write(struct.pack('<3d', *p))
    for t in faces: f.write(struct.pack('<B3i', 3, *t))
print(f"wrote {out}: {N} verts")
