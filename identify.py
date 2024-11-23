import snappy

M=snappy.Triangulation("l6a5.tri")
print(M)
M.num_tetrahedra()
snappy.LinkExteriors.identify(M)
