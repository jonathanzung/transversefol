import snappy
import spherogram as sp

L=sp.Link(braid_closure=[-4,3,2,1, 5,4,3,2,1,1,2,3,4,5])
M=L.exterior()
print(M.identify())
print(M.volume())
