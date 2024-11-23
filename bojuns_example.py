import snappy
import spherogram
from functools import reduce
import operator


def find_s3_slope(M):
	for ss in M.short_slopes()[0]:
		M.dehn_fill(ss,0)
		G=M.fundamental_group()
		if len(G.generators())==0:
			return ss
	return None

M=snappy.Manifold("o9_19364")
print(M.volume())
L=M.exterior_to_link()
L.view()
print(find_s3_slope(M))
print(L.exterior().volume())
print(L.alexander_polynomial(multivar=False,method='wirtinger'))

S=[1,1,2,1,2,2,2,1,2,1,1,1,2,1,2,2,1,2,1,1,1,2,2,1,2,2,2,1,2,2]



def word_to_braid(w):
	C, Id = spherogram.RationalTangle(1), spherogram.IdentityBraid(1)
	x = C | Id
	y = Id | C
	items = [x,y]
	gens = [items[x-1] for x in w]

	phi = reduce(operator.mul, gens[1:], gens[0])
	L = phi.denominator_closure()
	print(L.knot_floer_homology())
	E = L.exterior()
	print(E.volume())

word_to_braid(S)

eg=[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,3,2,3,2,3,2]
