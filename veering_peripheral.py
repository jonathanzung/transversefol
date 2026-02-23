import veering
import regina
import veering.veering_tri
import veering.taut as taut
import veering.transverse_taut
import itertools

import snappy
import snappy.snap.peripheral.peripheral as periph
import sys
sys.path.append("/home/jonathan/Dropbox/repo/Veering/scripts")
sys.set_int_max_str_digits(0)
import boundary_triangulation

"""
from snappy.snap import t3mlite as t3m
from snappy.snap.peripheral import link, dual_cellulation

def peripheral_curve_from_snappy(dual_cell, cusp_indices, snappy_data):
	ncusps = max(map(max,cusp_indices)) + 1
	D = dual_cell
	T = D.dual_triangulation
	M = T.parent_triangulation
	data = snappy_data
	weights = [len(D.edges)*[0] for i in range(ncusps)] #We are going to return the weights on the edges of the dual cell decomposition of \partial M
	print("ntetrahedra", len(M.Tetrahedra))
	for tet_index, tet in enumerate(M.Tetrahedra): #for each tetrahedron
		for x in t3m.EdgeFacePairs:
			pass
		for vert_index, V in enumerate(t3m.ZeroSubsimplices): #for each vertex in the tetrahedron
			triangle = tet.CuspCorners[V] #a triangle in the cusp triangulation
			current_cusp = cusp_indices[tet_index][vert_index]
			#print(cusp_indices[tet_index][vert_index])
			sides = triangle.oriented_sides()
			for tri_edge_index, tet_edge in enumerate(link.TruncatedSimplexCorners[V]): #for each edge in the corresponding truncation triangle
				#print((triangle.index,tri_edge_index))
				tet_face_index = t3m.ZeroSubsimplices.index(tet_edge ^ V)
				side = sides[tri_edge_index]
				global_edge = side.edge()
				if global_edge.orientation_with_respect_to(side) > 0:
					print()
					print(global_edge.__dir__())
					print(global_edge.face_index)
					print(global_edge.sides)
					print(global_edge.index)
					print(vert_index,tet_face_index)

					dual_edge = D.from_original[global_edge]
					#print(data[tet_index])
					weight = data[tet_index][4*vert_index + tet_face_index]
					weights[current_cusp][dual_edge.index] = -weight

	# Sanity check
	total_raw_weights = sum([sum(abs(x) for x in row) for row in data])
	assert 2*sum(abs(w) for x in weights for w in x) == total_raw_weights
	return weights
"""




#isosig="iLMzMPcbcdefghhhhhhhxxqdl_12211002"
isosig="gvLQQcdeffeffffaafa_201102"
get_peripheral_weights(isosig)
