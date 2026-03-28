import sys
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering/scripts")
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering")

import veering
from veering import file_io
import regina
import veering.veering_tri
import veering.taut as taut
import veering.transverse_taut
import numpy as np
import ast
import json
import snappy
from snappy.snap import t3mlite as t3m
from snappy.snap.peripheral import link, dual_cellulation
from collections import defaultdict
from collections import Counter
import itertools
import boundary_triangulation
import ast
import time
import os

sys.set_int_max_str_digits(0)


#veering_knots_with_data = file_io.parse_data_file("data/knot_hom_census_with_data.txt")

CACHE_PATH="/home/jonathan/Dropbox/jonathan/transversefol/cache"

def sign(perm):
	l=4
	ret = 1
	for i in range(4):
		for j in range(i+1, 4):
			if perm[j] > perm[i]:
				ret = -ret
	return ret

"""
degen = bun2.degeneracy_slopes()
fibre = bun2.fibre_slopes()
find_s3_slope(M)
meridian = [x.filling for x in M.cusp_info()]

print(M)

for i in range(M.num_cusps()):
	print("cusp "+str(i))
	if i==0:
		A=np.linalg.inv(np.transpose(np.array([fibre[i],meridian[i]])))
	else:
		A=np.linalg.inv(np.transpose(np.array([fibre[i],degen[i]])))

	print("degen_slope " + str(np.matmul(A,np.array(degen[i]))))
	print("fiber_slope " + str(np.matmul(A,np.array(fibre[i]))))
	print("s3_slope " + str(np.matmul(A,np.array(meridian[i]))))

"""

"""
s=bun.snappy_string()
M=snappy.Manifold(s)
tri=regina.SnapPeaTriangulation(s)
angle=[1 for i in range(M.num_tetrahedra())] #flipper arranges that all the tetrahedra are flattened in the same way
"""

def btpoles(bt):
	return [ttpoles(tt) for tt in bt.torus_triangulation_list]

def ttpoles(tt):
	return [list(set([a for lu in L.ladder_unit_list for a in pole_labels(lu)])) for L in tt.ladder_list]

def pole_labels(lu):
	if lu.is_on_left():
		assert len(lu.right_vertices)==1
		return [face_label(lu,v) for v in list(lu.right_vertices)]
	else:
		return []
		#assert len(lu.left_vertices)==1
		#return [face_label(lu,v) for v in list(lu.left_vertices)]

def btrungs(bt):
	return [ttrungs(tt) for tt in bt.torus_triangulation_list]

def ttrungs(tt):
    return [list(set([a for lu in L.ladder_unit_list for a in rung_labels(lu)])) for L in tt.ladder_list]

def rung_labels(lu):
	if lu.is_on_left():
		return [face_label(lu,v) for v in list(lu.left_vertices)]
	else:
		return [face_label(lu,v) for v in list(lu.right_vertices)]

def face_label(lu, face_num):
	tet = lu.vt.tri.tetrahedron( lu.tet_num )
	triangle = tet.triangle(face_num)
	vert = tet.vertex(face_num)
	triangle_num = triangle.index()

	#names of the vertices of this cusp triangle
	vertex_names = list(lu.verts_C.keys()) 
	#missing_vertex is the vertex at this cusp
	missing_vertex = [x for x in range(4) if not x in lu.verts_C.keys()]
	face_index=tet.faceMapping(2,face_num).inverse()[missing_vertex[0]]
	return (triangle_num, face_index)


def tt_all_edges(tt):
	return list(set(a for L in tt.ladder_list for lu in L.ladder_unit_list for a in all_labels(lu)))

def bt_all_edges(bt):
	return [tt_all_edges(tt) for tt in bt.torus_triangulation_list]


def all_labels(lu):
	return [face_label(lu,v) for v in list(lu.left_vertices) + list(lu.right_vertices)]


def get_peripheral_weights(isosig):
	tri, angle = taut.isosig_to_tri_angle(isosig)
	assert tri.isOriented() #important so that snappy doesn't change these orientations

	M=snappy.Triangulation(tri)
	cusp_indices, data = M._get_cusp_indices_and_peripheral_curve_data()

	#print(M._to_string())

	#cusp_indices tells you which cusp each vertex lives in

	ncusps = M.num_cusps()
	ntets = M.num_tetrahedra()
	assert len(data)==ntets*4
	assert len(cusp_indices)==ntets

	for row in data:
		for j in range(4):
			assert row[4*j+j]==0

	#see cython/core/triangulation.pyx
	#The rows congruent to 0 mod 4 are meridian
	# 1 mod 4 are left handed meridian
	# 2 mod 4 are longitude
	# 3 mod 4 are left handed longitude

	meridian_data = [data[i] for i in range(0, len(data), 4)]
	longitude_data = [data[i] for i in range(2, len(data), 4)]

	meridian_dict=defaultdict(lambda: 0)
	longitude_dict=defaultdict(lambda: 0)

	#print("ntets", len(tri.tetrahedra()))
	for tet in tri.tetrahedra():
		assert len(meridian_data[tet.index()])==16
		assert len(longitude_data[tet.index()])==16
		assert sum(meridian_data[tet.index()])==0
		assert sum(longitude_data[tet.index()])==0
		for v in range(4):
			#check that at each face of the cusp triangulation, the number
			#of times the meridian enters equals the number of times it exits
			assert sum([meridian_data[tet.index()][4*v+f] for f in range(4)])==0
			assert sum([longitude_data[tet.index()][4*v+f] for f in range(4)])==0
		for v in range(4):
			for f in range(4):
				if v != f:
					#we want to get the fth face of this tetrahedron
					find = tet.triangle(f).index()
					#print(tet.triangle(f).embeddings())
					#print(tet.triangleMapping(f))
					k = tet.triangleMapping(f).inverse()[v]
					assert tet.triangleMapping(f)[3]==f
					assert tet.triangleMapping(f)[k]==v
					assert k in [0,1,2]
					#check which way the orientation is going
					orient = sign(tet.triangleMapping(f))

					if orient == 1:
						meridian_dict[(find,k)]=meridian_data[tet.index()][4*v + f]
						longitude_dict[(find,k)]=longitude_data[tet.index()][4*v + f]
				else:
					assert meridian_data[tet.index()][4*v+f]==0
					assert longitude_data[tet.index()][4*v+f]==0
	
	for tet in tri.tetrahedra():
		for v in range(4):
			tmp_merid = 0
			tmp_long = 0
			for f in range(4):
				if v != f:
					k = tet.triangleMapping(f).inverse()[v]
					find = tet.triangle(f).index()
					tmp_merid += sign(tet.triangleMapping(f)) * meridian_dict[(find,k)]
					tmp_long += sign(tet.triangleMapping(f)) * longitude_dict[(find,k)]
			assert tmp_merid == 0
			assert tmp_long == 0
	
		
	return ([(key,val) for key,val in meridian_dict.items()], [(key,val) for key,val in longitude_dict.items()])

def intersection_number(x,y):
	return x[0]*y[1] - x[1] * y[0]

def flatten(lists):
	return [x for l in lists for x in l]

"""
def find_longitudes(vt):
	fans = veering.transverse_taut.edge_side_face_collections(vt.tri,vt.angle)

	relations = [([x[0] for x in f1], [x[0] for x in f2]) for (f1,f2) in fans]


	l = sorted(flatten(flatten(relations)))
	n = max(l)+1 #number of tetrahedra
	assert len(l) == (max(l)+1)*3
	#println(n)
	p = sage.all.MixedIntegerLinearProgram()
	w = p.new_variable(integer=True, nonnegative=True)


	M=np.zeros((len(relations),n))
	for (i,(l1,l2)) in enumerate(relations):
		for j in l1:
			M[i,j] += 1
		for j in l2:
			M[i,j] -= 1
	
	for i in range(len(relations)):
		p.add_constraint(sum([M[i,j]*w[j]]) == 0)
	
	for i in range(n):
		p.add_constraint(w[i]<= 100)
	
	print(p.polyhedron().integral_points())
"""


def compute_prep(isosig, longitude=None, draw_bt=False):
	x = taut.isosig_to_tri_angle(isosig)
	vt = veering.veering_tri.veering_triangulation(*x)


	fans = veering.transverse_taut.edge_side_face_collections(vt.tri,vt.angle)
	face_coorientations=veering.transverse_taut.convert_tetrahedron_coorientations_to_faces(vt.tri,vt.coorientations)
	#print("fans=" + str(fans), file=info_file)
	#print("face_coorientations=OffsetArrays.Origin(0)(" + str(face_coorientations) + ")", file=info_file)

	args = {'style':'ladders', 'draw_boundary_triangulation':True, 'draw_triangles_near_poles': False, 'ct_depth':-1, 'ct_epsilon':0.03, 'global_drawing_scale': 4, 'delta': 0.2, 'ladder_width': 10.0, 'ladder_height': 30.0, 'draw_labels': True}
	
	bt=boundary_triangulation.generate_boundary_triangulation(vt.tri,vt.angle,args=args)
	for tt in bt.torus_triangulation_list:
		for L in tt.ladder_list:
			L.ladder_origin = complex(0,0)
			L.calc_verts_C(args=args)
	
	if draw_bt:
		fname = "batch/" + isosig + ".pdf"
		bt.generate_canvases(args=args)#this does some important setup... I think the for loop above does the same
		bt.draw(fname,args=args)

	"""
	for tt in bt.torus_triangulation_list:
		print(tt.find_sideways(tt.tet_faces[0]))
		print(tt.ladder_holonomy)
		print(tt.sideways_holonomy)
		print(tt.sideways_once_holonomyi)
		print(tt.sideways_index_shift)
	"""
	
	"""
	weights = {}
	for e,v in bt.vertical_weights().items():
		weights[e] = [0,v]
	#print(bt.rungs())
	
	for runglist in bt.rungs():#iterate over the different boundary components
		for rung in runglist:
			weights[rung][0]=1
	"""
	
	rungs = btrungs(bt)
	alledges = bt_all_edges(bt)
	poles = btpoles(bt)


	"""
	print("firstrungs = " + str(first_rungs), file=info_file)

	print("rungs = " + str(btrungs(bt)), file=info_file)

	print("alledges = " + str(bt_all_edges(bt)), file=info_file)
	"""

	"""
	print("weights=[", file=info_file)
	for e in weights.keys():
		print(str(e) + "=>" + str(weights[e]) + ",", file=info_file)
	print("]", file=info_file)
	"""

	top_bot_embeddings = veering.transverse_taut.top_bottom_embeddings_of_faces(vt.tri, vt.angle, vt.coorientations)
	top_bot_pairs = [(x.simplex().index(),y.simplex().index()) for (x,y) in zip(*top_bot_embeddings)]
	meridian_dict, longitude_dict= get_peripheral_weights(isosig)
	
	mdict = dict(meridian_dict)
	ldict = dict(longitude_dict)
	

	"""
	for cusp in poles:
		for pole in cusp:
			#print(pole)
			print(sum([face_coorientations[e[0]]*mdict[e] for e in pole]),sum([face_coorientations[e[0]]*ldict[e] for e in pole]))
		print()
	print("---")
	"""

	degeneracy = [
			(-sum(face_coorientations[e[0]]*ldict[e] for pole in cusp[::2] for e in pole), #intersection with the longitude tells you how many meridians
			sum(face_coorientations[e[0]]*mdict[e] for pole in cusp[::2] for e in pole))
	for cusp in poles
	]
	#print("top_bot_pairs = " + str(top_bot_pairs), file=info_file)

	#find_longitudes(vt)

	"""
	if longitude is not None:
		l=[0 for i in range(vt.tri.countTriangles())]
		for i,j in longitude: #i is a tetrahedron number, j is the face index
			l[vt.tri.tetrahedron(i).triangle(j).index()] += 1
		print("longitude=OffsetArrays.Origin(0)(" + str(l) + ")", file=info_file)
	else:
		print("longitude=nothing", file=info_file)
	info_file.close()
	"""


	# For each tetrahedron: list of (triangle_index, orientation_sign) for its 4 faces.
	# orientation_sign is ±1 in the Regina/SnaPPy convention (NOT the veering convention).
	# To get the boundary operator in the veering co-orientation convention, multiply
	# each sign by face_coorientations[triangle_index].
	tet_faces = [
		[(tet.triangle(f).index(), sign(tet.triangleMapping(f))) for f in range(4)]
		for tet in vt.tri.tetrahedra()
	]

	data = {"fans": fans,
		"face_coorientations": face_coorientations,
		"poles": poles,
		"rungs": rungs,
		"alledges": alledges,
		"top_bot_pairs": top_bot_pairs,
		#"preferred_longitude": "nothing" if longitude == None else str(longitude),
		"meridian_dict": meridian_dict,
		"longitude_dict": longitude_dict,
		"degeneracy": degeneracy,
		"tet_faces": tet_faces,
		}

	"""
	if len(isosig.split("_")) >= 3:
		filling_slopes = ast.literal_eval(isosig.split("_")[2])
		#print(filling_slopes)
		assert len(filling_slopes) == len(degeneracy)
		d["prong_counts"] = [abs(intersection_number(a,b)) for a,b in zip(filling_slopes, degeneracy)]
	"""

	return data
	

def get_prep(isosig, draw_bt=False):
	cache_filename = os.path.join(CACHE_PATH, isosig + ".json")
	try:
		with open(cache_filename) as cachefile:
			#print("Retrieving from cache")
			data = {k: ast.literal_eval(v) for k,v in json.load(cachefile).items()}

	except FileNotFoundError:
		#print("Computing")
		data = compute_prep(isosig,draw_bt=draw_bt)
		with open(cache_filename,'w') as info_file:
			json.dump(dict((k,str(v)) for k,v in data.items()),
				info_file, indent=4)
	return data

if __name__ == "__main__":
	if False:
		n=2
		N=500
		veering_census_with_data = file_io.parse_data_file("veering_census_with_data.txt")
		f=open(f"batch/{n}cusp_manifest.txt","w")
		f2=open(f"batch/{n}cusp_manifest_data.txt","w")
		two_cusp_isosigs = []
		print("isosigs = [", file=f)
		i=0
		for _line in veering_census_with_data:
			line = _line.split()

			cusps = ast.literal_eval(line[7])
			if len(cusps)==n and line[1][0]=='F':
				print(line)
				two_cusp_isosigs.append(line[0])
				#x = taut.isosig_to_tri_angle(line[0])
				#v = veering.veering_tri.veering_triangulation(*x)
				#prepare_example(v, isosig = line[0])
				print("\"" + str(line[0])+"\",",file=f)
				print(line, file=f2)
				i += 1
				if i >= N:
					break
		print("]",file=f)
		f.close()	
		f2.close()

	if len(sys.argv) > 1:
		isosig = sys.argv[1]
		data = get_prep(isosig)
		print(json.dumps({k: str(v) for k, v in data.items()}))
	else:
		get_prep("gvLQQcdeffeffffaafa_201102")
