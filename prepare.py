import veering
import regina
import veering.veering_tri
import veering.taut as taut
import veering.transverse_taut
#import sage
import numpy as np
import ast


import sys
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering/scripts")
sys.set_int_max_str_digits(0)
import boundary_triangulation


from veering import file_io

#veering_knots_with_data = file_io.parse_data_file("knot_hom_census_with_data.txt")

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

def flatten(lists):
	return [x for l in lists for x in l]

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


def prepare_example(vt, isosig="test", longitude=None):
	info_file = open("batch/" + isosig + ".txt",'w')

	fname = "batch/" + isosig + ".pdf"
	fans = veering.transverse_taut.edge_side_face_collections(vt.tri,vt.angle)
	print("fans=" + str(fans), file=info_file)
	print("face_coorientations=OffsetArrays.Origin(0)(" + str(veering.transverse_taut.convert_tetrahedron_coorientations_to_faces(vt.tri,vt.coorientations)) + ")", file=info_file)

	args = {'style':'ladders', 'draw_boundary_triangulation':True, 'draw_triangles_near_poles': False, 'ct_depth':-1, 'ct_epsilon':0.03, 'global_drawing_scale': 4, 'delta': 0.2, 'ladder_width': 10.0, 'ladder_height': 30.0, 'draw_labels': True}
	
	bt=boundary_triangulation.generate_boundary_triangulation(vt.tri,vt.angle,args=args)
	bt.generate_canvases(args=args)#this does some important setup...
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
	
	first_rungs = [runglist[0][0] for runglist in bt.rungs()]
	print("firstrungs = " + str(first_rungs), file=info_file)

	print("rungs = " + str(bt.rungs()), file=info_file)

	print("alledges = " + str(bt.all_edges()), file=info_file)

	"""
	print("weights=[", file=info_file)
	for e in weights.keys():
		print(str(e) + "=>" + str(weights[e]) + ",", file=info_file)
	print("]", file=info_file)
	"""

	top_bot_embeddings = veering.transverse_taut.top_bottom_embeddings_of_faces(vt.tri, vt.angle, vt.coorientations)
	top_bot_pairs = [(x.simplex().index(),y.simplex().index()) for (x,y) in zip(*top_bot_embeddings)]

	print("top_bot_pairs = " + str(top_bot_pairs), file=info_file)

	#find_longitudes(vt)

	if longitude is not None:
		l=[0 for i in range(vt.tri.countTriangles())]
		for i,j in longitude: #i is a tetrahedron number, j is the face index
			l[vt.tri.tetrahedron(i).triangle(j).index()] += 1
		print("longitude=OffsetArrays.Origin(0)(" + str(l) + ")", file=info_file)
	else:
		print("longitude=nothing", file=info_file)
	info_file.close()

def prepare_by_isosig(isosig):
	x = taut.isosig_to_tri_angle(isosig)
	v= veering.veering_tri.veering_triangulation(*x)
	prepare_example(v, isosig = isosig)

if __name__ == "__main__":
	if False:
		n=3
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

	#prepare_by_isosig("eLMkbcddddedde_2100")
	#prepare_by_isosig("gvLQQcdeffeffffaafa_201102")
	#prepare_by_isosig("gLLAQcdecfffhsermws_122201")
	#prepare_by_isosig("fLLQcbecdeepuwsua_20102")	
	#prepare_by_isosig("fLLQcbeddeehhbghh_01110")
	#prepare_by_isosig("gLLPQbefefefhhhhhha_011102")
	#prepare_by_isosig("gLLPQcdfefefuoaaauo_022110")
	#prepare_by_isosig("jLvvQQQbhigghihgixaxxvvvvcc_102222010")
	#prepare_by_isosig("hLLzQkcdefggfghspadgsg_1222011")
	#prepare_by_isosig("iLLALQcbcdefghhhtsfxjoxqt_20102120")
	#prepare_by_isosig("gLLPQccdfeffhggaagb_201022")
	#prepare_by_isosig("iLLAMMccdecffghhhsermstqs_12220120")
	#prepare_by_isosig("jLLLzQQbeeeghihiixxaavvvvcv_211120000")
	#prepare_by_isosig("iLLLQPcbeegefhhhhhhahahha_01110221")
	for isosig in sys.argv[1:]:
		prepare_by_isosig(isosig)
