from flipper import *
from flipper.kernel import *
from collections import defaultdict
from unionfind import unionfind
from functools import reduce
import operator
import snappy
import snappy.snap.peripheral.peripheral as periph
import veering
import veering.veering_tri
import regina
import sys
import numpy as np
import random



sys.path.append("/home/jonathan/repo/Veering/scripts")
sys.set_int_max_str_digits(0)
import boundary_triangulation
import prepare

class DehnTwist():
	def __init__(self, i, ):
		"""A Dehn twist between columns i and i+1"""
		self.i=i
		self.edges=[]
	
	def push(self,e):
		self.edges.append(e)

def ncomponents(braidgens):
	"""Given a sequence of braid generators, label the different components"""
	n=max(braidgens)+2
	u=unionfind(n)
	r=list(range(n))
	for i in braidgens:
		r[i],r[i+1]=r[i+1],r[i]
	for i,j in enumerate(r):
		u.unite(i,j)
	ret = [0 for i in range(n)]
	for (i,g) in enumerate(u.groups()):
		for j in g:
			ret[j]=i
	return len(u.groups())

def make_triangulation(braidgens):
	n=max(braidgens)+2 #braid index
	v=flipper.kernel.Vertex(0)
	top_edges = [[] for i in range(n)]
	dehn_twists = [[] for i in range(n)]
	edge_count=[0]
	u=unionfind(10*len(braidgens))
	triangles=[]
	def new_edge():
		edge_count[0]+=1
		return edge_count[0]-1
	def new_triangle(e1,e2,e3):
		triangles.append((e1,e2,e3))
		return (e1,e2,e3)
	def glue_edges(e1,e2):
		u.unite(e1,e2)
	def representatives():
		d=dict()
		comps = list(filter(lambda x: len(x)!=1 or x[0] < edge_count[0], map(tuple,u.groups())))
		for i,x in enumerate(comps):
			for j in x:
				d[j]=i
		n_unique = len(comps)
		return n_unique,d

	def islastedge(i,count):
		return not ((i in braidgens[count+1:]) or (i-1 in braidgens[count+1:]))
	def isfirstbridge(i,count):
		return not (i in braidgens[:count])
	def islastbridge(i,count):
		return not (i in braidgens[count+1:])

	all_twists = []

	for count, i in enumerate(braidgens):
		bridgei=None
		bridgej=None
		lastdt = None
		nextdt = None

		leftdt = None
		rightdt = None

		if not isfirstbridge(i,count):
			"""then we've finished a Dehn twist"""
			lastdt = dehn_twists[i][-1]
		if islastbridge(i,count):
			dehn_twists[i].append(None)
		else:
			nextdt = DehnTwist(i)

		if i>0 and len(dehn_twists[i-1])!=0:
			leftdt = dehn_twists[i-1][-1]
		if i < n-1 and len(dehn_twists[i+1])!=0:
			rightdt = dehn_twists[i+1][-1]

		if len(top_edges[i])==0:
			bridgei = new_edge()
			top_edges[i].append(bridgei)
			if nextdt != None:
				nextdt.push(bridgei)
		elif islastedge(i,count):
			bridgei = top_edges[i][-1]
		else:
			bridgei = new_edge()
			ei=new_edge()
			ti=new_triangle(ei,top_edges[i][-1],bridgei)
			top_edges[i].append(ei)
			if lastdt != None:
				lastdt.push(bridgei)
			if nextdt != None:
				nextdt.push(bridgei)
				nextdt.push(ei)
			if leftdt != None:
				leftdt.push(ei)

		if len(top_edges[i+1])==0:
			bridgej = new_edge()
			top_edges[i+1].append(bridgej)
			if nextdt != None:
				nextdt.push(bridgej)
		elif islastedge(i+1,count):
			bridgej = top_edges[i+1][-1]
		else:
			bridgej = new_edge()
			ej=new_edge()
			ti=new_triangle(ej,top_edges[i+1][-1],bridgej)
			top_edges[i+1].append(ej)
			if lastdt != None:
				lastdt.push(bridgej)
			if nextdt != None:
				nextdt.push(bridgej)
				nextdt.push(ej)
			if rightdt != None:
				rightdt.push(ej)

		if nextdt != None:
			dehn_twists[i].append(nextdt)
		if lastdt != None:
			all_twists.append(lastdt)
		glue_edges(bridgei,bridgej)
	
	n_unique, labeller = representatives()

	edges = []
	for i in range(n_unique):
		#e=Edge(v,v,i)
		#edges.append([e,e.reversed_edge])
		edges.append([i, i+n_unique])
	
	def reverse_edge(k):
		if k< n_unique:
			return k+n_unique
		else:
			return k-n_unique
	
	def get_edge(i):
		return edges[labeller[i]].pop()

		
	_triangles = [(get_edge(i), get_edge(j), get_edge(k)) for i,j,k in triangles]
	for l in edges:
		assert len(l)==0
	u2=unionfind(2*n_unique)
	edge_vertex_labeller=[0 for i in range(2*n_unique)]
	for t in _triangles:
		for i in range(3):
			u2.unite(reverse_edge(t[i]), (t[(i+1) % 3]))
	
	for i,g in enumerate(u2.groups()):
		v=flipper.kernel.Vertex(i)
		for j in g:
			edge_vertex_labeller[j]=v
	
	assert len(u2.groups()) == ncomponents(braidgens)
	
	edges = []
	for i in range(n_unique):
		e=Edge(edge_vertex_labeller[i],edge_vertex_labeller[i+n_unique],i)
		edges.append([e,e.reversed_edge])

	triangles = [(get_edge(i), get_edge(j), get_edge(k)) for i,j,k in triangles]
	triangulation = Triangulation([Triangle(t) for t in triangles])
	
	encoded_twists=[]
	#for dt in all_twists:
	for column in dehn_twists:
		for dt in column:
			if dt==None:
				continue
			e=set([labeller[i] for i in dt.edges])
			weights=[1 if i in e else 0 for i in range(n_unique)]
			assert triangulation.lamination(weights).is_curve()
			encoded_twists.append(triangulation.lamination(weights).encode_twist())
	
	phi = reduce(operator.mul, encoded_twists[1:], encoded_twists[0])
	
	return triangulation, phi

def enumerate_braids(braid_index,maxlength):
	if maxlength==0:
		yield []
	else:
		for t in enumerate_braids(braid_index,maxlength-1):
			for j in range(1,braid_index):
				yield t + [j]


def random_braid(braid_index, length):
	return random.choices(range(1,braid_index), k=length)


def siddhi(m):
	return [2,1,3,2,3,3,2,1,3,2,3,3,3,2,1,3,2,3,3,2,3] + reduce(operator.add,[[2,3,2] for i in range(m)],[])

thebestknot = [1,2,2,1,1,2,2,2,2,2,2,2]
bojun=[1,1,2,1,2,2,2,1,2,1,1,1,2,1,2,2,1,2,1,1,1,2,2,1,2,2,2,1,2,2]

siddhi_mom = [1,2,3]*7 + [4,3,2,2,3,4]

small_examples=[
[1,1,1],	
[1,1,1,1,1],
[1,1,1,1,1,1,1],	
[1,1,1,2,1,1,1,2],
[1,1,1,1,1,1,1,1,1],	
[1,1,1,1,1,2,1,1,1,2],	 
[1,1,1,1,2,1,1,1,2,2],
[1,1,1,2,2,1,1,2,2,2],
[1,1,1,1,1,1,1,1,1,1,1],
[1,1,2,2,1,3,2,2,2,3,3],
[1,2,2,1,1,2,2,2,2,2,2,2],
[1,2,2,2,2,1,1,2,2,2,2,2],
[1,2,2,2,2,2,2,1,1,2,2,2],
[1,1,1,2,2,1,1,2,2,2,2,2],
[1,1,1,2,2,2,2,1,1,2,2,2],
[1,2,2,1,1,1,1,2,2,2,2,2],
[1,1,1,2,2,2,1,1,1,2,2,2]]

def sivek(n):
	return reduce(operator.add,[[4,3,2,1] for i in range(5*n-1)],[]) + [3,2,1,3,2,1]


def find_s3_slope(M,unfilled_index, n=20):
	current_fillings=[x.filling for x in M.cusp_info()]
	if M.num_cusps()==1:
		slope_candidates = M.short_slopes()[0]
	else:
		slope_candidates = [(p,q) for p in range(-n,n) for q in range(0,n)]

	for ss in slope_candidates:
		M.dehn_fill(ss,unfilled_index)
		G=M.fundamental_group()
		if len(G.generators())==0:
			return ss
	assert False
	return None

def find_s3_slope2(M):
	ret = []
	for i in range(M.num_cusps()):
		M.dehn_fill((0,0),i)
	for ss in itertools.product(*M.short_slopes(length=6)):
		M.dehn_fill(list(ss))
		G=M.fundamental_group()
		print(ss)
		if len(G.generators())==0:
			ret.append(ss)
	return ret


def analyze(bword, name="test"):
	print()
	print("braid word: ", bword)
	tau,phi = make_triangulation([x-1 for x in bword])
	print(str(len(tau.triangles)) + " triangles")
	#print(phi.nielsen_thurston_type())
	#print(phi.invariant_lamination())
	print(phi.stratum())
	bun=phi.bundle()
	s=bun.snappy_string()
	M=snappy.Manifold(bun)
	print(M)
	print(M.volume())
	print(M.identify())
	degen = bun.degeneracy_slopes()
	fibre = bun.fibre_slopes()

	current_fillings = [x.filling for x in M.cusp_info()]
	unfilled_index=0
	for i in range(len(current_fillings)):
		if current_fillings[i] == (0.0,0.0):
			unfilled_index=i
			break


	try:
		find_s3_slope(M,unfilled_index)
		meridian = [x.filling for x in M.cusp_info()]


		for i in range(M.num_cusps()):
			print("cusp "+str(i))
			if i==unfilled_index:
				A=np.linalg.inv(np.transpose(np.array([fibre[i],meridian[i]])))
			else:
				A=np.linalg.inv(np.transpose(np.array([fibre[i],degen[i]])))

			print("degen_slope " + str(np.matmul(A,np.array(degen[i]))))
			print("fiber_slope " + str(np.matmul(A,np.array(fibre[i]))))
			print("s3_slope " + str(np.matmul(A,np.array(meridian[i]))))
	except:
		print("failed to find s3 slope")

	M=snappy.Manifold(bun)
	tri=regina.SnapPeaTriangulation(s)
	angle=[1 for i in range(M.num_tetrahedra())] #flipper arranges that all the tetrahedra are flattened in the same way
	longitude =[(x[0].label,x[1](3)) for x in bun.immersion.values()]

	v=veering.veering_tri.veering_triangulation(tri,angle)
	prepare.prepare_example(v, isosig=name, longitude=longitude)

	info_file = open("batch/" + name + ".info.txt",'w')
	print(phi.stratum(), file=info_file)
	print(bun.snappy_string(), file=info_file)

	info_file.close()
	return bun

#for bword in small_examples:
#	try:
#		analyze([x-1 for x in bword])
#	except Exception as e:
#		print(e)

#bun=analyze(siddhi_mom,name="siddhi_mom")
#bun=analyze(siddhi(2),fname="siddhi_vt.pdf")

#bun = analyze(bojun,name="bojun")
bun=analyze(thebestknot,name="pretzel")


#n=2
#bun = analyze(siddhi(n),name="siddhi"+str(n))





#print(veering.transverse_taut.edge_side_face_collections(tri,angle))

#bword=[2-x for x in bojun]
#bword = [x-1 for x in thebestknot]
#for maxlength in range(15,20):
#	for i in range(10):
#		try:
#			analyze(random_braid(4,maxlength))
#			pass
#		except Exception as e:
#			print(e)



