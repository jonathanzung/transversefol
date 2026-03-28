import snappy
import csv
import sqlite3
import subprocess
import itertools

import os
import sys
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering/scripts")
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering")
sys.setrecursionlimit(100000)

import veering
import veering.taut
import sage
import sage.all
import ast

from queue import PriorityQueue
from itertools import count
import pathlib



sys.path.insert(0,"/home/jonathan/Dropbox/jonathan/transversefol")
import prepare

def basic_hash(manifold, digits=6):
	return to_str_at_prec(manifold.volume(), digits) + " " + repr(manifold.homology())

def to_str_at_prec(x, d):
	return ('%.' + repr(d) + 'f') % x

def generate_census():
	subprocess.run(["rm", "veering_census.sq3"]) 
	conn = sqlite3.connect("veering_census.sq3")

	cur = conn.cursor()
	reader = csv.reader(open("veering_census_with_data.txt"), delimiter=" ")
	data = []
	for i, row in enumerate(reader):
		isosig = row[0].split("_")[0]
		M=snappy.Manifold(isosig)


		cusps = M.cusp_info('is_complete').count(True)
		H = M.homology()
		betti = H.betti_number()
		torsion = [c for c in H.elementary_divisors() if c != 0]
		h=basic_hash(M)
		data.append((row[0], isosig, int(i), float(M.volume()), int(cusps), int(betti), str(torsion), h))

	cur.execute("CREATE TABLE census (name TEXT, triangulation TEXT, id INTEGER, volume REAL, cusps INT, betti INT, torsion TEXT, hash INT)")
	cur.executemany("INSERT INTO census VALUES(?,?,?,?,?,?,?,?)", data)
	conn.commit()

#generate_census()


_DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "veering_census.sq3")
VeeringDB=snappy.database.ManifoldTable(table='census', mfld_hash=basic_hash, db_path=_DB_PATH)

def find(isosig):
	reader = csv.reader(open("veering_census_with_data.txt"), delimiter=" ")
	for row in reader:
		if isosig.split("_")[0] == row[0].split("_")[0]:
			yield row

def fill_cusps(M):
	is_complete = M.cusp_info('is_complete')
	filled_cusps = [i for i in range(M.num_cusps()) if not is_complete[i]]
	if len(filled_cusps) == M.num_cusps() and len(filled_cusps) > 0:
		filled_cusps = filled_cusps[0:len(filled_cusps)-1]
	return M.filled_triangulation(filled_cusps)

def iterate_candidate_flows2(M, max_drill=2, max_segments=6, max_tets=20):
	pq = PriorityQueue()
	#store how many drillings we did and priority = num_tetrahedra

	#counter to disambiguate priorities
	unique = count()
	M.set_peripheral_curves('fillings')
	for k in range(max_drill):
		for to_drill in itertools.combinations(range(M.num_cusps()), k):
			N=M.copy()
			for cusp_ind in to_drill:
				N.dehn_fill((0,0), cusp_ind)
			N=fill_cusps(N)
			pq.put(((N.num_tetrahedra(), next(unique)), (N, k)))

	while not pq.empty():
		priority, (M,ndrill) = pq.get()

		try:
			M.volume()
			for x in VeeringDB.siblings(M):
				try:
					if M.is_isometric_to(x):
						yield M, x.name(), ndrill 
				except RuntimeError as e:
					pass
					#print(e)
		except ValueError as e:#in case M is not hyperbolic
			print(e, file=sys.stderr)

		if ndrill < max_drill:
			for curve in M.dual_curves(max_segments=max_segments):
				N=M.copy().drill(curve)
				N=fill_cusps(N)
				if N.num_tetrahedra() < max_tets:
					pq.put(((N.num_tetrahedra(), next(unique)),(N, ndrill+1)))


"""
yields N (drilled manifold), x (veering isosig), n (number of drillings)
"""


def iterate_candidate_flows(M, count=10, max_drill=2, maxlength=3):
	short_words = [(i.word,i.length.real()) for i in M.length_spectrum_alt(count=count, bits_prec=200)]

	#print("short words:", short_words)

	for i in range(max_drill+1):
		for words in itertools.combinations(short_words,i):
			#print(words)
			if sum(i[1] for i in words) < maxlength:
				try:
					N=M.drill_words([i[0] for i in words], bits_prec=200)
					for x in VeeringDB.siblings(N):
						if N.is_isometric_to(x):
							yield N, x.name(), len(words)
				except snappy.geometric_structure.geodesic.check_away_from_core_curve.ObjectCloseToCoreCurve as e:
					pass
					#print(e)
				except snappy.drilling.exceptions.DrillGeodesicError as e:
					pass
					#print(e)
				except RuntimeError as e:
					print(e, file=sys.stderr)


L=[
"s137(5,4)",
"s460(6,1)", 
"s593(6,1)", 
"s614(5,1)", 
"s753(6,-1)",
"s956(4,1)",
"v1333(5,-1)",
"v3045(4,-1)", 
"t06114(5,1)", 
"t08155(5,-1)", 
"o9_12518(6,-1)",
"o9_12544(4,-3)",
"o9_13679(6,-1)",
"o9_14675(1,-5)",
"o9_15066(5,-1)",
"o9_22743(7,1)", 
"o9_30634(6,1)", 
"o9_36699(7,1)"]


def pA_flows(MM, count=10, max_drill=2, maxlength=3, max_segments=6, max_tets=20, method='geodesic'):
	D=dict()
	seen=set()

	if method=='geodesic':
		it = iterate_candidate_flows(MM, count=count, max_drill=max_drill, maxlength=maxlength)
	elif method=='combinatorial':
		it = iterate_candidate_flows2(MM, max_drill=max_drill, max_segments = max_segments, max_tets = max_tets)


	for N,isosig,n in it:
		tri, angle = veering.taut.isosig_to_tri_angle(isosig)
		assert tri.isOriented() #important to do it this way, because veering fixes the orientation, and then we can pass to snappy without any problems
		M=snappy.Manifold(tri)
		try:
			isoms = N.is_isometric_to(M, return_isometries=True)
		except RuntimeError as e:
			print(e, file=sys.stderr)
			continue
		assert len(isoms) >= 1
		for isom in isoms:

			#compute the filling slopes
			filling_slopes = [(0,0) for i in range(n)]
			for i in range(n):
				filling_slopes[isom.cusp_images()[i]] = isom.cusp_maps()[i]*sage.all.vector([1,0])
			for i in range(n):
				if filling_slopes[i][0] < 0:
					filling_slopes[i] = (-filling_slopes[i][0], -filling_slopes[i][1])

			""" Some checks
			Mtmp = snappy.Manifold(tri)
			assert n == Mtmp.num_cusps()
			Mtmp.dehn_fill(filling_slopes)
			assert Mtmp.is_isometric_to(MM)
			"""

			s = isosig + "_" + str(filling_slopes).replace(" ", "")
			if s in seen:
				continue
			seen.add(s)

			data = prepare.get_prep(isosig)
			degeneracy = data["degeneracy"]

			prong_counts = [abs(prepare.intersection_number(a,b)) for a,b in zip(filling_slopes, degeneracy)]

			if all(x >= 2 for x in prong_counts):
				yield s



#M=snappy.OrientableClosedCensus[0]
#M=snappy.Manifold("K9_48(0,1)")
#M=snappy.Manifold(L[4])

"""
name = "K12n242(6,1)"
M=snappy.Manifold(name)


with open("batch/" + str(name) + "_pAflows.txt","w") as f:
	isosigs = list(dict.fromkeys(filter(honest_pA, pA_flows(M, count=6, ndrill=3, maxlength=3))))
	for isosig in isosigs:
		print(isosig, file=f)
"""
#print(pA_flows(snappy.Manifold("L6a4"), maxdepth=1))

#import pandas as pd
#import csv
#df = pd.read_csv('/home/jonathan/Downloads/conjecture_data/floer/final_data/QHSpheres.csv')

#print(df)

if __name__ == "__main__":
	"""
	print(sys.argv)
	start=ast.literal_eval(sys.argv[1])
	finish=ast.literal_eval(sys.argv[2])
	for M in snappy.OrientableClosedCensus[start:finish]:
		fname = "hodgson_weeks_pA/" + str(M) + "_pAflows.txt"
		if pathlib.Path(fname).is_file():
			with open(fname, "r") as f:
				if len(f.readlines()) > 0:
					continue
		try:
			with open(fname,"w") as f:
				for isosig in pA_flows(M, count=6, max_drill=3, max_segments=6, maxlength=3, max_tets=20, method='combinatorial'):
					print(isosig, file=f)
					break
		except Exception as e:
			print(e)
	"""

	"""	
	for name in L:
		print(name)
		fname = "dunfield_list/" + name + "_pAflows.txt"
		M=snappy.Manifold(name)
		with open(fname,"w") as f:
			for isosig in pA_flows(M, count=6, max_drill=4, max_segments=6, maxlength=3, max_tets=20, method='combinatorial'):
				print(isosig, file=f)
	"""
