import snappy
import csv
import sqlite3
import subprocess
import itertools

import os
import sys
sys.setrecursionlimit(100000)

import veering
import veering.taut
import sage
import sage.all
import ast

from queue import PriorityQueue
from itertools import count
import pathlib

from . import prepare

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


_DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "find_pA", "veering_census.sq3")
VeeringDB=snappy.database.ManifoldTable(table='census', mfld_hash=basic_hash, db_path=_DB_PATH)

def find(isosig):
	reader = csv.reader(open("veering_census_with_data.txt"), delimiter=" ")
	for row in reader:
		if isosig.split("_")[0] == row[0].split("_")[0]:
			yield row

#M is a manifold, and drilled_cusps is a subset of the unfilled cusps whech
#represent curves that have bene drilled.
#We want to eliminate some filled cusps from the underlying ideal tetrahedron.
#We need to return the new indices of the drilled cusps.
#We always need at least one cusp, so  we need a special handler if we are about
#to fill all the cusps.
def fill_unnecessary_cusps(M, drilled_cusps):
	is_complete = M.cusp_info('is_complete')
	filled_cusps = [i for i in range(M.num_cusps()) if not is_complete[i]]
	if len(filled_cusps) == M.num_cusps() and len(filled_cusps) > 0:
		filled_cusps = filled_cusps[0:len(filled_cusps)-1]

	# Each drilled cusp is complete (is_complete=True), so never in filled_cusps.
	# After removing filled_cusps, its new index decreases by the count of
	# filled cusps with smaller index.
	new_drilled_cusps = [d - sum(1 for c in filled_cusps if c < d) for d in drilled_cusps]
	return M.filled_triangulation(filled_cusps), new_drilled_cusps

def iterate_candidate_flows2(M, max_drill=2, max_segments=6, max_tets=20):
	pq = PriorityQueue()
	#store how many drillings we did and priority = num_tetrahedra

	#counter to disambiguate priorities
	unique = count()

	M.set_peripheral_curves('fillings') #We will maintain the invariant that all filled slopes are meridional.

	is_complete = M.cusp_info('is_complete')
	filled_cusps = [i for i in range(M.num_cusps()) if not is_complete[i]]

	for k in range(max_drill):
		for _drilled_cusps in itertools.combinations(filled_cusps, k):
			drilled_cusps = list(_drilled_cusps)
			N=M.copy()
			for cusp_ind in drilled_cusps:
				N.dehn_fill((0,0), cusp_ind)
			N,drilled_cusps =fill_unnecessary_cusps(N, drilled_cusps)
			pq.put(((N.num_tetrahedra(), next(unique)), (N, drilled_cusps)))

	while not pq.empty():
		priority, (M,drilled_cusps) = pq.get()

		try:
			M.volume()
			for x in VeeringDB.siblings(M):
				try:
					if M.is_isometric_to(x):
						yield M, x.name(), drilled_cusps 
				except RuntimeError as e:
					pass
					#print(e)
		except ValueError as e:#in case M is not hyperbolic
			print(e, file=sys.stderr)

		if len(drilled_cusps) < max_drill:
			for curve in M.dual_curves(max_segments=max_segments):
				N=M.copy().drill(curve)
				new_drilled_cusps = drilled_cusps + [N.num_cusps()-1]
				N, new_drilled_cusps=fill_unnecessary_cusps(N, new_drilled_cusps)
				#the newly drilled cusps has the last available index.
				if N.num_tetrahedra() < max_tets:
					pq.put(((N.num_tetrahedra(), next(unique)),(N, new_drilled_cusps)))


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
							yield N, x.name(), list(range(M.num_cusps(), N.num_cusps()))
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


def pA_flows(MM, count=10, max_drill=2, maxlength=3, max_segments=6, max_tets=20, method='geodesic', return_isom=False, return_prong_counts=False):
	D=dict()
	seen=set()

	if method=='geodesic':
		it = iterate_candidate_flows(MM, count=count, max_drill=max_drill, maxlength=maxlength)
	elif method=='combinatorial':
		it = iterate_candidate_flows2(MM, max_drill=max_drill, max_segments = max_segments, max_tets = max_tets)


	for N,isosig,drilled_cusps in it:
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

			#Transform the filling slopes on N to filling slopes on M
			filling_slopes = [(0,0) for i in range(N.num_cusps())]
			for i in drilled_cusps:
				filling_slopes[isom.cusp_images()[i]] = isom.cusp_maps()[i]*sage.all.vector([1,0])
			for i in range(len(filling_slopes)):
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

			if all(pc >= 2 for pc, s in zip(prong_counts,filling_slopes) if s != (0,0)):
				#print("degeneracy", degeneracy)
				if return_isom or return_prong_counts:
					extras = {}
					if return_prong_counts:
						extras["prong_counts"] = [int(pc) for pc in prong_counts]
					if return_isom:
						# The cusps of N that are not drilled correspond 1-to-1 with MM's cusps.
						# Use isom (N→M) to build the M→MM cusp mapping:
						# for each cusp j of M, find its preimage i in N; if i is an original
						# (non-drilled) cusp, it maps to cusp i of MM with map isom.cusp_maps()[i].inverse().
						drilled_set = set(drilled_cusps)
						NtoM_imgs = list(isom.cusp_images())
						MtoN = {NtoM_imgs[i]: i for i in range(N.num_cusps())}
						cusp_images_MtoMM = []
						cusp_maps_MtoMM   = []
						for j in range(M.num_cusps()):
							i = MtoN[j]
							if i not in drilled_set:
								cusp_images_MtoMM.append(i)
								cusp_maps_MtoMM.append(isom.cusp_maps()[i].inverse())
							else:
								cusp_images_MtoMM.append(None)
								cusp_maps_MtoMM.append(None)
						extras["isom"] = (cusp_images_MtoMM, cusp_maps_MtoMM)
					yield s, extras
				else:
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

