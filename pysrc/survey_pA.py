import snappy
import ast
import sys

sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering/scripts")
sys.path.insert(0,"/home/jonathan/Dropbox/repo/Veering")

import veering
from veering import file_io
import regina
import veering.veering_tri
import veering.taut as taut
import veering.transverse_taut

tot=0
i=0
for M in snappy.OrientableClosedCensus:
	try:
		with open("hodgson_weeks_pA/" + str(M) + "_pAflows.txt","r") as f:
			tot+=1
			tmp = f.readlines()
			if len(tmp) > 0:
				i+=1
				raw_isosig, face_coors, fillings= tmp[0].strip().split("_")
				fillings = ast.literal_eval(fillings)
				isosig = raw_isosig + "_" + face_coors

				tri, angle = taut.isosig_to_tri_angle(isosig)
				assert tri.isOriented() #important so that snappy doesn't change these orientations
				N=snappy.Manifold(tri)
				N.dehn_fill(fillings)

				try:
					assert N.is_isometric_to(M)
				except:
					i-=1
			print(i,tot)
	except FileNotFoundError:
		pass

print(tot-i)
print(i/tot)
