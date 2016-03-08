
""" 
This piece of code use two corresponding layers of material points coordinates and generate hex 8 nodes mesh.
"""


import xml.etree.cElementTree as ET
import numpy as np
#from gias.common import fieldvi
import subprocess
#from fieldwork.field import geometric_field
import copy

def multipleTwoArrays(A, B):
    Q = []
    for i in xrange(len(A)):
        q = A[i]*B[i]
        Q.append(q)
    return np.array(Q)
    

febFile = 'C:\Users\Mousa Kazemi\Phd_project\cartilage_growing\Utah_01.feb'

# subject 02 file names 
# tibia_mat_points_file = '/hpc/mkaz010/THESIS_PAPERS/PAPER_1/22_Embedding_cart_on_bone/tibia_subchondral_nodes_material_points_corrected.asc'
# tibia_mat_coords = np.loadtxt('/hpc/mkaz010/THESIS_PAPERS/PAPER_1/3_Calc_Cart_Thickness/02_HR_material_coord_tibia.asc', dtype = float)
# print len(tibia_mat_coords)
# tibia_mat_normals = np.loadtxt('/hpc/mkaz010/THESIS_PAPERS/PAPER_1/3_Calc_Cart_Thickness/02_HR_material_normals_tibia.asc', dtype = float)
# tibia_cart_thickness = np.loadtxt('/hpc/mkaz010/THESIS_PAPERS/PAPER_1/3_Calc_Cart_Thickness/02_HR_cart_thickness.asc', dtype = float)
# tibCartInCoords = tibia_mat_coords
# tibCartOutCoords = multipleTwoArrays(tibia_cart_thickness,tibia_mat_normals)+tibia_mat_coords
# print len(tibCartOutCoords)
# np.savetxt('/people/mkaz010/MAP_Project/MAP/Data/Paper_1_WFs/Hex8FromNodals/outputs/02_HR_tibCartOutCoords.asc', tibCartOutCoords )
#cart_base_nodal_co = np.loadtxt ('/people/mkaz010/MAP_Project/MAP/Data/Paper_1_WFs/PCA_2/Input/cartilage_thickness/femur_articular_mat_points_01.asc', dtype = float)
#mat_normals = np.loadtxt('/people/mkaz010/MAP_Project/MAP/Data/Paper_1_WFs/PCA_2/Input/cartilage_thickness/femur_material_normals_01.asc', dtype = float)
#cart_thickness = np.loadtxt('/people/mkaz010/MAP_Project/MAP/Data/Paper_1_WFs/PCA_2/Output/femur_cart_thickness_mean.txt', dtype = float)

#articular_surface_nodal_co = []
##for i in xrange(0, len(mat_normals)):
	##for j in xrange(0, len(cart_thickness)):
		##new_node = mat_normals[i]*cart_thickness[j]
		##cart_articulating_surface_node_coor.append(new_node)
		
#for i in xrange(0, len(mat_normals)):
	#new_node = mat_normals[i]*cart_thickness[i]+cart_base_nodal_co[i]
	#articular_surface_nodal_co.append(new_node)
		
SubBoneNodalNum = []
SubBoneNodalCoords = []
ArticularNodalNum = []
tree = ET.parse(febFile)
root = tree.getroot()
for child in root[2]:
	print child.tag, child.attrib

subSet = [] # element and nodal numbers of subcondral layer of cartilage mesh
artSet = [] # element nodal numbers of articular layer of cartilage mesh
subSetTib = []
artSetTib = []
artSet0 = []
subElemSet = root[2]
for child in subElemSet[1]:
	if subElemSet[1].attrib ['elset'] == 'Part18':
		t = child.attrib ['id']+ ','+child.text
		t1 = [float(word) for word in t.split(',')]
		subSet.append(t1)
#subSet = np.array(subSet)

# Subcondral nodal numbers
subSetNodNums = []
for child in subElemSet[10]:
	if subElemSet[10].attrib ['name'] == 'Nodeset01_fem_sub_cart':
		t0 = child.attrib ['id']
		subSetNodNums.append(t0)
subSetNodNums = [int(x) for x in subSetNodNums] 
#subSetNodNums = np.array(subSetNodNums)

for child in subElemSet[3]:
	if subElemSet[3].attrib ['elset'] == 'Part20':
		t = child.attrib ['id']+ ','+child.text
		t2 = [int(word) for word in t.split(',')]
		artSet.append(t2)
#artSet = np.array(artSet)

# Articular nodal numbers
artSetNodNums = []
for child in subElemSet[11]:
	if subElemSet[11].attrib ['name'] == 'Nodeset02_fem_art_cart':
		t0 = child.attrib ['id']
		#print t0
		artSetNodNums.append(t0)
artSetNodNums = [int(x1) for x1 in artSetNodNums] 
#artSetNodNums = np.array(artSetNodNums)

cartilageFiberDirect = []
for ind1 in xrange (len(subSetNodNums)):
	b = [subSetNodNums[ind1], artSetNodNums[ind1]]
	cartilageFiberDirect.append (b)	
#print np.shape(cartilageFiberDirect)

#cartilageFiberNodOrder = []
#for ind2 in xrange (len(cartilageFiberDirect)):
	#for ind3 in xrange (len(subSet)):
		#if 	(cartilageFiberDirect[ind2][0] in subSet[ind3][1:]) and (cartilageFiberDirect[ind2][1] in artSet[ind3][1:]):
			#print [(subSet[ind3], artSet[ind3]), (subSet[ind3].index (cartilageFiberDirect[ind2][0]), artSet[ind3].index(cartilageFiberDirect[ind2][1]))]

for child in subElemSet[4]:			# Tibial subchondral element-set
	if subElemSet[4].attrib ['elset'] == 'Part21':  
		t = child.attrib ['id']+ ','+child.text
		t1 = [float(word) for word in t.split(',')]
		subSetTib.append(t1)
subSetTib = np.array(subSetTib)

for child in subElemSet[6]:			#Tibial articular element-set
	if subElemSet[6].attrib ['elset'] == 'Part23': 
		t = child.attrib ['id']+ ','+child.text
		t2 = [int(word) for word in t.split(',')]
		artSetTib.append(t2)
artSetTib = np.array(artSetTib)

"""
In order to calculate the cartilage thickness and then changing that after fitting, first, we need to extract the direction of the cartilage 
fibers (from subcondral bone layer to articular surface. Then, we need to extract the coordinates of the nodes. By having the new cartilage thickness, 
we could re-calculate the thickness and replace the coordinates.  

# identify the direction of the nodals at cartilage thickness
# Note: node number are in element numbers and includes the repeatation of nodal numbers
#[subElemNum, artElemNum, subElemNum(n1),subElemNum(n2),artElemNum(n3),artElemNum(n4)]
"""
# 1) extract nodal numbers and their orders from subcondral layer to articular one
# femur
cartNodalDirect = [] 
for i in xrange (len(subSet)):
	d0 = [subSet[i][0], artSet[i][0], subSet[i][1],subSet[i][4], artSet[i][1], artSet[i][4]] 
	cartNodalDirect.append(d0)
	d1 = [subSet[i][0], artSet[i][0], subSet[i][2],subSet[i][3], artSet[i][2], artSet[i][3]] 
	cartNodalDirect.append(d1)
	d2 = [subSet[i][0], artSet[i][0], subSet[i][5],subSet[i][8], artSet[i][5], artSet[i][8]] 
	cartNodalDirect.append(d2)
	d3 = [subSet[i][0], artSet[i][0], subSet[i][6],subSet[i][7], artSet[i][6], artSet[i][7]] 
	cartNodalDirect.append(d3)
cartNodalDirect = np.array(cartNodalDirect, dtype = int)
subSetNodNumSorted = np.sort(np.unique (cartNodalDirect[:,2]))  # Subcondral cartilage layer sorted nodal numbers
artSetNodNumSorted = np.sort(np.unique (cartNodalDirect[:,-1])) # articular surface sorted nodal numbers 

print subSetNodNumSorted.shape
print artSetNodNumSorted.shape

# tibia
cartNodalDirectTib = [] 
for i in xrange (len(subSetTib)):
	d0 = [subSetTib[i][0], artSetTib[i][0], subSetTib[i][4],subSetTib[i][1], artSetTib[i][4], artSetTib[i][1]] 
	cartNodalDirectTib.append(d0)
	d1 = [subSetTib[i][0], artSetTib[i][0], subSetTib[i][3],subSetTib[i][2], artSetTib[i][3], artSetTib[i][2]] 
	cartNodalDirectTib.append(d1)
	d2 = [subSetTib[i][0], artSetTib[i][0], subSetTib[i][7],subSetTib[i][6], artSetTib[i][7], artSetTib[i][6]] 
	cartNodalDirectTib.append(d2)
	d3 = [subSetTib[i][0], artSetTib[i][0], subSetTib[i][8],subSetTib[i][5], artSetTib[i][8], artSetTib[i][5]] 
	cartNodalDirectTib.append(d3)
cartNodalDirectTib = np.array(cartNodalDirectTib, dtype = int)
subSetNodNumSortedTib = np.sort(np.unique (cartNodalDirectTib[:,2]))  # Subcondral cartilage layer sorted nodal numbers
artSetNodNumSortedTib = np.sort(np.unique (cartNodalDirectTib[:,-1])) # articular surface sorted nodal numbers 

print subSetNodNumSortedTib.shape
print artSetNodNumSortedTib.shape


# 2) extract the coordinates and calculate the thickness of cartilage at template mesh


nodCoords = []
for child in subElemSet[0]:
	#check
	if subElemSet[0].tag == 'Nodes':
		t = child.attrib ['id']+ ','+child.text
		t1 = [float(word) for word in t.split(',')]
		nodCoords.append(t1)		
nodCoords = np.array(nodCoords)

# femur
subCartNodCoords = []
artCartNodCoords = []
for n1 in xrange (len(nodCoords)):
	for n2 in xrange(len(subSetNodNumSorted)):
		if nodCoords[n1][0] ==subSetNodNumSorted[n2]:
			subCartNodCoords.append(nodCoords[n1])
	for n3 in xrange(len(artSetNodNumSorted)):
		if nodCoords[n1][0] ==artSetNodNumSorted[n3]:
			artCartNodCoords.append(nodCoords[n1])
subCartNodCoords = np.array(subCartNodCoords)
artCartNodCoords = np.array(artCartNodCoords)
print subCartNodCoords.shape
print artCartNodCoords.shape

# tibia
subCartNodCoordsTib = []
artCartNodCoordsTib = []
for n11 in xrange (len(nodCoords)):
	for n22 in xrange(len(subSetNodNumSortedTib)):
		if nodCoords[n11][0] == subSetNodNumSortedTib[n22]:
			subCartNodCoordsTib.append(nodCoords[n11])
	for n33 in xrange(len(artSetNodNumSortedTib)):
		if nodCoords[n11][0] ==artSetNodNumSortedTib[n33]:
			artCartNodCoordsTib.append(nodCoords[n11])
subCartNodCoordsTib = np.array(subCartNodCoordsTib)
artCartNodCoordsTib = np.array(artCartNodCoordsTib)
print subCartNodCoordsTib.shape
print artCartNodCoordsTib.shape
#3) calculate cart thickness wrt element numbers

def cartilageThickness( subNodCoords, artCartCoords):
	nodalUnitVector = []
	cartThickness = []
	#check
	if len (subNodCoords) != len(artCartCoords):
		print "different in Number of articular nodals & subcondral nodes"
	else:
		for v1 in xrange(len(subNodCoords)):
			diff = artCartCoords[v1][1:4] - subNodCoords [v1][1:4]
			diffMag = np.sqrt(diff[0]**2+diff[1]**2+diff[2]**2)
			thickness = [subNodCoords [v1][0], artCartCoords[v1][0], diffMag]
			cartThickness.append(thickness)
			unitV = diff/diffMag # unit vector of cartilage thickness
			unitVector = [subNodCoords [v1][0], artCartCoords[v1][0], unitV[0], unitV[1],unitV[2]]
			nodalUnitVector.append(unitVector)
	cartThickness = np.array(cartThickness) # including subcondral & articular nodal numbers and its corresponding cart thickness
	nodalUnitVector = np.array(nodalUnitVector) # including subcondral 
	return cartThickness, nodalUnitVector	
# femur	
T_Fem, U_Fem = cartilageThickness (subCartNodCoords, artCartNodCoords)
#tibia 
T_Tib, U_Tib = cartilageThickness (subCartNodCoordsTib, artCartNodCoordsTib)

		 
# 4) reconstruct 3 layers cartilage mesh and their nodal coordinates using another thickness field 

def cartReconstructFemur(subCartOldNodCoords, cartThicknesses, unitVectors):
	ReconstCartNodalDirect=[]
	nodalsCoords = []
	NodalCoordsOrdered = []
	nodalsCoordsReshape = []
	NodalCoordsOrdered = []
	
	"""
	# columns of subCartOldNodCoords (nod number, x, y, z)
	# columns of unitVectors (subNod number, artNod number, x,y,z) 
	# column of cartThicknesses (subNod number, artNod number, thickness magnetude)
	"""
		
	#check consistancy	
	if len(subCartOldNodCoords) !=len(unitVectors) and len(subCartOldNodCoords)!=len(cartThicknesses):
		
		
		print "inconsistancy in dataset lengthes"
	else:				
		for f in xrange (len(subCartOldNodCoords)):
			CartLayers =3
			thicknessCoef = np.dot(cartThicknesses[f][2], unitVectors[f][2:5])/CartLayers
			correctCoef = 1.05 # 5 percent increase of thickness		
			n0 = subCartOldNodCoords[f][1:4]
			n1 = n0+ thicknessCoef*correctCoef
			n2 = n1+ thicknessCoef*correctCoef
			n3 = n2+ thicknessCoef*correctCoef
			new = [unitVectors[f][0], unitVectors[f][1], n0[0], n0[1], n0[2] , n1[0], n1[1], n1[2], n2[0], n2[1], n2[2] , n3[0], n3[1], n3[2]] # [node number of(n0) , node number of (n3),n0,n1,n2,n3]
			ReconstCartNodalDirect.append(new)
	ReconstCartNodalDirect = np.array(ReconstCartNodalDirect)
	
	# adding four nodal numbers of the cartilage thickness along its fiber direction (from sub to art)
	
	for d1 in xrange (len(ReconstCartNodalDirect)):
		for d2 in xrange (len(cartNodalDirect)):
			if ReconstCartNodalDirect[d1][0] == cartNodalDirect[d2][2] and ReconstCartNodalDirect[d1][1] == cartNodalDirect[d2][5]:
				tem1 = np.hstack((cartNodalDirect[d2][2:], ReconstCartNodalDirect[d1][2:]))
				nodalsCoords.append(tem1)
	nodalsCoords = np.array(nodalsCoords)

	for r in xrange(len(nodalsCoords)):
		s1 = [nodalsCoords[r][0], nodalsCoords[r][4], nodalsCoords[r][5], nodalsCoords[r][6]]
		nodalsCoordsReshape.append(s1)
		s2 = [nodalsCoords[r][1], nodalsCoords[r][7], nodalsCoords[r][8], nodalsCoords[r][9]]
		nodalsCoordsReshape.append(s2)
		s3 = [nodalsCoords[r][2], nodalsCoords[r][10], nodalsCoords[r][11], nodalsCoords[r][12]]
		nodalsCoordsReshape.append(s3)
		s4 = [nodalsCoords[r][3], nodalsCoords[r][13], nodalsCoords[r][14], nodalsCoords[r][15]]	
		nodalsCoordsReshape.append(s4)
	nodalsCoordsReshape = np.array(nodalsCoordsReshape)
	print nodalsCoordsReshape.shape

	## get rid of the duplicated rows 
	#g = copy.deepcopy(nodalsCoordsReshape)
	#g = g[np.argsort(g[:,0]), :]

	##NodalCoordsOrdered = []
	##g = array([[5,3],[1,2],[1,7],[2,9],[7,1],[1,3], [2,2],[3, 4]])
	##g = g[np.argsort(g[:,0]), :]
	##print 'sorted', g
	#NodalCoordsOrdered.append(g[0])
	#for i in xrange (1, len(g)):
		#if np.all(g[i-1][0] != g[i][0]):
			#NodalCoordsOrdered.append(g[i])
	##print NodalCoordsOrdered
	##NodalCoordsOrdered = np.array(NodalCoordsOrdered)
	
	# get rid of the duplicated rows 
	g = copy.deepcopy(nodalsCoordsReshape)
	g = g[np.argsort(g[:,0]), :]
	NodalCoordsOrdered = []
	NodalCoordsOrdered.append(g[0])
	for i in xrange (len(g)-1):
		#for j in xrange (i+1,len(g)):
			if g[i][0] != g[i+1][0]:
				#	NodalCoordsOrdered.append(g[i])
					NodalCoordsOrdered.append(g[i+1])
					
	return np.array(NodalCoordsOrdered)
	
D_Fem = cartReconstructFemur (subCartNodCoords, T_Fem,U_Fem)
print D_Fem.shape


def cartReconstructTibia(subCartOldNodCoords, cartThicknesses, unitVectors):
	ReconstCartNodalDirect=[]
	nodalsCoords = []
	NodalCoordsOrdered = []
	nodalsCoordsReshape = []
	NodalCoordsOrdered = []
		
	#check consistancy	
	if len(subCartOldNodCoords) !=len(unitVectors) and len(subCartOldNodCoords)!=len(cartThicknesses):	
		print "inconsistancy in dataset lengthes"
	else:				
		for f in xrange (len(subCartOldNodCoords)):
			CartLayers =3
			thicknessCoef = np.dot(cartThicknesses[f][2], unitVectors[f][2:5])/CartLayers
			correctCoef = 1.05 # 5 percent increase of thickness		
			n0 = subCartOldNodCoords[f][1:4]
			n1 = n0+ thicknessCoef*correctCoef
			n2 = n1+ thicknessCoef*correctCoef
			n3 = n2+ thicknessCoef*correctCoef
			new = [unitVectors[f][0], unitVectors[f][1], n0[0], n0[1], n0[2] , n1[0], n1[1], n1[2], n2[0], n2[1], n2[2] , n3[0], n3[1], n3[2]] # [node number of(n0) , node number of (n3),n0,n1,n2,n3]
			ReconstCartNodalDirect.append(new)
	ReconstCartNodalDirect = np.array(ReconstCartNodalDirect)
	print ReconstCartNodalDirect.shape
	print ReconstCartNodalDirect[0]
	
	# adding four nodal numbers of the cartilage thickness along its fiber direction (from sub to art)
	
	for d1 in xrange (len(ReconstCartNodalDirect)):
		for d2 in xrange (len(cartNodalDirectTib)):
			if ReconstCartNodalDirect[d1][0] == cartNodalDirectTib[d2][2] and ReconstCartNodalDirect[d1][1] == cartNodalDirectTib[d2][5]:
				tem1 = np.hstack((cartNodalDirectTib[d2][2:], ReconstCartNodalDirect[d1][2:]))
				nodalsCoords.append(tem1)
	nodalsCoords = np.array(nodalsCoords)

	for r in xrange(len(nodalsCoords)):
		s1 = [nodalsCoords[r][0], nodalsCoords[r][4], nodalsCoords[r][5], nodalsCoords[r][6]]
		nodalsCoordsReshape.append(s1)
		s2 = [nodalsCoords[r][1], nodalsCoords[r][7], nodalsCoords[r][8], nodalsCoords[r][9]]
		nodalsCoordsReshape.append(s2)
		s3 = [nodalsCoords[r][2], nodalsCoords[r][10], nodalsCoords[r][11], nodalsCoords[r][12]]
		nodalsCoordsReshape.append(s3)
		s4 = [nodalsCoords[r][3], nodalsCoords[r][13], nodalsCoords[r][14], nodalsCoords[r][15]]	
		nodalsCoordsReshape.append(s4)
	nodalsCoordsReshape = np.array(nodalsCoordsReshape)
					
	# get rid of the duplicated rows 
	g = copy.deepcopy(nodalsCoordsReshape)
	g = g[np.argsort(g[:,0]), :]
	NodalCoordsOrdered = []
	NodalCoordsOrdered.append(g[0])
	for i in xrange (len(g)-1):
		#for j in xrange (i+1,len(g)):
			if g[i][0] != g[i+1][0]:
				#	NodalCoordsOrdered.append(g[i])
					NodalCoordsOrdered.append(g[i+1])
					
	return np.array(NodalCoordsOrdered)

D_Tib = cartReconstructTibia (subCartNodCoordsTib, T_Tib,U_Tib)
print D_Tib.shape


# 6) update the Feb file
#i = 0
#for child in root.findall(".//Nodes/node"):
	#if i >=len(D_Fem):
		#break
	#else:
		#old_txt =  child.text 
		#print i
		#new_txt = ','.join(str(x) for x in D_Fem[i,1:])
		#child.text = new_txt
		#child.set('updated','yes')
		#i +=1
for i in D_Fem[:,0]:
	for child in root.findall(".//Geometry/Nodes/node[%d]"%i):
		if str(D_Fem[i,0]) == child.attrib["id"]:
			#print child.attrib['id'], child.text
			#old_txt =  child.text 
			new_txt = ','.join(str(x) for x in D_Fem[i,1:])
			child.text = new_txt
			child.set('updated','yes')
				
for i in D_Tib[:,0]:
	for child in root.findall(".//Geometry/Nodes/node[%d]"%i):
		if str(D_Tib[i,0]) == child.attrib["id"]:
			#print child.attrib['id'], child.text
			#old_txt =  child.text 
			new_txt = ','.join(str(x) for x in D_Tib[i,1:])
			child.text = new_txt
			child.set('updated','yes')
				
#print old_txt
#print new_txt


	
## 6) update the Feb file
#i = 0
#for child in root.findall(".//Nodes/node"):
	#if i >=len(D_Tib):
		#break
	#else:
		#old_txt =  child.text 
		#print i
		#new_txt = ','.join(str(x) for x in D_Tib[i,1:])
		#child.text = new_txt
		#child.set('updated','yes')
		#i +=1
#print old_txt
#print new_txt
# 7) 
	

#temp01 = []
#for i in xrange(len(temp)):
	#for j in xrange(len(temp)):
		#if temp[i][2:5] !=temp[j][2:5]:
			#temp01.append(temp[i])
#temp01 = np.array(temp01)
#print temp01.shape
			
#for NodeSet in root[1].findall ('NodeSet'):
	##print NodeSet.attrib['name']
	#name = NodeSet.get('name')
	##print name
	#if name == 'Nodeset01_tib_cart_out':
		#NodalNum = NodeSet.attrib['id']
		#ArticularNodalNum.append(int(NodalNum))
		#print NodalNum
#ArticularNodalNum = np.array(ArticularNodalNum)
		
		#for child in NodeSet:
			#print child.attrib.values()

	

## extract nodal numbers from feb file

#for ch in root[2[2]:   
	#NodalNum = ch.attrib['id']
	#ArticularNodalNum.append(int(NodalNum))
#ArticularNodalNum = np.array(ArticularNodalNum)

#for ch in root[1][3]:
	#NodalNum = ch.attrib['id']
	#SubBoneNodalNum.append(int(NodalNum))
#SubBoneNodalNum = np.array(SubBoneNodalNum)
#j = 0
#for Nodes in root[1].findall('Nodes'):
	#for node in Nodes:
		##print int(node.attrib ['id'])
		#if int(node.attrib['id']) == SubBoneNodalNum[j]:
			#print node.attrib['id']
			#t = node.text
			#t1 =[float(word) for word in t.split(',')]
			#SubBoneNodalCoords.append(t1)
			#j +=1
			#print j
			
#SubBoneNodalCoords = np.array(SubBoneNodalCoords)
#np.save('/people/mkaz010/MAP_Project/MAP/Data/Paper_1_WFs/Hex8FromNodals/outputs/OKnee_tibia_SubBoneNodalCoords.npy', SubBoneNodalCoords)
##np.savetxt('/people/mkaz010/MAP_Project/MAP/Data/Paper_1_WFs/Hex8FromNodals/outputs/OKnee_tibia_SubBoneNodalCoords.xyz', SubBoneNodalCoords)

#print SubBoneNodalCoords.shape		
#print SubBoneNodalNum.shape 
#print ArticularNodalNum.shape 

#temp1 = []
#temp2 = []
#NewCoords = []
#for i in xrange (len(SubBoneNodalNum.T)):
	#a = [SubBoneNodalNum[i], tibCartInCoords[i][0],tibCartInCoords[i][1],tibCartInCoords[i][2]]
	#temp1.append(a)
#temp1 = np.array(temp1)
#print temp1.shape

#for i in xrange (len(ArticularNodalNum.T)):
	#a = [ArticularNodalNum[i], tibCartOutCoords[i][0],tibCartOutCoords[i][1],tibCartOutCoords[i][2]]
	#temp2.append(a)
#temp2 = np.array(temp2)
#print temp2.shape

#NewCoords = np.concatenate ((temp1, temp2), axis = 0)
#temp3 = NewCoords[np.argsort(NewCoords[:,0]),:]
#print temp3.shape

#j = 0
#for Nodes in root.iter('Nodes'):
	#for node in Nodes:
		#print node.attrib.values()
		#print 'before', node.text
		#node.text = str(temp3[j][1:4]).strip('[]')
		#print 'after', node.text
		#j = j+1
		#print j
		
#ET.tostring(tree, encoding='UTF-8', xml_declaration=True)
		
tree.write('C:\Users\Mousa Kazemi\Phd_project\cartilage_growing\femtib_cart_hex8_01.feb', encoding='utf-8', xml_declaration=True)

#subprocess.call(['/people/mkaz010/preview-1.17.0/preview.lnx64', '/people/mkaz010/MAP_Project/MAP/Data/Paper_1_WFs/Hex8FromNodals/outputs/femurtib_cart_hex8.feb'])

##for child in root[1][1]: #  Geometry
	##print child.attrib['id']  # extract node numbers

##geom = root[2]

##nodes = geom[0]
##node_numbers = [int(n.attrib['id']) for n in nodes]
##node_coords = [[float(x) for x in n.text.split(',')] for n in nodes]

##node_set_1 = geom[17]
##node_set_1_name = node_set_1.attrib['name']
##node_set_1_nodes = [int(n.attrib['id']) for n in node_set_1]
##node_set_1_coords = [node_coords[i] for i in node_set_1_nodes]

##node_set_2 = geom[18]
##node_set_2_name = node_set_2.attrib['name']
##node_set_2_nodes = [int(n.attrib['id']) for n in node_set_2]
##node_set_2_coords = [node_coords[i] for i in node_set_2_nodes]
