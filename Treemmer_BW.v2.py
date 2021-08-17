#Treemmer

#   This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# You have to install ete3 http://etetoolkit.org/
# and joblib https://pythonhosted.org/joblib/ to run Treemmer

# Original code written by Fabrizio Menardo.
# Edited by Bryan Wee to maximise phenotypic diversity for GWAS analyses - 20180308
# v2 - 20180313 - Bug fixes - a contraPair sometimes gets deleted with each iteration
#	Included new function get_dmin()

from joblib import Parallel, delayed
from ete3 import Tree
import sys
import random
import operator
import argparse



############################################################			define arg type float 0 < X > 1		###############################################################

def restricted_float(x):
	x = float(x)
	if x < 0.0 or x > 1.0:
		raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
	return x

##########################################			FIND LEAVES NEIGHBORS OF A LEAF (2 NODE OF DISTANCE MAX) and calc DISTANCE			#######################

def find_N(t,leaf):
	dlist ={}
	parent= leaf.up
	dist_parent=leaf.dist
	flag=0

	if arguments.verbose==3:
		print "leaf findN at ieration:	" + str(counter)
		print leaf
		print "parent findN at ieration:	" + str(counter)
		print parent
		print parent.get_children()

	sister_flag=0
	for n in range(0,len(parent.get_children())):				##this for loop start from parent and climb up max two nodes, if it finds leaves calculate the distances,
		if parent.is_root():
			flag=1
			break
		if arguments.verbose==3:
			print "children	" + str(n)
			print parent.children[n]

		if (parent.children[n].is_leaf()):						# search at one node of distance
			if (parent.children[n] != leaf):
				DIS = leaf.get_distance(parent.children[n])
				dlist.update({leaf.name + "," +parent.children[n].name : DIS})
				flag=flag+1
				if arguments.verbose==3:
					print leaf.name + "," +parent.children[n].name + str(DIS) + "have one node of distance"
		else:
			if flag == 0:
				if arguments.verbose==3:					#going up, search at two nodes of distance
					print "going up, brother is node"

				temp_dlist={}
				for nn in range(0,len(parent.children[n].get_children())):
					if (parent.children[n].children[nn].is_leaf()):
						DIS = leaf.get_distance(parent.children[n].children[nn])
						temp_dlist.update({leaf.name + "," +parent.children[n].children[nn].name : DIS})
						sister_flag=sister_flag +1


	if ((sister_flag==1) and (flag==0)):						#collect results at two nodes of distance only if there are no leaves that are closer
		dlist.update(temp_dlist)
		if arguments.verbose==3:
			print str(temp_dlist) + "	are not sister taxa, but neighbours first is leaf, second is upper neighbor"



	if (flag == 0):			#### this means that the leaf has no neighbors at one node of dist
		parent=parent.up 		#### therefore I climb the tree down towards the root of one more step and look for leaves
		multi_flag=0

		if arguments.verbose==3:
			print "going down"
			print "gran parent"
			print parent
		temp_dlist={}
		for n in range(0,len(parent.get_children())):		#this for loop start from gran parent and climb up max one nodes, if it finds leaves calculate the distances,
			if parent.is_root():
				break
			if (parent.children[n].is_leaf()):
				DIS = leaf.get_distance(parent.children[n])
				multi_flag = multi_flag+1
				temp_dlist.update({leaf.name + "," +parent.children[n].name : DIS})
		if multi_flag==1:					# this is to deal with polytomies
			dlist.update(temp_dlist)
			if arguments.verbose==3:
				print leaf.name + "," +parent.children[n].name + str(DIS) + "	are not sister taxa, but neighbours first is leaf, second is neighbor of downstair (towards root)"

	#print dlist #BW: print out distances between pairs of leaves
	return dlist

##########################################		BW: GET DMIN VALUES AFTER REMOVING CONTRAPAIRS			#######################
	#BW: contraPair - a pair of taxa that differ in phenotype

def get_dmin(dlist, list_IN): #pulls out list of pairs that are not CONTRAPAIRS
	#Remove contrapairs from last round
	for contraPair in list_IN: #This loop removes all the contraPairs found so far before. This is so that the minimum d value does not take into account the contraPairs
		#print "Keeping contrapair: " + contraPair
		#print "Checking size of dlist before deleting contrapair: " + str(len(dlist)) #BW Working
		dlist.pop(contraPair, None) #delete pair from dlist if in list_IN
		#print "Checking size of dlist after deleting contrapair: " + str(len(dlist)) #BW Working
	min_val = min(dlist.itervalues()) #BW: i think min_value is the minimum distance bewteen two pairs

	d_min = {}
	pairs_to_remove_from_d_min = []
	for k, v in dlist.iteritems(): #BW: iterates over all paris (k) and v is the distance between the two leaves
		#print str(k) + ' -- ' + str(v) + "    TESTING" #BW test
		if v == min_val:
			d_min.update({k:v})

	#BW: This loop strips out pairs with the same phenotype
	#BW: d_min are the number of pairs that share the minimum distance.
	for key in d_min:
		key1 = [key2.strip("''") for key2 in key.split(",")] #BW this splits the pair names and then removes quotation marks from the taxa names
		if key1[0][-3:] == key1[1][-3:]: #This ifelse loop checks if the last 3 characters of the taxon name matches.
			next#BW if suffix matches (i.e. same ), leave in d_min
		else:
			pairs_to_remove_from_d_min.append(key)

	for i in pairs_to_remove_from_d_min:
		#print "Before " + str(len(d_min))
		del d_min[i]
		#print "After " + str(len(d_min))

	return (d_min,pairs_to_remove_from_d_min)

##########################################		IDENTIFY  LEAF TO PRUNE			#######################

def find_leaf_to_prune(dlist,list_IN):					#parse the list with all neighbor pairs and distances, find the closest pair and select the leaf

	d_min={}

	#BW This function was added to get the list of pairs with the smallest distance between, after taking out pairs with contrasting phenotypes
	#Sometimes, the all the pairs with the minimum distances are contrapairs. Therefore an empty dictionary is returned which does not allow the next function to prune a random pair. To prevent this, everytime d_min is returned empty, it gets fed back into the get_dmin() function to get the new list pairs with mininimum distance values.
	while not d_min:
		(d_min,list_IN) = get_dmin(dlist,list_IN)

	list_OUT = list_IN

	min_val = min(d_min.itervalues()) #BW: min_value is the minimum distance bewteen two pairs
	print "Min value in the list of pairs: " + str(min_val)
	find_leaf_to_prune.min_val = min_val #define min_val to call outside this function
	#END BW EDIT

	pair = str(random.choice(list(d_min)))
	#BW TODO: It might be possible to delete all of the pairs with minimum tree values in each iteration instead of just deleting one per iteration
	print "Deleting pair: " + pair #Test

	if arguments.prune_random: #BW: Prune random leaf flag used, then
		pair= str(random.choice(list(dlist)))
	pair=pair.split(",")
	leaf1 = t.search_nodes(name=pair[0])[0]
	leaf2 = t.search_nodes(name=pair[1])[0]

	if  (leaf1.dist > leaf2.dist):
		if (arguments.leaves_pair == 1):
			leaf_to_prune = leaf2.name
			dist = leaf2.dist
		if (arguments.leaves_pair == 0):
			leaf_to_prune = leaf1.name
			dist = leaf1.dist

	if  (leaf1.dist < leaf2.dist):
		if (arguments.leaves_pair == 1):
			leaf_to_prune = leaf1.name
			dist = leaf1.dist
		if (arguments.leaves_pair == 0):
			leaf_to_prune = leaf2.name
			dist = leaf2.dist

	if  ((leaf1.dist == leaf2.dist) or (arguments.leaves_pair ==2)):
		leaf_to_prune = random.choice(list(pair))			#this select the leaf at random within the couple
		dist = leaf1.dist

	#print "list at end of loop: , " + str(list_OUT) #BW: Check the list


	return (leaf_to_prune,dist,list_OUT)

##########################################				PRUNE LEAF FROM TREE			#######################

def prune_t(leaf_to_prune,tree):

	G = tree.search_nodes(name=leaf_to_prune)[0]
	parent= G.up
	dist_parent=G.dist

	if (len(parent.get_children()) == 2):


		if  parent.children[0] != G:
			parent.children[0].dist = parent.children[0].dist + parent.dist

		if parent.children[1] != G:
			parent.children[1].dist = parent.children[1].dist + parent.dist

	G.detach()

	if (len(parent.get_children()) == 1):
		parent.delete()		# after pruning the remaining branch will be like this ---/---leaf_name. I delete useless node keeping the b length


	return tree

####################################################################	calculate Tree length ##########################################################3

def calculate_TL(t):
	tree_length=0
	for n in t.traverse():
		tree_length=tree_length+n.dist
	tot_TL = tree_length
	return(tot_TL)

##########################################		PRUNE LEAF FROM MATRIX		#######################

def prune_dist_matrix(dlist,leaf_to_prune):
	key_del=[]
	for k, v in dlist.iteritems():

		(one,two)=k.split(",")
		if ((one == leaf_to_prune) or (two == leaf_to_prune)):
			key_del.append(k)

	for KK in key_del:
		del dlist[KK]
	return dlist

##########################################		parallel loop		#######################

def parallel_loop(i):
	n=i
	while n < len(leaves):
		N_list=find_N(t,leaves[n])
		n=n+arguments.cpu    			#n of  threads
		if N_list:
			DLIST.update(N_list)
	return DLIST

##########################################		write output with stop option		#######################

def write_stop(t,output1,output2):
	F=open(output1,"w")
	F.write(t.write())
	F.close()
	leaves = t.get_leaves()
	list_names=[]
	for leaf in leaves:
		list_names.append(leaf.name)
	F=open(output2,"w")
	F.write("\n".join(list_names))
	F.close()


######   SOFTWARE STARTS

parser = argparse.ArgumentParser()

parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-X','--stop_at_X_leaves', metavar='0-n_leaves', default='0', help='stop pruning when the number of leaves =  X', type =int, nargs='?')
parser.add_argument('-MDIST','--stop_at_MDIST', metavar='0-1000', default='0', help='stop pruning when minimum distance between any pair =  MDIST', type =float, nargs='?') #BW added this option
parser.add_argument('-RTL','--stop_at_RTL', metavar='0-1', default='0', help='stop pruning when the relative tree length falls below RTL', type =restricted_float,nargs='?')
parser.add_argument('-r','--resolution', metavar='INT', default=1,help='number of leaves to prune at each iteration (default: 1)',type =int, nargs='?')
parser.add_argument('-p','--solve_polytomies',help='resolve polytomies at random (default: FALSE)',action='store_true',default =False)
parser.add_argument('-pr','--prune_random',help='prune random leaves (default: FALSE)',action='store_true',default =False)
parser.add_argument('-lp','--leaves_pair', metavar='0,1,2', default=2,help='After the pair of leaves with the smallest distance is dentified Treemmer prunes: 0: the longest leaf\n1: the shortest leaf\n2: random choice (default: 2)',type =int, nargs='?')
parser.add_argument('-np','--no_plot',help='do not load matplotlib and plot (default: FALSE)',action='store_true',default =False)
parser.add_argument('-fp','--fine_plot',help='when --resolution > 1, plot RTL vs n leaves every time a leaf is pruned  (default: FALSE => plot every X leaves (X = -r))',action='store_true',default =False)
parser.add_argument('-c','--cpu', metavar='INT', default=1,help='number of cpu to use (default: 1)',type =int, nargs='?')
parser.add_argument('-v' ,'--verbose', metavar='0,1,2', default='1', help='0: silent, 1: show progress, 2: print tree at each iteration, 3: only for testing (findN), 4: only for testing (prune_t) (default: 1)', type =int, nargs='?',choices=[0,1,2,3,4])


arguments = parser.parse_args()


if ((arguments.stop_at_RTL > 0) and (arguments.stop_at_X_leaves > 0)):
	raise argparse.ArgumentTypeError("-X and -RTL are mutually exclusive options")

t = Tree(arguments.INFILE,format=1)

counter =0
output=[]
stop=0

TOT_TL=calculate_TL(t)
TL=TOT_TL
ori_length = len(t)
if arguments.solve_polytomies:
	t.resolve_polytomy()

if arguments.verbose > 0:													# print progress on standard output
	print "N of taxa in tree is : "+ str(len(t))

	if arguments.solve_polytomies:
		print "\nPolytomies will be solved at random"
	else:
		print "\nPolytomies will be kept"
	if arguments.prune_random:
		print "\nA random leaf is pruned at each iteration, you don't really need Treemmer to do this"
	if arguments.stop_at_X_leaves:
		print "\nTreemmer will reduce the tree to" + str(arguments.stop_at_X_leaves) + " leaves"
	if arguments.stop_at_MDIST:
		print "\nTreemmer will reduce the tree until the minimum dist between two leaves with non contrasting pairs is " + str(arguments.stop_at_MDIST)
	else:
		if arguments.stop_at_RTL:
			print "\nTreemmer will reduce the tree to" + str(arguments.stop_at_RTL) + " of the original tree length"
		else:
			print "\nTreemmer will calculate the tree length decay"

	print "\nTreemmer will prune " + str(arguments.resolution) + " leaves at each iteration"
	print "\nTreemmer will use " + str(arguments.cpu) + " cpu(s)"
x=[]
y=[]


output.append ('1	' + str(len(t))) 						#append first point to the output with RTL = 1 (before starting pruning)################################

length=len(t)
x.append(length)
y.append(1)

list_IN = [] #BW: Initialise list_IN


while (len(t) > 3):								#################### Main loop ################################
	counter = counter +1
	leaves = t.get_leaves()
	DLIST={}

	if arguments.verbose > 0:
		print "\niter		" + str(counter)
	if arguments.verbose > 1:
		print "calculating distances"

	DLIST = Parallel(n_jobs=arguments.cpu)(delayed(parallel_loop)(i) for i in range(0,arguments.cpu))
	result = {}

	for d in DLIST:								#when running in parallel DLIST is updated in a weird way, it is a dict of dicts, this for loop merge them all in one
		result.update(d)

	DLIST=result

	if arguments.verbose > 1:
		print DLIST
		print "\npruning big deal\n"

	for r in range (1,arguments.resolution+1):

		if ((len(DLIST)==0) or (len(t)<4)):
			break

		(leaf_to_p,dist,list_OUT) = find_leaf_to_prune(DLIST,list_IN) #BW: leaf_to_p is the leaf/taxon that will be pruned
		leaf_to_prune = t.search_nodes(name=leaf_to_p)[0]
		t = prune_t(leaf_to_p,t)
		TL= calculate_TL(t)
		DLIST=prune_dist_matrix(DLIST,leaf_to_p)
		rel_TL=TL/TOT_TL

		list_IN = list_OUT #BW

		#################################  		OUTPUT 		##########################################################

		if (arguments.fine_plot):										# plot point in rtld after every leaf independently of -r
			output.append (str(rel_TL) + '	' + str(len(t)))
			length=len(t)
			x.append(length)
			y.append(rel_TL)

		if arguments.stop_at_X_leaves:										# if stop criterium is met (X) ==> output
			if arguments.stop_at_X_leaves >= len(t):
				output1=arguments.INFILE+"_trimmed_tree_X_" + str(arguments.stop_at_X_leaves)
				output2=arguments.INFILE+"_trimmed_list_X_" + str(arguments.stop_at_X_leaves)
				write_stop(t,output1,output2)
				stop=1
				break

		if arguments.stop_at_MDIST:										# BW added: if stop criterium is met min dist ance is > M
			print "Min distance value argument is : " + str(arguments.stop_at_MDIST)
			print "Min distance value is: " + str(find_leaf_to_prune.min_val)
			if arguments.stop_at_MDIST <=  find_leaf_to_prune.min_val: #####BW get min dist value from somewhere len(t):
				output1=arguments.INFILE+"_trimmed_tree_MDIST_" + str(arguments.stop_at_MDIST)
				output2=arguments.INFILE+"_trimmed_list_MDIST_" + str(arguments.stop_at_MDIST)
				write_stop(t,output1,output2)
				stop=1
				break


		if arguments.stop_at_RTL:										# if stop criterium is met (RTL) ==> output
			if arguments.stop_at_RTL >= rel_TL:
				output1=arguments.INFILE+"_trimmed_tree_RTL_" + str(arguments.stop_at_RTL)
				output2=arguments.INFILE+"_trimmed_list_RTL_" + str(arguments.stop_at_RTL)
				write_stop(t,output1,output2)
				stop=1
				break

		if arguments.verbose > 1:										# print progress on standard output

			print "\n ITERATION RESOLUTION:	" + str(r)
			print "leaf to prune:\n" + str(leaf_to_p) + "	" + str(dist)
			print "\n new tree"
			print t
			print "\nRTL :	" + str(rel_TL) + " N_seq:	" +str(len(t))
			print "\nnew matrix\n"
			print DLIST

	if (stop ==1):
		print "\nRTL :	" + str(rel_TL) + " N_seq:	" +str(len(t))
		break

	if not (arguments.fine_plot):											# normal plot (with -fp = FALSE)
		output.append (str(rel_TL) + '	' + str(len(t)))
		length=len(t)
		x.append(length)
		y.append(rel_TL)


	if arguments.verbose==1:
		print "\nRTL :	" + str(rel_TL) + " N_seq:	" +str(len(t)) + " Min_value: " + str(find_leaf_to_prune.min_val)



if stop == 0:														# create file for plot of rltd
	F=open(arguments.INFILE+"_res_"+ str(arguments.resolution) + "_LD","w")
	F.write("\n".join(output))


#################################################  make plot  directly ##########################################################
	if not arguments.no_plot:
		import numpy as np
		import matplotlib.pyplot as plt
		from matplotlib.ticker import MaxNLocator
		ax = plt.figure().gca()
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		plt.scatter(x, y, s= 2, c= 'black')
		plt.xlim(ori_length,0)
		plt.ylim(-0.02,1.02)
		plt.xlabel('Number of leaves')
		plt.ylabel('Relative tree length')
		#plt.savefig(arguments.INFILE+'_res_'+ str(arguments.resolution)+'_TLD.png')
		plt.savefig(arguments.INFILE+'_res_'+ str(arguments.resolution)+'_TLD.pdf')
