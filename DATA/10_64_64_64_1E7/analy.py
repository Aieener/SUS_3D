# analy.py
# A python program to analyze the SUS weighting function in order to reach the following goals:
# 1. plot the weight function
# 2. generate the normalized distribution for Z=1
# 3. extrapolate the N distribution for different Zs given by the user.
# Author: Yuding Ai
# Date: 2015 Oct 23

import math
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def PN():
	WF = [] # a list of my target Weighting function
	PN = [] # a list of number distribution

	with open("SUSWeight_function.txt","r") as file:
		for line in file:
			words = line.split()
			n = float(words[0]) #take the value
			WF.append(n); #append value into my WF list

	maxi = max(WF)
	if maxi > 500:
		for i in range(len(WF)):
			WF[i] = WF[i]-maxi +500
			PN.append(math.exp(WF[i]));

	PN = [float(i)/sum(PN) for i in PN]
	return WF,PN


def Pplot(PN,z):
	fig = plt.figure()	
	plt.plot(PN,'+b',markersize=3)
	Z = str(z)
	ylabel = 'P(N;Z='+ Z + ')'
	plt.ylabel(ylabel)
	plt.xlabel('N')
	title = 'P(N;Z='+ Z + ').png'
	fig.savefig(title, dpi=300, bbox_inches='tight')

def enlargePplot(PN,z):
	fig = plt.figure()	
	plt.plot(PN,'+b-',markersize=3,linewidth = 0.1)
	plt.xlim(8600,9600)
	plt.ylim(0,0.007)
	Z = str(z)
	ylabel = 'P(N;Z='+ Z + ')'
	plt.ylabel(ylabel)
	plt.xlabel('N')
	title = 'ENLP(N;Z='+ Z + ').png'
	fig.savefig(title, dpi=300, bbox_inches='tight')

def Wplot(WN):
	fig = plt.figure()	
	plt.plot(WN,'+r',markersize=1,)
	plt.ylabel('Weighting Function')
	plt.xlabel('N')
	title = 'WeightingFunc.png'
	fig.savefig(title, dpi=300, bbox_inches='tight')


def exploPN(W,z):
	P = [] # a list of number distribution
	for i in range(len(W)):
		W[i] = W[i] + i*math.log(z)

	maxi = max(W)
	if maxi > 500:
		for j in range(len(W)):
			W[j] = W[j]-maxi +500
			P.append(math.exp(W[j]));
	P = [float(k)/sum(P) for k in P]
	return P


def main():


	P = PN()[1] # take the P(N;z=1)
	W = PN()[0] # take the original weighting function 

	# Wplot(W)
	# Pplot(P,"1")
	# Pe = exploPN(W,4.44)
	# enlargePplot(Pe,4.44)

	for i in range(10):
		W = PN()[0] # take the original weighting function 		
		t = 3.83 + 0.02*i
		Pe = exploPN(W,t)
		# Pplot(Pe,t)
		enlargePplot(Pe,t)

main()
