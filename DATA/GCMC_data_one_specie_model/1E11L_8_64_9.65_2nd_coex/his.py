#Analysis # distribution for the 3-D Rods
#Author: Yuding Ai
#Date: 2015 July 28

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import math
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def his():
	N1 = [] # Ver
	N2 = [] # Hor
	N3 = [] # Up
	N = [] # Total number

	with open("dataplot.dat","r") as file:
		for line in file:
			words = line.split()

			n1 = float(words[2]) # Ver
			n2 = float(words[3]) # Hor
			n3 = float(words[4]) # Up
			ntot = n1+n2+n3
			N1.append(n1);
			N2.append(n2);
			N3.append(n3);
			N.append(ntot);


	fig1 = plt.figure()
	fig2 = plt.figure()
	fig3 = plt.figure()
	fig4 = plt.figure()
	fig5 = plt.figure()
	fig6 = plt.figure()
	ax1 = fig1.add_subplot(111)
	ax2 = fig2.add_subplot(111)
	ax3 = fig3.add_subplot(111)
	ax4 = fig4.add_subplot(111)
	ax5 = fig5.add_subplot(111)
	ax6 = fig6.add_subplot(111)
	numBins = 100
	# if N == "N1":
	# 	ax.hist(N1,numBins,color = 'blue', alpha = 0.8)
	# 	title = 'N1_#distribution.png'
	# 	fig.savefig(title, dpi=180, bbox_inches='tight')
	# elif N == "N2":
	# 	ax.hist(N2,numBins,color = 'red', alpha = 0.8)
	# 	title = 'N2_#distribution.png'
	# 	fig.savefig(title, dpi=180, bbox_inches='tight')
	# elif N == "N3":
	# 	ax.hist(N3,numBins,color = 'green', alpha = 0.8)
	# 	title = 'N3_#distribution.png'
	# 	fig.savefig(title, dpi=180, bbox_inches='tight')
	# else:
	# 	ax.hist(N1,numBins,color = 'blue', alpha = 0.6)
	# 	ax.hist(N2,numBins,color = 'red', alpha = 0.6)
	# 	ax.hist(N3,numBins,color = 'green', alpha = 0.6)
	# 	title = 'All_#distribution.png'
	# 	fig.savefig(title, dpi=720, bbox_inches='tight')

	ax1.set_title("Number Distribution for Vertical Rods")
	ax1.set_xlabel('Numbers')
	ax1.set_ylabel('Frequency')
	ax1.hist(N1,numBins,color = 'blue', alpha = 0.8, label='Vertical Rods')
	leg = ax1.legend()
	leg.get_frame().set_alpha(0.5)
	title = 'N1_#distribution.png'
	fig1.savefig(title, dpi=180, bbox_inches='tight')

	ax2.set_title("Number Distribution for Horizontal Rods")
	ax2.set_xlabel('Numbers')
	ax2.set_ylabel('Frequency')
	ax2.hist(N2,numBins,color = 'red', alpha = 0.8,label ='Horizontal Rods')
	leg = ax2.legend()
	leg.get_frame().set_alpha(0.5)
	title = 'N2_#distribution.png'
	fig2.savefig(title, dpi=180, bbox_inches='tight')

	ax3.set_title("Number Distribution for Up Rods")
	ax3.set_xlabel('Numbers')
	ax3.set_ylabel('Frequency')
	ax3.hist(N3,numBins,color = 'green', alpha = 0.8, label = 'Up Rods')
	leg = ax3.legend()
	leg.get_frame().set_alpha(0.5)
	title = 'N3_#distribution.png'
	fig3.savefig(title, dpi=180, bbox_inches='tight')

	ax4.set_title("Number Distribution for All")
	ax4.set_xlabel('Numbers')
	ax4.set_ylabel('Frequency')
	ax4.hist(N1,numBins,color = 'blue', alpha = 0.6,label = 'Vertical Rods')
	ax4.hist(N2,numBins,color = 'red', alpha = 0.6,label = 'Horizontal Rods')
	ax4.hist(N3,numBins,color = 'green', alpha = 0.6,label = 'Up Rods')
	leg = ax4.legend()
	leg.get_frame().set_alpha(0.5)
	title = 'All_#distribution.png'
	fig4.savefig(title, dpi=180, bbox_inches='tight')

	ax5.set_title("Total Number Distribution")
	ax5.set_xlabel('Numbers')
	ax5.set_ylabel('Frequency')
	ax5.hist(N,numBins,color = 'yellow', alpha = 0.8, label = 'Total Rods')
	leg = ax5.legend()
	leg.get_frame().set_alpha(0.5)
	title = 'Ntot_#distribution.png'
	fig5.savefig(title, dpi=180, bbox_inches='tight')

	ax6.set_title("Log (Total Number Distribution)")
	ax6.set_xlabel('Numbers')
	ax6.set_ylabel('Log(Frequency)')
	ax6.hist(N,numBins,color = 'pink', alpha = 0.8, label = 'Log (Total Rods)',log =True)
	leg = ax6.legend()
	leg.get_frame().set_alpha(0.5)
	title = 'LogNtot_#distribution.png'
	fig6.savefig(title, dpi=180, bbox_inches='tight')



his()



# def main():
# 	print "#=====================================#"
# 	print "# Welcome to the Distribution factory #"
# 	print "#=====================================#"

# 	check = True
# 	while (check):
# 		N = raw_input("Please tell me which specie distribution do you want for this time? \nType 'N1','N2','N3' or 'all'\n")
# 		his(N)
# 		again = raw_input("Do you want to checkout other distribution? y/n\n")
# 		if again == "n":
# 			check = False
# 	print "Thanks for trying on our Drawing tool!\nSee you soon! LOL"


# main()