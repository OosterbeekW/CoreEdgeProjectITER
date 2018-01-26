#!/usr/bin/env python

###########################################################################
# Python script to construct an edge stability diagram
#
# Guido Huijsmans, ITER Organization, 11/11/2011
#
# - each submitted job calculates all toroidal mode numbers
#
############################################################################
#
# Edited Wouter Oosterbeek --- 20/10/2017
#
###########################################################################
#
# Specific Construction for file locations:
#
#  main (called jobdir in code)
#   run.py (this script)	
#   HELENA (exefile and script to run - hel12a runhel)
#     namelist (with the inputfile of helena)	
#     output (the outputfile of helena goes here)
#   MISHKA (exefile and script to run - mishka1_fast runmis)
#     namelist (inputfile of mishka)
#     output (outputfile of mishka goes here)	  
#   hel_...mis_...(created) - called savedir in code: 
#				inputfiles in name + possible date
#	case_(created)  - casedir, file for every npcore,nbped,ntor  
#	   HELENA_ouput (created) Outputfiles of HELENA copied here 
#	   MISHKA_output (created) Outputfiles of MISHKA copied here	
#	datafile.txt(created) - with numbers->values	
#
##############################################################################


import matplotlib
matplotlib.use('Agg')
import random
import os
import sys
import shutil
import re
import fileinput
import getpass
import time
from datetime import datetime
import numpy
import matplotlib.pyplot as plt
import pylab
import Tkinter as tk
from matplotlib.backends.backend_pdf import PdfPages
import subprocess 


###########################################################
# Function definitions
###########################################################

def replace_parameter(file_name,var,value):
  for line in fileinput.FileInput(file_name, inplace=1):
    print re.sub(var+"[ ]*=[ ]*[0-9.\+\-E]*[ ]*",var+" = "+str(value),line),

def replace_name(file_name,old_var,new_var):
  for line in fileinput.FileInput(file_name, inplace=1):
    print line.replace(old_var,new_var),

def find_B(file_name):
  for all_lines in fileinput.FileInput(file_name):
    for line in all_lines.split("\n"):
      if re.search("A,B,C" , line):
        abc_str = line.split(":")
        B = float(abc_str[1].strip().split(" ")[2])
  return B	
 
def find_normbeta(file_name):
  for all_lines in fileinput.FileInput(file_name):
    for line in all_lines.split("\n"):
      if re.search("NORM. BETA" , line):
        abc_str = line.split(":")
        B = float(abc_str[1].strip())
  return B	

def find_Betapedtot(file_name):
  for all_lines in fileinput.FileInput(file_name):
    for line in all_lines.split("\n"):
      if re.search("POLOIDAL BETA" , line):
        abc_str = line.split(":")
        Betapedtot = float(abc_str[1].strip())
  return Betapedtot	
  
def make_structure():
	try:	
		os.mkdir("HELENA")
		os.mkdir("HELENA/namelist")
		os.mkdir("HELENA/mapping")
		os.mkdir("HELENA/output")
		os.mkdir("HELENA/plot")
		os.mkdir("HELENA/allout")
		os.mkdir("MISHKA")
		os.mkdir("MISHKA/namelist")
		os.mkdir("MISHKA/output")
		os.mkdir("MISHKA/plot")
		os.mkdir("MISHKA/allout")
		os.mkdir("Results")
	except OSError:
		pass

def launch_case(job_dir,case_dir,bped,pcore,Ntor,helrun,helfile,misfile):
# Replace betap_ped and p_core values in helena input file
	#time.sleep(1) #4 #Give time to job to finish copying all files
	
	os.chdir(job_dir+"/HELENA/namelist")
	replace_parameter(helfile,"BETAP_PED",bped)
	replace_parameter(helfile,"COREP",pcore)
	# A value of B has to be chosen so that HELENA converges over the range of bped/pcore values
	# (4.5*pcore**3-7.9*pcore**2+4.7*pcore-0.35): negative for x<0.09 ; Sometimes B too large for large pcore
	# (1.24*pcore**3-3.8*pcore**2+3.9*pcore-0.3): negative for x<0.09
	# (0.06*pcore**3-0.22*pcore**2+0.31*pcore-0.01)
	replace_parameter(helfile, " B", (0.06*pcore**3-0.22*pcore**2+0.31*pcore-0.01)) 
	#FCUR is chosen to adapt the bootstrap current to the pedestal pressure
	#helena input file ITERhelcur: (3*bped**2+5.29*bped-0.0295)
	replace_parameter(helfile,"FCUR", (12*bped-0.57) )
	#print("Core pressure and pedestal beta parameters: " + str(pcore) + " " + str(bped))
	
	os.chdir(job_dir+"/MISHKA/namelist")

	Ntor = Ntor*-1
	replace_parameter(misfile,"NTOR",Ntor)
	

	run_dir = "/tmp/"+case_dir
	final_dir = "hel_"+helrun+"_"+helfile+"_mis_"+misfile+"__"+dt
	jobname = "CoreEdge_"+dt+"_b"+str(bped)+"p"+str(pcore)+"n"+str(Ntor)
# Change names in general PBS file to case-specific names
	os.chdir(job_dir)	
	os.system("cp job_master job")
	replace_name("job","workdir",run_dir)
	replace_name("job","tbc_helfile",helfile)
	replace_name("job","tbc_misfile",misfile)
	replace_name("job","tbccase",case_dir)
	replace_name("job","finaldir",final_dir)
	replace_name("job","tbcname",jobname)

	print("Name: "+jobname)
#Copy of the job PBS file to save directory
	os.system("cp job "+save_dir+"/"+case_dir+"/job")
#Copy changed files to case_dir folder (so that they can't be overwritten)
	os.system("cp HELENA/namelist/"+helfile+" "+save_dir+"/"+case_dir+"/"+helfile)
	os.system("cp MISHKA/namelist/"+misfile+" "+save_dir+"/"+case_dir+"/"+misfile)
	#os.system("pwd")
	#print("cp HELENA/namelist/"+helfile+" "+save_dir+"/"+case_dir+"/"+helfile)
	#print("cp HELENA/namelist/"+misfile+" "+save_dir+"/"+case_dir+"/"+misfile)
#Send job to be executed in cluster	
	#print(" ")
	os.chdir(job_dir+"/MISHKA/namelist")
	os.system("sed -n 8p "+misfile)
	os.chdir(job_dir)
	#os.system("df -h /tmp/")
	#os.system("qsub -q batch job")
	os.system("qsub  job")

	#time.sleep(4)
	
	print "submitted : ",ij,",",ip," ; ",bped,",",pcore,",",Ntor

def plot_realworld(job_dir,r,jpar,p0,jboot,bped,pcore,nbped,npcore,ip,ik): 
# Plot real profile, output from helena
	cmap = matplotlib.cm.gist_rainbow( float(ip*nbped+ij)/float(nbped*npcore) )
	cur_dir = os.getcwd()
	os.chdir(save_dir+"/case_bped"+str(ik)+"_pcore"+str(ip)+"_Ntor2/HELENA_output")
	
	plotfile = open('profilesout.TXT')

	r = []
	jpar = []
	p0 = []
	jboot = []

	for line in plotfile:
		data = line.split(" ")
		data = filter(None, data)
		r.append(float(data[0]))
		jpar.append(float(data[1]))
		p0.append(float(data[2]))
		data[3] = data[3].strip()
		jboot.append(float(data[3]))	

	plotlbl = "b"+str(bped)+"p"+str(pcore)
	plt.figure(1)
	plt.plot(r,jpar,label = plotlbl, color=cmap)
	plt.figure(2)
	plt.plot(r,p0,label = plotlbl, color=cmap)
	plt.figure(3)
	plt.plot(r,jboot,label = plotlbl, color=cmap)


	os.chdir(cur_dir)

###########################################################
start_time = time.time()
copytime = 0

os.system("cp /home/ITER/oosterw/CoreEdgeStability/runall.py /home/ITER/oosterw/CoreEdgeStability/runall.save" )

# Choose helena script to run and helena,mishka input files
helrun = "helwo"
helfile = "ITERhelcur" #"ITERhelcur"
misfile = "test" #"test"
numplot = 1 # plotting real data profiles: 1 is on, otherwise (0) is off
arch = 1 # copy files to ~/Archive when set to 1

# Number of points for both variables
nbped = 10 #21 #21 #6 #13 #19  # Every loop (one value for corep) submits nbped*len(Ntor) jobs   
npcore = 23 #31 #31 #41   
# toroidal mode numbers to be checked
Ntor = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32] #Should be integers (not eg 2.0) #[2,3,5,8,11,14,15,17,20,23,25,26,29,30,32] #[30] #[2,3,5,10,20,30]
print("Toroidal mode numbers: "+str(Ntor))
# Current directory as job_dir 
job_dir = os.path.join(os.getcwd()) 

print("\n job_dir  : " + job_dir)
print(str(npcore)+" pcore loops, "+str(nbped)+" bped values, with "+str(len(Ntor))+" toroidal modes")
# Starting values of both variables and steps to make
bped_start = 0.05 #0.05 #0.08 #beta_poloidal = 0.5562
bped_step   = 0.005 #0.005
pcore_start = 0.1 #[round(elem,2) for elem in numpy.linspace(4.0,1,nbped)*0.1] #0.1
pcore_step =  [0.09, 0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05] #numpy.round( numpy.array([5.9,5.2,4.5,4.0,3.6,3.3,2.9])*1.3*(1.0/(npcore-1))*0.4 , 2 )
# For bped 0.05 - 0.005 : numpy.round(numpy.array([5.9,5.2,4.5,4.0,3.6,3.3,2.9,2.7,2.5,2.3,2.1,2.0,1.85,1.7,1.6,1.5,1.4,1.3,1.25,1.17,1.1])*1.3*(1.0/(npcore-1)), 2 )
# for npcore=31 --> ([ 0.26,  0.23,  0.2 ,  0.17,  0.16,  0.14,  0.13,  0.12,  0.11,  0.1 ,  0.09,  0.09,  0.08,  0.07,  0.07,  0.06,  0.06,  0.06,  0.05,  0.05,  0.05])
### 		Until *1.1 all points equilibria exist, higher gives no existing equilibria for low beta_ped
#[round(elem,2) for elem in numpy.linspace(mult_fact_pcorestep,1,nbped)*1.5] #[0.7]*nbped #0.05 # # x,1 multiplication dependent on scan of parameters
print pcore_step
if len(pcore_step) != nbped:
	print("pcore_step size does not match number of bped steps")
	sys.exit()

dt = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
#try:
#	os.mkdir("betaped_"+bped_start+"_"+bped_end+"_pcore_"+pcore_start+"_"+pcore_end)
#except Exception, e: print(e)
#	#sys.exit(1)

os.chdir(job_dir+"/Results")
# Create direction in which results of this run are saved (unique name through datetime)
new_dir = "hel_"+helrun+"_"+helfile+"_mis_"+misfile+"__"+dt
os.mkdir(new_dir)
if arch == 1:
	os.mkdir("/home/ITER/oosterw/Archive/Results/"+new_dir)
save_dir = job_dir+"/Results/"+new_dir
print("New dir made at "+new_dir)

os.chdir(save_dir)
# Create two files to write data to: one full data, one only for most unstable cases
datafile = open("fulldata.txt",'w') 
uf_data = open("data.txt",'w')
datafile.write("HELENA: "+helrun+" "+helfile+"\n"+"MISHKA: "+misfile+"\n"+dt+"\n"+str(pcore_step)+"\n")
uf_data.write("HELENA: "+helrun+" "+helfile+"\n"+"MISHKA: "+misfile+"\n"+dt+"\n"+str(pcore_step)+"\n")
datalostr = "{0:>5}\t{1:>6}\t{2:>9}\t{3:>6}\t{4:>6}\t{5:>5}\t\t{6:>11} , {7:>11}\t\t{8:>15}\n"
datafile.write(datalostr.format("ibeta","ipcore","beta_ped","pcore","betaN","Ntor","Ev Real","Ev Imag","Relative change"))
uf_data.write(datalostr.format("ibeta","ipcore","beta_ped","pcore","betaN","Ntor","Ev Real","Ev Imag","Relative change"))
datalo = "{0:>5}\t{1:>6}\t{2:>9}\t{3:>6}\t{4:>6}\t{5:>5}\t\t{6:>+1.4e} , {7:>+1.4e}\t\t{8:>1.2e}"

# Make empty array to put values of stability plot in
x = []
y = []
z = [] #[None for i in range(npcore)] for j in range(nbped)
bN_plot = []
Ntor_plot = [] 
#After n loops (i is increased per loop), program goes to sleep: Dont submit too many jobs
njobs=55
i=1

# Set parameters for plots
if numplot==1:
	jparplot = plt.figure(1)
	presplot = plt.figure(2)
	jbootplot = plt.figure(3)
	#add1 = jparplot.subplot(111)
	#add2 = presplot.subplot(111)
	#add3 = jbootplot.subplot(111)
	#pdf = PdfPages("ThreeGraphs.pdf")
	r = []
	jpar = []
	p0 = []
	jboot = []

#Create variables to check change in eigenvalue	
old_stab = 0.0
old_stab2 = 0.0
for ip in range(npcore):       # loop over core pressure parameter
	i=0	
	for ij in range(nbped):	# loop over pedestal poloidal beta
		pcore = pcore_start + ip*pcore_step[ij]	
		if (pcore <= 0):
			print("Core pressure parameter smaller than zero")
			break

		bped = bped_start + ij*bped_step
		if (bped <= 0):
			print("Beta ped parameter smaller than zero")
			break
				
		for intor in range(len(Ntor)):
			os.chdir(save_dir)
# Create new folder for every bped-pcore combination
			case_dir = "case_bped"+str(ij)+"_pcore"+str(ip)+"_Ntor"+str(Ntor[intor])		
			try:
				os.mkdir(case_dir)
			except OSError:
				sys.exit(1)
# Create directories for output files
			os.chdir(case_dir)
			os.mkdir("HELENA_output")
			os.mkdir("MISHKA_output")			
# Run helena and mishka
			launch_case(job_dir,case_dir,bped,pcore,Ntor[intor],helrun,helfile,misfile)
# Check if not too many runs (too many jobs submitted)
			i = i + 1
			if (i%njobs)==0:
				jobname = "CoreEdge_"+dt+"_b"+str(bped)+"p"+str(pcore)+"n"+str(-1*Ntor[intor])
				print("--- Enough jobs submitted. Just relax and take it easy ;)")
				for remaining in range(60*60, 0, -1):
					check_completed = subprocess.check_output("qselect -u oosterw -N "+jobname+" -s C | wc -l", shell=True)
					check_queued = subprocess.check_output("qselect -u oosterw -N "+jobname+" -s Q | wc -l", shell=True)
					check_running = subprocess.check_output("qselect -u oosterw -N "+jobname+" -s R | wc -l", shell=True)
					if (int(check_queued) == 0) and (int(check_running) == 0):
						break
					else:
						sys.stdout.write("\r")
						sys.stdout.write("{:2d} seconds relaxing, while waiting for jobs to finish.".format(remaining))
						sys.stdout.flush()
						time.sleep(1)
				print("\n")
				
# #Hold on for some minutes, to let all qsub finish
	# for remaining in range(int(0.1*60), 0, -1):
		# sys.stdout.write("\r")
		# sys.stdout.write("{:2d} seconds remaining.".format(remaining))
		# sys.stdout.flush()
		# time.sleep(1)
	# print(" ")	

	stab_crit = 1.0
	print("\n\n")
# Another loop over bped and Ntor after all submitted jobs have finished
	for ik in range(nbped):
		pcore = pcore_start + ip*pcore_step[ik]	
		if (pcore <= 0):
			print("Core pressure parameter smaller than zero")
			break
		bped = bped_start + ik*bped_step
		if (bped <= 0):
			print("Beta ped parameter smaller than zero")
			break
			
		#Define start parameters for determining and checking eigenvalues from created files 'Eigenvalues.txt'
		ev0_old = 0.0	
		ev0 = None
		ev1 = None
		for intor in range(len(Ntor)):
# Wait for the job to be completed (not running or in queue) before trying to read files
			jobname = "CoreEdge_"+dt+"_b"+str(bped)+"p"+str(pcore)+"n"+str(-1*Ntor[intor])
			#print(jobname)
			for holdon in range(1,int(3600),1):
				check_completed = subprocess.check_output("qselect -u oosterw -N "+jobname+" -s C | wc -l", shell=True)
				check_queued = subprocess.check_output("qselect -u oosterw -N "+jobname+" -s Q | wc -l", shell=True)
				check_running = subprocess.check_output("qselect -u oosterw -N "+jobname+" -s R | wc -l", shell=True)

				if (int(check_completed) > 0) and (int(check_queued) == 0) and (int(check_running) == 0):
					sys.stdout.flush()
					print("\nJob "+jobname+" is completed. (b"+str(ik)+"p"+str(ip)+"n"+str(-1*Ntor[intor])+")") #If job is already done, this gives whiteline, if you had  to wait not...
					time.sleep(1)
					break	
				elif (int(check_completed) == 0) and (int(check_queued) == 0) and (int(check_running) == 0):
					sys.stdout.flush()
					print("\nJob "+jobname+" not found, (b"+str(ik)+"p"+str(ip)+"n"+str(-1*Ntor[intor])+"). Continue to try to open eigenvalues file.")
					if os.path.isfile(save_dir+"/"+case_dir+"/MISHKA_output/Eigenvalues.txt"):
						print("Eigenvalues.txt found and opened.\n")
					time.sleep(1)
					break
				else:
					sys.stdout.write("\r")
					sys.stdout.write("Job "+jobname+" not yet completed. Wait for {:2d} seconds.".format(holdon))
					sys.stdout.flush()
					time.sleep(1)	

# Obtain final eigenvalue from mishka output
			case_dir = "case_bped"+str(ik)+"_pcore"+str(ip)+"_Ntor"+str(Ntor[intor])
			os.chdir(save_dir+"/"+case_dir+"/MISHKA_output")
			if os.path.isfile("Eigenvalues.txt"):
				#print("Eigenvalues.txt exists!")
				evfile = open('Eigenvalues.txt','r')
				#print("Eigenvalues.txt nicely opened")
			else:
				datafile.write("{0:>5}\t{1:>6}\t{2:>9}\t{3:>6}\t{4:>5}\t\t{5:>32}\n".format(ik,ip,bped,pcore,Ntor[intor]," Eigenvalues.txt not found"))
				print(" Eigenvalues.txt not found")
				continue
			evlist = evfile.readlines()
			evstr = evlist[1]
			relchange = evlist[0].strip()

			evstr = evstr.replace("(","")
			evstr = evstr.replace(")","")
			evstr = evstr.strip()
	
			evvalues = evstr.split(',')

			ev0_check = float(evvalues[0])
			ev1_check = float(evvalues[1])
			
			if (ev0 is None) or (ev0_check > ev0_old):
				ev0_old = ev0_check
				ev0 = ev0_check
				ev1 = ev1_check
				Ntor_save = Ntor[intor]
				relCh_save = relchange
				betaN = find_normbeta(save_dir+"/"+case_dir+"/HELENA_output/fort.20")
				
			# Write output to datafile	
			os.chdir(save_dir)
			datafile.write(datalo.format(ik,ip,bped,pcore,betaN,Ntor[intor],float(ev0_check),float(ev1_check),float(relchange))) 
			#datafile.write("\n")
				
# If stability variable didn't change something has probably gone wrong (HELENA not converging, often Starting Value not found error: try changing B)
			if ev0_check == old_stab:
				print(case_dir+" : stab_crit did NOT change!")
				datafile.write("\t NO CHANGE \n")
				continue
			else:
				old_stab = ev0_check
				datafile.write("\n")
			#Remove large postscript files	
			os.system("rm "+save_dir+"/"+case_dir+"/HELENA_output/PCUBEQ")
			os.system("rm "+save_dir+"/"+case_dir+"/MISHKA_output/CASPLOT")
			
			if Ntor[intor] == 2:
				# Make plot of real-world Helena data 
				if numplot == 1:
					plot_realworld(job_dir,r,jpar,p0,jboot,bped,pcore,nbped,npcore,ip,ik)
			else: #Make a link to HELENA output in case for n=2, delete other HELENA output ( to save disk space)
				# os.chdir(save_dir+"/"+case_dir)
				# print "ln -s "+save_dir+"/case_bped"+str(ik)+"_pcore"+str(ip)+"_Ntor2/HELENA_ouput/ "+save_dir+"/"+case_dir+"/HELENA_output_link"
				# os.chdir(save_dir)
				# os.system("ln -s "+save_dir+"/case_bped"+str(ik)+"_pcore"+str(ip)+"_Ntor2/HELENA_ouput/ "+save_dir+"/"+case_dir+"/HELENA_output_link")		
				os.system("rm -r "+save_dir+"/"+case_dir+"/HELENA_output")
			evfile.close()
# Stability defined by eigenvalue (positive or negative)
		os.chdir(save_dir)
		if ev0 == None:
			datafile.write("{0:>5}\t{1:>6}\t{2:>9}\t{3:>6}\t     \t\t{4:>32}\n".format(ik,ip,bped,pcore," No eigenvalues found"))
			uf_data.write("{0:>5}\t{1:>6}\t{2:>9}\t{3:>6}\t     \t\t{4:>32}\n".format(ik,ip,bped,pcore," No eigenvalues found"))
			continue
		else:			
			stab_crit = ev0
			uf_data.write(datalo.format(ik,ip,bped,pcore,betaN,Ntor_save,float(ev0),float(ev1),float(relCh_save))) 
			if stab_crit == old_stab2:
				uf_data.write("\t NO CHANGE \n")
				continue
			elif (stab_crit > 0) and (abs(float(relCh_save)) > 1e-6):
				if stab_crit < 1e-4:
					old_stab2 = stab_crit
					print("Very small growth rate")
					uf_data.write("\t SMALL VALUE TAKEN AS STABLE \n")
					stab_crit = -1e-4
				else:	
					old_stab2 = stab_crit
					print("Not converged")
					uf_data.write("\t DID NOT CONVERGE \n")
					continue
			else:
				old_stab2 = stab_crit
				uf_data.write("\n")	
			
# Save values into arrays (add if stable, or loop over pcore completed)
			x.append(pcore) 
			y.append(bped)
			z.append(stab_crit)
			bN_plot.append(betaN)
			Ntor_plot.append(Ntor_save)
			
		print("stab_crit:" + str(stab_crit))# + "\n")
	
# Copy folders to Archive	
	
	if arch == 1:
		copytime_start = time.time()
		save_folder = "hel_"+helrun+"_"+helfile+"_mis_"+misfile+"__"+dt				
		cpstr = 'cp -r '+save_dir+'/case* /home/ITER/oosterw/Archive/Results/'+save_folder+'/'
		rmstr = 'rm -r '+save_dir+'/case* '
		os.system(cpstr)
		os.system(rmstr)
		print("Folders are copied to ~/Archive")
		copytime = copytime + (time.time()-copytime_start)
	
	#rsync -zvhr --progress --checksum --links iter_6:[directory WITHOUT end slash] [directory to where to copy it]
	#os.system("rsync -zvhr --progress --checksum --links iter_6:"+save_dir+"/"+case_dir+" /home/ITER/oosterw/Archive/Results/"+"hel_"+helrun+"_"+helfile+"_mis_"+misfile+"__"+dt)
	
# Check stability: stable for ev<0, unstable ev>0
# #	 Break out of the loop when stable
		# if (stab_crit <= 0): 
			# os.chdir(save_dir+"/"+case_dir+"/HELENA_output")
			# Bpt = find_Betapedtot("fort.20")
			# print("Poloidal beta, total and pedestal only: "+str(Bpt)+" , "+str(bped))
			# #bped_start = bped
			# break
			# #continue
			# #print("Stable solution")
	os.chdir(save_dir)

# # Unstable for all bped-pcore combinations
	# if ((ik==nbped-1) and (stab_crit>0)):
		# print("\n No convergence to stability reached \n")
		# datafile.write("No stability reached for pcore = "+str(pcore)+"\n")
		# continue	

os.chdir(save_dir)
# Write x,y,z to datafile
xyzfile = open("plotdata.txt",'w') 
layout = "{0:>4}\t{1:>4}\t{2:>8}\t{3:>8}\t{4:>3}\n"
for l in range(len(x)):
	xyzfile.write(layout.format(x[l],y[l],z[l],bN_plot[l],Ntor_plot[l]))

# Make plots

if numplot == 1:
	plt.figure(1)
	plt.legend(loc=2,bbox_to_anchor=(1.3,1.0))
	jparplot.savefig("jpar.eps")
	plt.figure(2)
	plt.legend(loc=2)
	presplot.savefig("pressure.eps")
	plt.figure(3)
	plt.legend(loc=10)
	jbootplot.savefig("jbootstrap.eps")

	# pdf.savefig(jparplot)
	# pdf.savefig(presplot)
	# pdf.savefig(jbootplot)


fig = plt.figure()
#plot = plt.plot(x,y)
#plot = plt.contourf(X, Y, Z, 10, cmap=plt.cm.bone, origin=origin)
plot = plt.scatter(x,y,c=z,marker='o')
#for i in range(len(z)):
#    plt.annotate(s='{:.1e}'.format(z[i]),xy=(x[i],y[i]))
plt.colorbar()
plt.xlabel("Core pressure parameter")
plt.ylabel("Beta poloidal pedestal")
plt.title(dt)
plt.grid(b=True, which="major", linestyle="-")
plt.grid(b=True, which="minor", linestyle="--")
fig.savefig("StabilityGraph.eps")

#plt.figure()
#CP = plt.contourf(x,y,z)
#plt.show
#fig.savefig(ContourStability.ps)

#Close all opened files
xyzfile.close()
datafile.close()
uf_data.close()

# Copy folders to Archive	
if arch == 1:
	save_folder = "hel_"+helrun+"_"+helfile+"_mis_"+misfile+"__"+dt				
	cpstr = 'cp -r '+save_dir+'/* /home/ITER/oosterw/Archive/Results/'+save_folder+'/'
	os.system(cpstr)
	print("Final files are copied to ~/Archive , not removed")

os.chdir(job_dir)

#os.system("mv -v stdout.txt "+save_dir+"/")

print("--- Duration: %s seconds --- Of which %s seconds for copying" %(time.time() - start_time) , copytime )
