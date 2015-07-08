import os
import sys
import ConfigParser
import string
import time


def RunSuperMCWithArgs(args = '',instances = 8):
	"""Runs ./superMC with [instances] terminals at once with the correct arguments"""
	i = 1
	while i < instances:
		os.system('gnome-terminal -e \"./superMC.e ' + args+"\"")
		i+=1

	time.sleep(6)
	os.system('./superMC.e ' + args)

	UpdateParameters(args)

def UpdateParameters(args):
	
	fout = open('ModParameters.dat','w')
	fin = open('parameters.dat','r')

	separateArgs = args.split()
	argNames = [arg[0:arg.index('=')] for arg in separateArgs]

	for line in fin:
		writeThis = line
		words = line.split()
		for i in range(len(argNames)):
			if len(words) > 0:
				if words[0] == argNames[i]:
					writeThis = separateArgs[i] + '\n'

		fout.write(writeThis)

	fout.close()
	fin.close()


def CutAndMoveData(plotUtilityDirectory,destination,args,name,cutFlag = True):
	"""Cuts the data using the pythonPlotUtility. Requires the location be specified"""

	if name is not None:
		folderName = name
	else:
		folderName = string.replace(args,' ','_')

	destinationFolder = destination + '/' + folderName

	os.system('mkdir -p ' + destinationFolder)
	centralityCutParams = ['Npart','total_entropy']

	#Make the database
	os.system('python '+plotUtilityDirectory + '/EbeCollectorShell_minbisaEcc.py data')

	if cutFlag:
		#Cut the data
		for param in centralityCutParams:
			os.system('python ' + plotUtilityDirectory + '/minbias_CentralityCut.py data/minbiasEcc_sn.db -type ' + param)
			#os.system('mv eccnStatistics_ed_'+param+'.dat ' + destinationFolder)
			os.system('mv eccnStatistics_sd_'+param+'.dat ' + destinationFolder)
			os.system('mv centralityCut_'+param+'.dat ' + destinationFolder)

		os.system('python ' + plotUtilityDirectory + '/minbias_ImpactParameterCut.py data/minbiasEcc_sn.db')
		#Move the files over
		#os.system('mv eccnStatistics_ed_b.dat ' + destinationFolder)
		os.system('mv eccnStatistics_sd_b.dat ' + destinationFolder)
		os.system('mv centralityCut_b.dat ' + destinationFolder)

	os.system('mv data/minbiasEcc_sn.db ' + destinationFolder)
	os.system('mv ModParameters.dat ' + destinationFolder);

if __name__ == "__main__":

	Config = ConfigParser.ConfigParser()
	Config.read('BatchCollectionSettings.txt')

	cutFlag = True
	for arg in sys.argv:
		if arg == '-noCut':
			cutFlag = False
			sys.argv.remove(arg)

	locationOfPlotUtility = Config.get('FolderParams','LocationOfPlotUtility')
	destination = Config.get('FolderParams','Destination')


	superMCArgs = Config.items('SuperMCArgs')

	instances = Config.get('RunTimeParams','Instances')
	names = Config.items('FileNames')

	for i in range(len(superMCArgs)):
		RunSuperMCWithArgs(superMCArgs[i][1],int(instances))
		if i < len(names):
			CutAndMoveData(locationOfPlotUtility,destination,superMCArgs[i][1],names[i][1],cutFlag)
		else:
			CutAndMoveData(locationOfPlotUtility,destination,superMCArgs[i][1],None,cutFlag)
		os.system('rm -fr data/*')

