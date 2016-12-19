from glob import glob
import os

xmlfilelist = glob('*.xml')

for xml in xmlfilelist:
	#loadCommand = 'module load BEAST'
	commandText = 'sbatch --exclude=/home/sbell23/gizmod.txt --time=7-0 --wrap=\"java -jar /home/sbell23/build/beast.jar %s"'%xml
	os.system(commandText)
