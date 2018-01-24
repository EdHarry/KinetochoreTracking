from subprocess import call
import sys

imarisPath = sys.argv[1]
id = sys.argv[2]

imarisPath = imarisPath.replace('~w',' ')

call(["/Users/Ed/Documents/MATLAB/maki/imarisICE/startImaris.sh", imarisPath, id])








