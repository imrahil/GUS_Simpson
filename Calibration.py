# THIS IS A ConceptFORGE PRODUCT.
# GPL LICENSE

# Calibration.py
# Author: Niggle (aka Nigel Haslock)
# Creation Date: May 30 2014

# Usage:	Calibration [-i {points filename}] [-o {calibration filename}]
#				-i specifies the name of a file containing points data. If not specified, the program
#					will look for Points.GUS in the current directory. If that fails, the program
#					will prompt for a file.
#				-o specifies the name of a file to hold the calibration results. If not specified, the program
#					will use Calibration.GUS in the current directory.
#
# File formats
#	Point data files, such as Points.GUS, are plain text.
#
#	Each line may be blank, a comment or a set of three values.
#
#	Comment lines must have '#' or ';' as the first non-whitespace character.
#
#	Values lines must have three numeric values separated by ',' or whitespace.
#		e.g. 137.1 137.2, 137.3
#
#	There must be at least six value lines in the file.
#
#	Calibration result files consist of two value lines and some comment lines
#	The first value line reports the calculated shoulder height for each arm
#	The second value line reports the calculated bolt to bolt distance for 
#		each arm when the endstop is triggered.

# All lengths are in mm
POINTS=[]
SIZE=250.0			# Distance between shoulder pivot bolts (i.e. the vertical bolts through the bed
BED_Z=79.5			# Approximate vertical distance between arm pivot bolts when the extruder is touching the bed
MAX_ARM_LENGTH=282	# Approximate straight line distance between arm pivot bolts when the endstop is triggered
HEIGHT_COMPENSATION=0.1

from scipy.optimize import leastsq
import numpy.linalg
import math, random, copy, sys
import os
import stat
import time

def num(s):
#    try:
#        return int(s)
#    except ValueError:
        return float(s)
		
# Parse commandline 
pointfile="Points.GUS"
calibfile="Calibration.GUS"

for i in range(len(sys.argv)):
	if sys.argv[i]=="-i" and len(sys.argv)>i+1:
		pointfile=sys.argv[i+1]
	if sys.argv[i]=="-o" and len(sys.argv)>i+1:
		calibfile=sys.argv[i+1]

# Check if the calibration data file exists
cfileinfo = []
ifileinfo = []
try:
	fileStats = os.stat(calibfile)
	cfileinfo = time.ctime(fileStats[stat.ST_MTIME])
except:
	cfileinfo = 0
try:
	fileStats = os.stat(pointfile)
	ifileinfo = time.ctime(fileStats[stat.ST_MTIME])
except:
	print "Failed to find file " + pointfile
	quit()
	
if cfileinfo > ifileinfo:
	print "Calibration file is newer than point data file. Done."
	quit()
	

try:
	inpoints=open(pointfile)
except:
	inpoints=file(raw_input("Points File: "))
try:
	outcalib=open(calibfile,"w")
except:
	inpoints=file(raw_input("Calibration File: "),"w")
	
# Start point for calibration analysis
DEFAULT_VALUES=[BED_Z,BED_Z,BED_Z,MAX_ARM_LENGTH,MAX_ARM_LENGTH,MAX_ARM_LENGTH]
# Declaration of where the results are stored
SHOULDER_Z1,SHOULDER_Z2,SHOULDER_Z3,MAX_LENGTH_1,MAX_LENGTH_2,MAX_LENGTH_3=DEFAULT_VALUES

outcalib.write("# Calibration data derived from the following points data \n")

pointcount=0
for rawline in inpoints:
	outcalib.write("# " + rawline)
	line=rawline.upper()
	line=line.strip()
	if line == "":
		continue
	if line[0] == "#" or line[0] == ";":
		continue
	processed_line=""
	flagSpace = int(0)
	for letter in line:
		if letter.isdigit() or letter == ".":
			processed_line+=letter
			flagSpace = 0
		elif flagSpace == 0:
			flagSpace = 1
			processed_line += " "
		else:
			pass
	coords = [float(i) for i in processed_line.split()]
	if len(coords) != 3:
		print "Invalid line " + line + " in Points file"
	else:
		# Is there a way to have POINTS be an array of objects rather than an array of floats?
		coords[0] += HEIGHT_COMPENSATION;
		coords[1] += HEIGHT_COMPENSATION;
		coords[2] += HEIGHT_COMPENSATION;
		POINTS+=coords
		pointcount+=1
inpoints.close()

if pointcount<6:
	print "Insufficient valid points in Points file"
	outcalib.write("#\n#\n# Calibration aborted\n")
	quit()
	
# Convert arms lengths to cartesian coordinates using TRILATERATION
def getxyz(r1,r2,r3):
    d=SIZE*1.0
    i=SIZE/2.0
    j=SIZE*math.sqrt(3)/2.0
    x=(r1*r1-r2*r2+d*d)/(2*d)
    y=(r1*r1-r3*r3-x*x+(x-i)**2+j*j)/(2*j)
    zsq=(r1*r1-x*x-y*y)
    if zsq<0:
        # TRILATERATION failed so we force to 0
        print "Trilateration failure - calibration suspect"
        outcalib.write("# Trilateration failure - calibration suspect\n")
        zsq=0
    z=math.sqrt(zsq)
    return x,y,z

#GET VALUES FOR EACH POINT TO SEE HOW CLOSE TO THE PLANE IT IS
def equations(p):
    SHOULDER_Z1,SHOULDER_Z2,SHOULDER_Z3,MAX_LENGTH_1,MAX_LENGTH_2,MAX_LENGTH_3=p
    m=[]
    for i in range(len(POINTS)/3):
        R1,R2,R3=MAX_LENGTH_1-POINTS[i*3],MAX_LENGTH_2-POINTS[i*3+1],MAX_LENGTH_3-POINTS[i*3+2]
        X,Y,Z=getxyz(R1,R2,R3)
        d=SIZE*1.0
        i=SIZE/2.0
        j=SIZE*math.sqrt(3)/2.0
        q=[[0,0,SHOULDER_Z1,1],
           [d,0,SHOULDER_Z2,1],
           [i,j,SHOULDER_Z3,1],
           [X,Y,Z,1]]
        det=numpy.linalg.det(q)
        m.append(det**2)
    return m

# Derive calibration values from point data
SHOULDER_Z1,SHOULDER_Z2,SHOULDER_Z3,MAX_LENGTH_1,MAX_LENGTH_2,MAX_LENGTH_3=leastsq(equations, DEFAULT_VALUES)[0]
text="\n\n" + str(SHOULDER_Z1) + ", " + str(SHOULDER_Z2) + ", " + str(SHOULDER_Z3) + " \n"
print text
outcalib.write(text)
text=str(MAX_LENGTH_1) + ", " + str(MAX_LENGTH_2) + ", " + str(MAX_LENGTH_3) + "\n"
print text
outcalib.write(text)
outcalib.close()

print "done"

