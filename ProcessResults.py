#===============================================================================
# Imports
#===============================================================================

import commands, math, numpy, scipy, os, pdb, shutil, sys, pickle
from scipy import stats

#===============================================================================
# Parameters
#===============================================================================

# File path to directory
base_path = '/User/username/directory/droplet-analysis/'
# Note: Must match sigma_fluid from the 'SetupSystem.py' file
sigma = 0.3045
# Note: Must match sigma_wall from the 'SetupSystem.py' file
sigma_wall = 0.2705
# Number of files that will be grouped into a single section for processing
fileSectionSize = 2000
# Side length of a density sampling box in terms of sigma
# Note: This must be set to be optimal so a section of files, not the total
# number of files being analyzed. If the box side length is too large or small
# for a section of files, the results will not be correct
# This is also used as the bin width for centering each file by density peak
density_length = 0.20
# The fraction of the total system fluid height (spacing between innermost wall
# roughness) that will be fit with a linear least squares regression
# Note: The first 2 sigmas of height above the wall surfaces will not be
# considered due to unpredictable fluid-surface interactions. However, the
# correct fraction of total height on both side of the fluid will be considered
dsigma = 0.10
dsigmaCircle = 1.0
# Time step used to calculate the velocity of the center of mass
# Note: Must match dt from the 'SetupSystem.py' file
# Units are in picoseconds
timestep = 0.01
# Frequency of writing system state to .xtc file. Must match .mdp settings. 
sampleFrequency = 10
# Full output contaning all distribution profiles
fullResults = False

#===============================================================================
# Subroutines
#===============================================================================


def extract_wall(atomsX, atomsY, atomsZ, numAtoms):
    """
    Given two arrays containing the Y and Z coordinates of the wall atoms for a
    single system snapshot and the number of atoms in the system, this will
    obtain all required wall information in addition to the number of fluid and
    wall atoms in the system
    """

    # Determines the minimum and maximum Y and Z values of the walls
    # Wall sigma additional length must be considered in Y for any
    # fluid greater than the last wall center but within the system box
    maxWallY = (max(atomsY) + sigma_wall)
    minWallY = min(atomsY)
    maxWallZ = max(atomsZ)
    minWallZ = min(atomsZ)
    maxWallX = max(atomsX)
    minWallX = min(atomsX)
    # Determines wall thickness based on the thickness of the bottom wall
    # Middle Z value between the walls is determined
    maxThicknessZ = ((maxWallZ + minWallZ)/2)
    
    # Stores all Z values for wall atoms below the midpoint between the walls
    wallThicknessBottom = list()
    # For every wall atom Z coordinate in the system
    for atom in atomsZ:
        # If the Z value of the atom is below the wall mid plane Z value
        if (atom < maxThicknessZ):
            # Adds the Z value to the appropriate list
            wallThicknessBottom.append(atom)
    # Thickness is temporarily equal to the maximum bottom wall Z value
    thicknessBottom = max(wallThicknessBottom)
    # Thickness of the wall is equal to thickness minus minWallZ
    thicknessBottom -= minWallZ
    
    # Stores all Z values for wall atoms above the midpoint between the walls
    wallThicknessTop = list()
    # For every wall atom Z coordinate in the system
    for atom in atomsZ:
        # If the Z value of the atom is above the wall mid plane Z value
        if (atom > maxThicknessZ):
            # Adds the Z value to the appropriate list
            wallThicknessTop.append(atom)
    # Thickness is temporarily equal to the minimum top wall Z value
    thicknessTop = min(wallThicknessTop)
    # Thickness of the wall is equal to maxWallZ minus thickness
    thicknessTop = maxWallZ - thicknessTop

    # Assigns thickness to be equal to the thickness of the thickest wall
    # If the top wall is the thickest wall
    if (thicknessTop > thicknessBottom):
        # The wall thickness equals the thickness of the top wall
        thickness = thicknessTop
    # If the bottom wall is thicker than or equal to the thickness of
    # the top wall
    else:
        # The wall thickness equals the thickness of the bottom wall
        thickness = thicknessBottom

    # Target center of the fluid is wall center in the Y dimension
    targetCenter = ((maxWallY + minWallY)/2)
    # Determines the number of wall atoms in the system
    numWal = len(atomsY)
    # Determines the number of fluid atoms in the system
    numFld = (numAtoms - numWal)
    # Determines the boxes between the walls in the system x dimension
    boxesX = int((maxWallX - minWallX)/(density_length*sigma))

    # Returns values representing important wall features
    return targetCenter, minWallY, maxWallY, minWallZ, maxWallZ, thickness, numWal, numFld, boxesX


def extract_fluid(fileAtoms, density_length, targetCenter, minWallY, maxWallY, minWallZ, maxWallZ):
    """
    Given the Y and Z atom positions for all of the atoms in one file in
    addition to various descriptive values, this centers the system based on
    density peak then center of mass in the Y direction before adding the atoms
    to a two dimensional matrix containing integers to represent the number of
    atoms contained within each density sampling box
    """

    # Matrix of density sampling boxes used before the reference
    global densityBoxes, positions, centeringDistance
    
    # Determines the fluid center of mass system length dimension position
    velocityTrack = numpy.average(fileAtoms[:,0])
    # Adds the center of mass position to a global list so that the center of
    # mass velocity can be determined
    positions.append(velocityTrack)
    
    # Initially centers by density maximum in Y to prevent the fluid from being
    # split apart by the periodic boundary conditions
    # Accumulator to store the Y values of atoms to find the average
    # One dimensional array stores the atoms in each Y cross section
    numcenters = numpy.ceil((maxWallY - minWallY)/(sigma*density_length))
    centerSections = numpy.zeros(numcenters, int)
    # Increment the atoms Y cross section density bin of each atom in the file
    bincenters = numpy.floor((fileAtoms[:,0] - minWallY)/(sigma*density_length))
    centerSections = numpy.histogram(bincenters, bins=numcenters, range=[0,numcenters])[0]
    
    # Determine the index of the peak density cross section and converts it back
    # to the position of the bin center
    centerDensityY = float(((numpy.argmax(centerSections))*(sigma*density_length)) + ((density_length*sigma)/2) + minWallY)
    # Determines the adjustment factor that must be added to all Y values
    adjustment = (targetCenter - centerDensityY)
    # Determines the space between the walls in the Y dimension
    betweenwall = (maxWallY - minWallY)
    # Adds the adjustment factor to the Y coordinate or every atom
    fileAtoms[:,0] += adjustment
    # Determines the fluid atoms below the lower wall Y limit
    toosmall = (fileAtoms[:,0] < minWallY)
    # Adds the spacing between the walls to the atoms below this limit
    fileAtoms[toosmall,0] += betweenwall
    # Determines the fluid atoms above the upper wall Y limit
    toobig = (fileAtoms[:,0] > maxWallY)
    # Subtracts the spacing between the walls from the atoms above this limit
    fileAtoms[toobig,0] -= betweenwall

    # Centers the fluid atoms by center of mass
    # Determines the average fluid atom Y location
    avY = numpy.average(fileAtoms[:,0])
    # Determines the adjustment factor that must be added to Y to center it 
    adjustment2 = (targetCenter - avY)
    # Adds the adjustment factor to the Y coordinate or every atom
    fileAtoms[:,0] += adjustment2
    # Determines the fluid atoms below the lower wall Y limit
    toosmall = (fileAtoms[:,0] < minWallY)
    # Adds the spacing between the walls to the atoms below this limit
    fileAtoms[toosmall,0] += betweenwall
    # Determines the fluid atoms above the upper wall Y limit
    toobig = (fileAtoms[:,0] > maxWallY)
    # Subtracts the spacing between the walls from the atoms above this limit
    fileAtoms[toobig,0] -= betweenwall

    # Determines the X and Z density sampling box index for all fluid atoms
    boxY = numpy.floor((fileAtoms[:,0] - minWallY)/(sigma*density_length)).astype(int)
    boxZ = numpy.floor((fileAtoms[:,1] - minWallZ)/(sigma*density_length)).astype(int)
    # Determines the dimensions of the density sampling box matrix
    boxSize = numpy.shape(densityBoxes)
    # Forms a density sampling box matrix for all fluid atoms of a file
    density = numpy.histogram2d(boxY, boxZ, bins=boxSize, range=[[0,boxSize[0]], [0,boxSize[1]]])[0]
    # Adds the matrix of file fluid atoms to the matrix of section fluid atoms
    densityBoxes += density
    
    # Void return as all atoms have been added to the respective element of
    # a two dimensional matrix used to represent density sampling boxes
    return

def process(densBoxes, density_length, dsigma, dsigmaCircle, thickness, snapshots, boxesX):
    """
    The density distribution is determined by dividing the number of atoms
    within each sampling box by the number of system snapshots contained in
    the subsample and the volume of the density sampling box in order to obtain
    a distribution that does not depend on the subsample size or density box
    side length. The cutoff density is determined based on half the peak
    sampling box density. Boxes with less than the cutoff density are set to
    zero density and ignored for further calculations. Of the boxes with
    sufficient density, the outermost points in the Y dimension are determined
    at each Z value and, if those points fall within a specified range of fluid
    height to be considered, they are added to a list of all points within that
    appropriate quadrant. These points will later be fit with a linear least
    squares regression in order to determine the contact angles. 
    """

    # Makes the density sampling box distribution independent of the number of
    # system snapshots within the subsample and volume of the density box
    denseBoxes = (densBoxes/(snapshots*boxesX*((sigma*density_length)**3)))
    # Obtains the two dimensional initial distribution
    collapsedBoxesInitial = denseBoxes.copy()
    # Obtains one dimensional initial distributions about each axis
    vertDistributionInitial = numpy.sum(collapsedBoxesInitial, axis=0)
    horizDistributionInitial = numpy.sum(collapsedBoxesInitial, axis=1)
    
    # Determines the index dimensions of the collapsed density matrix
    shape = numpy.shape(denseBoxes)
    ys = shape[0]
    zs = shape[1]
    # Determines the number of elements in the collapsed density matrix
    numElements = float(ys*zs)
    
    # Calculates the median modified density of a sampling box
    densityCutoff = numpy.mean(denseBoxes)
    # Values to be used in determining the modified density half way between
    # that of bulk fluid and bulk vapor
    fluidDensity = 0
    fluidBoxes = 0
    vaporDensity = 0
    vaporBoxes = 0
    # For every y index in the sampling box matrix
    for yDim in range(0,ys):
        # For every z value at that y value
        for zDim in range(0,zs):
            # If the sampling box is in a fluid region
            if(denseBoxes[yDim,zDim] > densityCutoff):
                # Increment the number of fluid sampling boxes
                fluidBoxes = fluidBoxes + 1
                # Add the density of the sampling box to the fluid accumulator
                fluidDensity = fluidDensity + denseBoxes[yDim,zDim]
            # Otherwise the sampling box is in a vapor region
            else:
                # Increment the number of vapor sampling boxes
                vaporBoxes = vaporBoxes + 1
                # Add the density of the sampling box to the vapor accumulator
                vaporDensity = vaporDensity + denseBoxes[yDim,zDim]
    # Set the threshold density equal to density half way between bulk fluid
    # and bulk vapor
    densityCutoff = ((vaporDensity/vaporBoxes)+(fluidDensity/fluidBoxes))/2

    
    # Set all density sampling boxes with below the cutoff density atoms to zero
    denseBoxes[denseBoxes < densityCutoff] = 0

    # Obtains one dimensional filtered distributions about each axis
    vertDistribution = numpy.sum(denseBoxes, axis=0)
    horizDistribution = numpy.sum(denseBoxes, axis=1)

    # Creates a list to store all outermost points that will be considered
    points = list()
    # For every Z value in the two dimensional system
    for z in range(zs):
        # Determine the Y box indices of all nonzero boxes at that Z value
        nonzeros = numpy.nonzero(denseBoxes[:,z])[0]
        # If there is at least one nonzero density sampling box at that Z value
        if (numpy.size(nonzeros) > 0):
            # Determines the maximum and minimum Y box index at that Z value
            yMin = numpy.min(nonzeros)
            yMax = numpy.max(nonzeros)
        # If there are no nonzero box values at that Z value
        else:
            # Sets the maximum and minimum box indices to infinite values
            yMax = float('-inf')
            yMin = float('inf')
        # If the current Z index contains at least one box to be considered
        if ((yMax != float('-inf')) and (yMin != float('inf'))):
            # Converts the current Z index to the position of the box center
            Z = float((z*sigma*density_length) + (0.5*sigma*density_length) + minWallZ)
            # Converts the outermost Y indices at that Z value back to position
            yl = float((yMin*sigma*density_length) + (0.5*sigma*density_length) + minWallY)
            yu = float((yMax*sigma*density_length) + (0.5*sigma*density_length) + minWallY)
            # Adds the points to a list of all outside fluid points
            points.append([yl, Z])
            points.append([yu, Z])
                
    # Lists to store the four quadrants of final points
    upper_right = list()
    upper_left = list()
    lower_right = list()
    lower_left = list()
    # Lists for all circular regression points
    circle_list = list()
    circlePoints_advancing = list()
    circlePoints_receding = list()
    
    # Adjusts dsigma to obtain an integer number of boxes of height to analyze
    dsigmaBoxes = int((dsigma*(maxWallZ - minWallZ - (2*thickness)))/(density_length*sigma))
    # Determines lower and upper cutoff Z values based on the number of boxes
    lowerCutoff = (minWallZ + thickness + (sigma*dsigmaBoxes*density_length) + (2*sigma))
    upperCutoff = (maxWallZ - thickness - (sigma*dsigmaBoxes*density_length) - (2*sigma))
    # For all boundary sampling box positions with sufficient fluid density
    for point in points:
        # Note the first two sigmas of fluid height is not considered due
        # to unpredictable behavior where the fluid contacts the surface
        # If the point is above the upper cutoff Z value and greater than the
        # center Y value
        if ((float(point[1]) > upperCutoff) and (float(point[1]) < (maxWallZ - thickness - (2*sigma))) and (float(point[0]) > targetCenter)):
            # Add the point to the upper right list
            upper_right.append(point)
        # If the point is above the upper cutoff Z value and less than the
        # center Y value
        elif ((float(point[1]) > upperCutoff) and (float(point[1]) < (maxWallZ - thickness - (2*sigma))) and (float(point[0]) < targetCenter)):
            # Add the point to the upper left list
            upper_left.append(point)
        # If the point is below the lower cutoff Z value and greater than the
        # center Y value
        elif ((float(point[1]) < lowerCutoff) and (float(point[1]) > (minWallZ + thickness + (2*sigma))) and (float(point[0]) > targetCenter)):
            # Add the point to the lower right list
            lower_right.append(point)
        # If the point is below the lower cutoff Z value and less than the
        # center Y value
        elif ((float(point[1]) < lowerCutoff) and (float(point[1]) > (minWallZ + thickness + (2*sigma))) and (float(point[0]) < targetCenter)):
            # Add the point to the lower left list
            lower_left.append(point)

    # Center of circular regression height range
    circleRegressionCenter = (minWallZ + ((maxWallZ - minWallZ)/2))
    # Upper and lower height range cutoffs for circular regression analysis
    circleLowerCutoff = (circleRegressionCenter - (0.5*dsigmaCircle*((maxWallZ - minWallZ) - (2*thickness))))
    circleUpperCutoff = (circleRegressionCenter + (0.5*dsigmaCircle*((maxWallZ - minWallZ) - (2*thickness))))
    # For all boundary sampling box positions with sufficient fluid density
    for point in points:
        # Note the first two sigmas of fluid height is not considered due
        # to unpredictable behavior where the fluid contacts the surface
        # If the point is within the circular regression height range and not
        # within two sigmas of either wall and is an advancing point
        if ((float(point[1]) < circleUpperCutoff) and (float(point[1]) > circleLowerCutoff) and (float(point[1]) < (maxWallZ - thickness - (2*sigma))) and (float(point[1]) > (minWallZ + thickness + (2*sigma))) and (float(point[0]) > targetCenter)):
            # Add the point to the list of advancing points
            circlePoints_advancing.append(point)
        # If the point is within the circular regression height range and not
        # within two sigmas of either wall and is a receding point
        if ((float(point[1]) < circleUpperCutoff) and (float(point[1]) > circleLowerCutoff) and (float(point[1]) < (maxWallZ - thickness - (2*sigma))) and (float(point[1]) > (minWallZ + thickness + (2*sigma))) and (float(point[0]) < targetCenter)):
            # Add the point to the list of advancing points
            circlePoints_receding.append(point)
    
    # Returns four lists of points that will each be fit with a linear
    # regression in order to determine contact angles for the system
    # in addition to distributions of fluid density
    return circlePoints_advancing, circlePoints_receding, upper_right, lower_right, upper_left, lower_left, densityCutoff, dsigmaBoxes, collapsedBoxesInitial, denseBoxes, horizDistributionInitial, horizDistribution, vertDistributionInitial, vertDistribution

#===============================================================================
# Main
#===============================================================================

# Adds the 'simulation' directory to the base_path to form the work_path
work_path = base_path + 'simulation/'
# Changes directories to the work_path
os.chdir(work_path)
# Stores center of mass positions
positions = list()
# Opens the 'analyze.gro' file
inFile = open('analyze.gro', 'r+')
# Reads past the file header
discard = inFile.readline()
# Obtains the number of atoms in the system from the second line
numAtoms = int(str(inFile.readline()).strip())
# Determines the total number of lines in a single file
linesFile = (numAtoms + 3)
# Stores the distance the fluid center of each file is shifted while centering it
centeringDistance = float(0)
# Generates lists to store sample section distributions
saveSampleSections = list()

# Lists to store the wall atoms Y and z coordinates respectively
wallAtomsX = list()
wallAtomsY = list()
wallAtomsZ = list()
# While the current atom location is less than the total number of atoms
for number in range(linesFile - 2):
    # Obtain the current line based on the atom index
    line = inFile.readline()
    # If the line in the first file represents a wall atom
    if ('WAL' in line):
        # Add the wall atoms coordinates to the appropriate list
        wallAtomsX.append(float(line[20:28]))
        wallAtomsY.append(float(line[28:36]))
        wallAtomsZ.append(float(line[36:44]))
# Analyze the wall atoms to obtain various descriptive system values
targetCenter, minWallY, maxWallY, minWallZ, maxWallZ, thickness, numWal, numFld, boxesX = extract_wall(wallAtomsX, wallAtomsY, wallAtomsZ, numAtoms)

# Generates a matrix of zeros to represent the section sampling box densities
densityBoxes = numpy.zeros([numpy.ceil((maxWallY - minWallY)/(sigma*density_length)), numpy.ceil((maxWallZ - minWallZ)/(sigma*density_length))], int)
# Generates a matrix of zeros to represent the total sampling box densities
totalDensityBoxes = numpy.zeros([numpy.ceil((maxWallY - minWallY)/(sigma*density_length)), numpy.ceil((maxWallZ - minWallZ)/(sigma*density_length))], int)
# Establishes a numpy matrix for all fluid atom positions of a file
fileAtoms = numpy.empty([numFld,2], dtype='float')

# Returns to the beginning of the file so it can be read again
inFile.seek(0, 0)
# For all of the lines in the first file before the first fluid atom
for number in range(numWal + 1):
    # Read the line but do nothing with it
    line = inFile.readline()

# Accumulator for the current file section
filesSection = 0
# Accumulator for the current file number
fileNum = 0
# Accumulator for the current line number of the current file
lineNum = 0
# While there is a line remaining in the system
while (line):
    # Read the current line of the file, which will always be a fluid atom
    line = inFile.readline()
    # Obtain the Y and Z coordinates of the fluid atom
    fileAtoms[lineNum,0] = numpy.float(line[28:36])
    fileAtoms[lineNum,1] = numpy.float(line[36:44])
    # Increment the line number as one more line is analyzed
    lineNum += 1
    # If the line number reaches the total number of fluid atoms
    if (lineNum > (numFld - 1)):
        # Increment file number as another file has been analyzed
        fileNum += 1
        # Reset line number as a new file is reached
        lineNum = 0
        # Increment the number of files analyzed in the current section
        filesSection += 1
        # Extract the fluid atoms when a complete snapshot has been read
        extract_fluid(fileAtoms, density_length, targetCenter, minWallY, maxWallY, minWallZ, maxWallZ)
        # Reset the array holding the fluid atoms of an individual file
        fileAtoms = numpy.empty([numFld, 2], float)
        # For the total number of non fluid lines after the fluid lines of this
        # file end and before the fluid lines of the next file begin
        for number in range(numWal + 3):
            # Read the line but do nothing with it
            line = inFile.readline()
        # If the number of files to be analyzed per section is reached
        if (filesSection > (fileSectionSize - 1)):
            # Process the accumulated density of the section to obtain points
            sectionCircle_advancing, sectionCircle_receding, upper_right, lower_right, upper_left, lower_left, densityThreshold0, numBoxesHigh0, distribution00, distribution0, horizDistribution00, horizDistribution0, vertDistribution00, vertDistribution0 = process(densityBoxes, density_length, dsigma, dsigmaCircle, thickness, fileSectionSize, boxesX)
            # Adds the section density boxes to the cumulative density boxes
            totalDensityBoxes += densityBoxes
            # Adds the sample section density distribution to the a list so it can be saved
            saveSampleSections.append(densityBoxes)
            # Resets the matrix used to represent the section density boxes
            densityBoxes = numpy.zeros([numpy.ceil((maxWallY - minWallY)/(sigma*density_length)), numpy.ceil((maxWallZ - minWallZ)/(sigma*density_length))], int)
            # Reset the number of files analyzed in the current group
            filesSection = 0

# Closes the file
inFile.close()

# Process the total density of all sections to obtain points cumulative points
totalCircle_advancing, totalCircle_receding, upper_right_total, lower_right_total, upper_left_total, lower_left_total, densityThreshold, numBoxesHigh, distribution0, distribution, horizDistribution0, horizDistribution, vertDistribution0, vertDistribution = process(totalDensityBoxes, density_length, dsigma, dsigmaCircle, thickness, fileNum, boxesX)

# Generate a file for processed data to be written to
d = open('Data.pkl','w+')

# Formats all required data so it can be saved
saveData = [saveSampleSections, totalDensityBoxes, targetCenter, minWallY, maxWallY, minWallZ, maxWallZ, thickness, boxesX, numWal, numFld, numAtoms, wallAtomsX, wallAtomsY, wallAtomsZ, positions, sigma, sigma_wall, fileSectionSize, density_length, dsigma, dsigmaCircle, timestep, sampleFrequency, fullResults, fileNum]
# Saves all required data so it can later be analyzed
pickle.dump(saveData,d)

# Closes the 'Results.txt' file
d.flush()
d.close()
