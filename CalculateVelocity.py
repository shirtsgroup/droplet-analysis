#===============================================================================
# Imports
#===============================================================================

import commands, math, numpy, scipy, os, pdb, shutil, sys, pickle
from scipy import stats
import matplotlib.pyplot as plt

#===============================================================================
# Parameters
#===============================================================================

# File path to directory containing the stored data
work_path = '/User/username/directory/superhydrophobic/'
# Starting index values for velocity calculations
startIndex = [1000, 10000, 20000]
# Time analyzed per segment beginning at each starting index value
boundLength = 5000
# Write results to a file or only display the results
# writeResults = True
writeResults = False

#===============================================================================
# Main
#===============================================================================

# Changes directories to the work_path
os.chdir(work_path)
# Opens the file containing the saved and processed data
pickleData = open('Data.pkl','r+')
# Loads processed data for analysis
saveSampleSections, totalDensityBoxes, targetCenter, minWallY, maxWallY, minWallZ, maxWallZ, thickness, boxesX, numWal, numFld, numAtoms, wallAtomsX, wallAtomsY, wallAtomsZ, positions, sigma, sigma_wall, fileSectionSize, density_length, dsigma, dsigmaCircle, timestep, sampleFrequency, fullResults, fileNum = pickle.load(pickleData)
# Closes the 'Data.pkl' file
pickleData.flush()
pickleData.close()

# List to hold all velocities calculated from endpoint values
checkVelocities = list()
# Calculated the velocity from the endpoints of each segment and averages them
# in order to provide a check that the reported velocity doesnt include non linear
# portions of the segment
# For all segments
for startingIndexCheck in startIndex:
    # Ending index for the segment
    endingIndex = startingIndexCheck + boundLength -1
    # Adds the velocity calculated from endpoint values to the list of values
    checkVelocities.append((positions[endingIndex]-positions[startingIndexCheck])/float(0.001*timestep*sampleFrequency*boundLength))
# Calculates check velocity
checkVelocity = numpy.average(checkVelocities)

# List to hold all velocities
velocities = list()
# Determines all velocity values
for positionIndex in range(1, (len(positions)-1)):
    # Add the velocity value to the list
    velocities.append((positions[positionIndex]-positions[positionIndex-1])/float(0.001*timestep*sampleFrequency))
# Converts the list of all velocity values to a numpy array
velocitiesAnalyze = numpy.array(velocities)

# List of all final velocity values and position values to be analyzed
velocitiesFinal = []
positionsFinal = []
# For all segments specified to be analyzed
for startingIndex in startIndex:
    # Determines the upper bound of the segment
    upperBound = startingIndex + boundLength
    # Add all velocity values observed within those bounds
    velocitiesFinal = numpy.append(velocitiesFinal, velocitiesAnalyze[startingIndex:upperBound])
    positionsFinal = numpy.append(positionsFinal, positions[startingIndex:upperBound])

# Calculates the average velocity value that will be reported
finalVelocity = numpy.average(velocitiesFinal)
# Calculates the standard deviation of all measured velocity values
standardError = scipy.stats.sem(velocitiesFinal)
# Calculates the total simulation time included in analysis
analysisTime = (boundLength*len(startIndex)*10*0.001*timestep)

#===============================================================================
# Results
#===============================================================================

# Plots the center of mass positions over time to allow for bounds to be identified
plt.figure(1)
plt.plot(positions)
plt.xlabel('Sample Number')
plt.ylabel('Center of Mass Position (nm)')
plt.title('Center of Mass Position Vs. Time')
if (writeResults is True):
    plt.savefig('poisitions.pdf')
else:                        
    plt.show()

# Plots all position values used in analysis
plt.figure(2)
plt.plot(positionsFinal, linestyle='None', marker='o')
plt.title('Position Values Used in Analysis')
plt.xlabel('Sample Number')
plt.ylabel('Center of Mass Position (nm)')
if (writeResults is True):
    plt.savefig('poisitionsFinal.pdf')
else:                        
    plt.show()

# Plots a histogram of the distribution of velocity values
plt.figure(3)
plt.hist(velocitiesFinal)
plt.title('Velocity Distribution')
if (writeResults is True):
    plt.savefig('velocityDistribution.pdf')
else:                        
    plt.show()

# Prints the average velocity and the standard error of the mean
print('Average Velocity: ' + str(finalVelocity) + ' meters per second')
print('Uncertainty: ' + str(standardError*1.96) + ' meters per second')
print('Time Analyzed: ' + str(analysisTime) + ' nanoseconds')
print('Percent Difference from Check Velocity: ' + str(100*(finalVelocity-checkVelocity)/finalVelocity) + '%')

# If results should be written
if (writeResults is True):
    # Writes contact angle parameters to the results file
    f = open('Velocity.txt','a')
    f.write('\n')
    f.write('\n')
    f.write('Velocity Calculation Results')
    f.write('\n')
    f.write('Starting Index Values: ' + str(startIndex))
    f.write('\n')
    f.write('Segment Length: ' + str(boundLength))
    f.write('\n')
    f.write('Simulation Time Included in Analysis: ' + str(analysisTime) + ' nanoseconds')
    f.write('\n')
    f.write('Average Velocity: ' + str(finalVelocity) + ' meters per second')
    f.write('\n')
    f.write('Standard Error of the Mean: ' + str(standardError) + ' meters per second')
    f.write('\n')
    f.write('Uncertainty: ' + str(1.96*standardError) + ' meters per second')
    f.write('\n')
    f.write('_________________________________________________________________________________________')
    f.write('\n')
    f.write('Velocity Check Data')
    f.write('\n')
    f.write('Average Check Velocity from the Segment Endpoints: ' + str(checkVelocity) + ' meters per second')
    f.write('\n')
    f.write('Percent Difference from Check Velocity: ' + str(100*(finalVelocity-checkVelocity)/finalVelocity) + '%')
    f.write('\n')
    
    # Closes the 'Results.txt' file
    f.flush()
    f.close()
