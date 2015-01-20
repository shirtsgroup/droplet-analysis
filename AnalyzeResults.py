#===============================================================================
# Imports
#===============================================================================

import commands, math, numpy, scipy, os, pdb, shutil, sys, pickle
from scipy import stats

#===============================================================================
# Parameters
#===============================================================================

# File path to the directory
base_path = '/User/username/directory/droplet-analysis/'

#===============================================================================
# Subroutines
#===============================================================================

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


def compute_contact_angle(points, position):
    """
    Computes the contact angle of a subset of points by fitting the points with
    a linear least squares regression and analyzing the slope
    """
    
    # Perform a linear regression and determine the slope of the line
    slope, intercept, r_value, p_value, standard_error = stats.linregress(points)
    
    # Inverts the top right and bottom left slopes so they will be positive for
    # hydrophobic angles and all slopes can be analyzed via the same process
    if ((position == 'Top Right') or (position == 'Bottom Left')):
        slope = -slope
    
    # Handles some common warnings and exceptions resulting from the slope
    # If the line is infinity or not a number
    if ((slope == float("inf")) | (math.isnan(slope))):
        # The contact angle is 90 degrees
        contact_angle = math.degrees(math.pi/2)
    # If the slope is negative infinity
    elif (slope == float("-inf")):
        # The contact angle is -90 degrees
        contact_angle = math.degrees(-math.pi/2)
    # For all other cases
    else:
        # The contact angle is the inverse tangent of the slope converted to
        # degrees. Inverse tangent 2 is used to obtain the correct sign
        contact_angle = math.degrees(math.atan2(slope, 1))
    
    # If the system is hydrophilic the angle will now be less than zero
    if (contact_angle < 0):
        # Inverts the angle to obtain the final contact angle
        contact_angle = -contact_angle
    # Otherwise the angle is hydrophobic and the angle is positive
    else:
        # The contact angle to be returned is the current angle subtracted
        # from 180 degrees in order to comply with typical contact angle
        # measurement conventions
        contact_angle = (180 - contact_angle)
    
    # Returns the computed contact angle
    return contact_angle

def fit_circle(points, guessAngle, advancing):
    """
    Computes the contact angle of an advancing or receding section of the fluid
    with a circular fit and returns the determined contact angle
    """
    
    # Global variables to be used in the regression
    global targetCenter, thickness, maxWallZ, minWallZ, maxWallY, minWallY
    # Stores all points for the system rotated 90 degrees so that a fit can be
    # performed without the need for a multivalued function. Absolute values
    # of all coordinates are used for consistency between advancing and
    # receding point sets
    # Lists to store y and z coordinates respectively (note this is not an
    # error, the system coordinate convention differs from typical horizontal
    # and vertical coordinate conventions)
    xValues = list()
    yValues = list()
    # For each point to be fit with a circular regression
    for point in points:
        # Adds the rotated horizontal coordinate to the appropriate list
        xValues.append(point[1] - (minWallZ + ((maxWallZ - minWallZ)/2)))
        # Adds the rotated vertical coordinate to the appropriate list
        # If the point belongs to a set for an advancing contact angle
        if (advancing):
            # Add the rotated vertical coordinate to the appropriate list
            yValues.append(point[0] - targetCenter)
        # If the point belongs to a set for a receding contact angle
        else:
            # Flip the point about the x axis so the curve is concave down
            # and add the rotated and fliped vertical coordinate to the
            # appropriate list
            yValues.append(targetCenter - point[0])

    # Converts lists to arrays for use in regression analysis
    xArray = numpy.asarray(xValues)
    yArray = numpy.asarray(yValues)
    # Initial guess parameter values for [dropHeight, a] where a = cos(theta)
    p0 = [((targetCenter - minWallY)/2), math.cos(math.radians(guessAngle))]

    # Performs a curve fit to determine the optimal parameter values required to
    # fit the data with a circle via a nonlinear least squared regression
    parameterValues, other = scipy.optimize.leastsq(circleFunctionError, p0, args=(xArray, yArray))
    
    # Determines the contact angle from optimized parameters
    a = parameterValues[1]
    contact_angle = math.degrees(math.acos(a))

    # Returns the computed contact angle
    return contact_angle

def circleFunctionError(parameters, x, y):
    """
    Circular fit objective function to be used in parameter optimization
    Note: This function is horizontally shifted by dropDiameter/2 in order to make sure
    the fit center corresponds to x=0 and the data point center
    """
    
    # Constants and perameters used in evaluating the circular fit function
    dropHeight = parameters[0]
    a = parameters[1]
    dropDiameter = ((maxWallZ - minWallZ) - (2*thickness))
    
    # Residual objective function for use in least squares minimization
    rval = ((dropHeight + (dropDiameter/(2*a))*(numpy.sqrt(1-a**2) - numpy.sqrt(1 - ((a - ((2*a*(x + (0.5*dropDiameter)))/dropDiameter)))**2))) - y)

    # Retun the residual evaluated at the provided x value with the given parameters
    return rval

#===============================================================================
# Main
#===============================================================================

# Adds the 'simulation' directory to the base_path to form the work_path
work_path = base_path + 'simulation/'
# Changes directories to the work_path
os.chdir(work_path)
# Lists to store all advancing and receding contact angles respectively
advancingResults = list()
recedingResults = list()
advancingCircleAngles = list()
recedingCircleAngles = list()
# Opens the file containing the saved and processed data
pickleData = open('Data.pkl','r+')
# Loads processed data for analysis
saveSampleSections, totalDensityBoxes, targetCenter, minWallY, maxWallY, minWallZ, maxWallZ, thickness, boxesX, numWal, numFld, numAtoms, wallAtomsX, wallAtomsY, wallAtomsZ, positions, sigma, sigma_wall, fileSectionSize, density_length, dsigma, dsigmaCircle, timestep, sampleFrequency, fullResults, fileNum = pickle.load(pickleData)
# Closes the 'Data.pkl' file
pickleData.flush()
pickleData.close()

# For all of the obtained sample section density distributions
for densityDistributionSection in saveSampleSections:
    # Process the accumulated density of the section to obtain points
    sectionCircle_advancing, sectionCircle_receding, upper_right, lower_right, upper_left, lower_left, densityThreshold0, numBoxesHigh0, distribution00, distribution0, horizDistribution00, horizDistribution0, vertDistribution00, vertDistribution0 = process(densityDistributionSection, density_length, dsigma, dsigmaCircle, thickness, fileSectionSize, boxesX)
    # Calculates all four contact angles from the points via linear regression
    angleUR = compute_contact_angle(upper_right, 'Top Right')
    angleLR = compute_contact_angle(lower_right, 'Bottom Right')
    angleUL = compute_contact_angle(upper_left, 'Top Left')
    angleLL = compute_contact_angle(lower_left, 'Bottom Left')
    # Add the calculated advancing contact angles to the list of all
    # advancing contact angles
    advancingResults.append(angleUR)
    advancingResults.append(angleLR)
    # Add the calculated receding contact angles to the list of all
    # receding contact angles
    recedingResults.append(angleUL)
    recedingResults.append(angleLL)
    # Calculate an advancing and receding contact angle via a circular
    # fit of the data
    advancingCircleAngle = fit_circle(sectionCircle_advancing, ((angleUR + angleLR)/2), True)
    recedingCircleAngle = fit_circle(sectionCircle_receding, ((angleUL + angleLL)/2), False)
    # Add the calculated angles to a list of all advancing and
    # receding angles respectively
    advancingCircleAngles.append(advancingCircleAngle)
    recedingCircleAngles.append(recedingCircleAngle)

# Process the total density of all sections to obtain points cumulative points
totalCircle_advancing, totalCircle_receding, upper_right_total, lower_right_total, upper_left_total, lower_left_total, densityThreshold, numBoxesHigh, distribution0, distribution, horizDistribution0, horizDistribution, vertDistribution0, vertDistribution = process(totalDensityBoxes, density_length, dsigma, dsigmaCircle, thickness, fileNum, boxesX)
# Calculates all four total contact angles from the points cumulative points
angleUR_total = compute_contact_angle(upper_right_total, 'Top Right')
angleLR_total = compute_contact_angle(lower_right_total, 'Bottom Right')
angleUL_total = compute_contact_angle(upper_left_total, 'Top Left')
angleLL_total = compute_contact_angle(lower_left_total, 'Bottom Left')
# Determines the advancing and receding contact angles
totalAdvancing = float((angleUR_total + angleLR_total)/2)
totalReceding = float((angleUL_total + angleLL_total)/2)
# Calculates a total advancing and receding contact angle via a circular fit
totalCircleAdvancing = fit_circle(totalCircle_advancing, totalAdvancing, True)
totalCircleReceding = fit_circle(totalCircle_receding, totalReceding, False)

# For the linear regression contact angles
# Determines the average advancing contact angle
advancingAvAngle = numpy.mean(advancingResults)
# Determines the standard deviation of the advancing contact angles
advancingStandardDeviation = numpy.std(advancingResults)
# Determines the standard error of the mean advancing contact angle
advancingStandardErrorMean = stats.sem(advancingResults)
# Determines the average receding contact angle
recedingAvAngle = numpy.mean(recedingResults)
# Determines the standard deviation of the receding contact angles
recedingStandardDeviation = numpy.std(recedingResults)
# Determines the standard error of the mean receding contact angle
recedingStandardErrorMean = stats.sem(recedingResults)
# Computes the average contact angle
avContactAngle = float((totalAdvancing + totalReceding)/2)
# Computes the standard error of the mean for the average contact angle
averageStandardErrorMean = float(((advancingStandardErrorMean**2) + (recedingStandardErrorMean**2))**0.5)/2
# Computes the average contact angle hysteresis based on the average advancing
# and receding contact angles
contactAngleHysteresis = (totalAdvancing - totalReceding)
# Determines the standard error of the mean based on the propagation of
# uncertainty from both the advancing and receding contact angle data
hysteresisStandardErrorMean = float((((advancingStandardErrorMean)**2) + ((recedingStandardErrorMean)**2))**0.5)

# For the circular regression contact angles
# Determines the average advancing contact angle
circleAdvancingAvAngle = numpy.mean(advancingCircleAngles)
# Determines the standard deviation of the advancing contact angles
circleAdvancingStandardDeviation = numpy.std(advancingCircleAngles)
# Determines the standard error of the mean advancing contact angle
circleAdvancingStandardErrorMean = stats.sem(advancingCircleAngles)
# Determines the average receding contact angle
circleRecedingAvAngle = numpy.mean(recedingCircleAngles)
# Determines the standard deviation of the receding contact angles
circleRecedingStandardDeviation = numpy.std(recedingCircleAngles)
# Determines the standard error of the mean receding contact angle
circleRecedingStandardErrorMean = stats.sem(recedingCircleAngles)
# Computes the average contact angle
circleAvContactAngle = float((totalCircleAdvancing + totalCircleReceding)/2)
# Computes the standard error of the mean for the average contact angle
circleAverageStandardErrorMean = float(((circleAdvancingStandardErrorMean**2) + (circleRecedingStandardErrorMean**2))**0.5)/2
# Computes the average contact angle hysteresis based on the average advancing
# and receding contact angles
circleContactAngleHysteresis = (totalCircleAdvancing - totalCircleReceding)
# Determines the standard error of the mean based on the propagation of
# uncertainty from both the advancing and receding contact angle data
circleHysteresisStandardErrorMean = float((((circleAdvancingStandardErrorMean)**2) + ((circleRecedingStandardErrorMean)**2))**0.5)

# Holds all velocity values determined from center of mass position
# analysis
velocities = list()
# Holds the center of mass position of the previous sample
comPrevious = targetCenter
# Determined the center of mass velocity at every sample
# For every center of mass position
for comPosition in positions:
    # Determines the velocity of the center of mass at that sample
    comVelocity = abs(float(comPosition - comPrevious)/float(0.001*timestep*sampleFrequency))
    # Sets the previous center of mass position to the current position
    comPrevious = comPosition
    # Adds the determined velocity value to a list of all velocity values
    velocities.append(comVelocity)

# Note: this reports the velocity magnitude and is not capable of
# determining velocity directions. Since the center of mass velocity
# is always in the positive length dimension this is typically not an
# issue. However, in the case of no center of mass acceleration,
# the small changes in center of mass position will all be calculated as
# positive. This will result in a center of mass velocity slightly greater
# than zero as movements in the negative direction will be calculated as
# positive. Also, due to the presence of strong outliers, the median
# center of mass velocity should always be reported and not the mean.
                                
# Writes contact angle parameters to the results file
f = open('Results.txt','a')
f.write('\n')
f.write('\n')
f.write('Contact Angle Calculation Parameters')
f.write('\n')
f.write('Sigma (fluid): ' + str(sigma) + ' nanometers')
f.write('\n')
f.write('Wall Sigma: ' + str(sigma_wall) + ' nanometers')
f.write('\n')
f.write('Side Length of a Sampling Box: ' + str(density_length) + ' Wall Sigma')
f.write('\n')
f.write('Fraction of Total Fluid Height Analyzed Above Each Fluid-Wall Contact Plane with Each Linear Regression: ' + str(dsigma))
f.write('\n')
f.write('Note: This height range begins two sigmas above the plane of wall-fluid interaction to avoid unpredictable behavior at the interface')
f.write('\n')
f.write('Fraction of Total Fluid Height Analyzed with Each Circular Regression: ' + str(dsigmaCircle))
f.write('\n')
f.write('Note: This height range is centered around the midplane of the walls and does not include the 2 sigma above each wall that are not considered in regression analysis')
f.write('\n')
f.write('Number of Sampling Boxes Fit with a Linear Regression to Obtain Each Contact Angle: ' + str(numBoxesHigh) + ' Boxes')
f.write('\n')
f.write('Note: This value plus the number of boxes comprising the two sigmas above each wall is equivalent to the number of boxes excluded from circular regression analysis')
f.write('\n')
f.write('Modified Density Threshold for a Density Box to be Considered Fluid: ' + str(densityThreshold))
f.write('\n')
f.write('Note: This is equal to one half the median adjusted density of a sampling box. The density distribution has been divided by the number of system samples in a sample section, the number of sampling boxes in the depth dimension of the system, and the volume of a density box so that the density distribution and cutoff threshold are roughly independent of contact angle calculation parameters.')
f.write('\n')
f.write('Number of Samples Analyzed: ' + str(fileNum) + ' Samples')
f.write('\n')
f.write('Number of Samples per Sample Section: ' + str(fileSectionSize) + ' Samples')
f.write('\n')
f.write('Number of Sample Sections Analyzed: ' + str(int(fileNum/fileSectionSize)) + ' Sections')
f.write('\n')
f.write('Note: Two advancing and two receding angles are obtained from each sample section via a linear regression and a single advancing and a single receding contact angle is obtained from each sample section via a circular regression')
f.write('\n')
f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes results to the results file so they can be analyzed
f.write('Linear Regression Contact Angle Data | Units are degrees')
f.write('\n')
f.write('Cumulative Advancing Contact Angle: ' + str(totalAdvancing))
f.write('\n')
f.write('Average Advancing Sample Section Contact Angle: ' + str(advancingAvAngle))
f.write('\n')
f.write('Advancing Angle Standard Deviation: ' + str(advancingStandardDeviation))
f.write('\n')
f.write('Advancing Angle Standard Error in the Mean: ' + str(advancingStandardErrorMean))
f.write('\n')
f.write('Cumulative Receding Contact Angle: ' + str(totalReceding))
f.write('\n')
f.write('Average Receding Sample Section Contact Angle: ' + str(recedingAvAngle))
f.write('\n')
f.write('Receding Angle Standard Deviation: ' + str(recedingStandardDeviation))
f.write('\n')
f.write('Receding Angle Standard Error in the Mean: ' + str(recedingStandardErrorMean))
f.write('\n')
f.write('Average Contact Angle: ' + str(avContactAngle))
f.write('\n')
f.write('Average Contact Angle Standard Error in the Mean: ' + str(averageStandardErrorMean))
f.write('\n')
f.write('Average Contact Angle Hysteresis: ' + str(contactAngleHysteresis))
f.write('\n')
f.write('Average Contact Angle Hysteresis Standard Error in the Mean: ' + str(hysteresisStandardErrorMean))
f.write('\n')
f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes the final results with the correct 95% confidence interval
f.write('Linear Regression Contact Angle Results with a 95% Confidence Interval')
f.write('\n')
f.write('Advancing Contact Angle: ' + str(totalAdvancing) + ' Degrees Plus or Minus ' + str(1.96*advancingStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('Receding Contact Angle: ' + str(totalReceding) + ' Degrees Plus or Minus ' + str(1.96*recedingStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('Contact Angle Hysteresis: ' + str(contactAngleHysteresis) + ' Degrees Plus or Minus ' + str(1.96*hysteresisStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('Average Contact Angle: ' + str(avContactAngle) + ' Degrees Plus or Minus ' + str(1.96*averageStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes advancing contact angles to the results file so they can be analyzed
f.write('Linear Regression Advancing Angles: ')
f.write('\n')
f.write('\n')
for angle in advancingResults:
    f.write(str(angle))
    f.write('\n')

f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes receding contact angles to the results file so they can be analyzed
f.write('Linear Regression Receding Angles: ')
f.write('\n')
f.write('\n')
for angle in recedingResults:
    f.write(str(angle))
    f.write('\n')

f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes results to the results file so they can be analyzed
f.write('Circular Regression Contact Angle Data | Units are degrees')
f.write('\n')
f.write('Cumulative Advancing Contact Angle: ' + str(totalCircleAdvancing))
f.write('\n')
f.write('Average Advancing Sample Section Contact Angle: ' + str(circleAdvancingAvAngle))
f.write('\n')
f.write('Advancing Angle Standard Deviation: ' + str(circleAdvancingStandardDeviation))
f.write('\n')
f.write('Advancing Angle Standard Error in the Mean: ' + str(circleAdvancingStandardErrorMean))
f.write('\n')
f.write('Cumulative Receding Contact Angle: ' + str(totalCircleReceding))
f.write('\n')
f.write('Average Receding Sample Section Contact Angle: ' + str(circleRecedingAvAngle))
f.write('\n')
f.write('Receding Angle Standard Deviation: ' + str(circleRecedingStandardDeviation))
f.write('\n')
f.write('Receding Angle Standard Error in the Mean: ' + str(circleRecedingStandardErrorMean))
f.write('\n')
f.write('Average Contact Angle: ' + str(circleAvContactAngle))
f.write('\n')
f.write('Average Contact Angle Standard Error in the Mean: ' + str(circleAverageStandardErrorMean))
f.write('\n')
f.write('Average Contact Angle Hysteresis: ' + str(circleContactAngleHysteresis))
f.write('\n')
f.write('Average Contact Angle Hysteresis Standard Error in the Mean: ' + str(circleHysteresisStandardErrorMean))
f.write('\n')
f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes the final results with the correct 95% confidence interval
f.write('Circular Regression Contact Angle Results with a 95% Confidence Interval')
f.write('\n')
f.write('Advancing Contact Angle: ' + str(totalCircleAdvancing) + ' Degrees Plus or Minus ' + str(1.96*circleAdvancingStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('Receding Contact Angle: ' + str(totalCircleReceding) + ' Degrees Plus or Minus ' + str(1.96*circleRecedingStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('Contact Angle Hysteresis: ' + str(circleContactAngleHysteresis) + ' Degrees Plus or Minus ' + str(1.96*circleHysteresisStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('Average Contact Angle: ' + str(circleAvContactAngle) + ' Degrees Plus or Minus ' + str(1.96*circleAverageStandardErrorMean) + ' Degrees')
f.write('\n')
f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes advancing contact angles to the results file so they can be analyzed
f.write('Circular Regression Advancing Angles: ')
f.write('\n')
f.write('\n')
for angle in advancingCircleAngles:
    f.write(str(angle))
    f.write('\n')

f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes receding contact angles to the results file so they can be analyzed
f.write('Circular Regression Receding Angles: ')
f.write('\n')
f.write('\n')
for angle in recedingCircleAngles:
    f.write(str(angle))
    f.write('\n')

f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Writes the center of mass velocity from one time step to the next
f.write('Fluid Center of Mass Average Velocity: ' + str(float(numpy.mean(velocities))) + ' Nanometers per Nanosecond')
f.write('\n')
f.write('Fluid Center of Mass Median Velocity: ' + str(float(numpy.median(velocities))) + ' Nanometers per Nanosecond')
f.write('\n')
f.write('Fluid Center of Mass Velocity Standard Deviation: ' + str(float(numpy.std(velocities))) + ' Nanometers per Nanosecond')
f.write('\n')
f.write('Fluid Center of Mass Velocity Standard Error in the Mean: ' + str(float(stats.sem(velocities))) + ' Nanometers per Nanosecond')

f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# If full results are desired
if (fullResults):
    # Writes simulation conditions to the results file
    inFile1 = open('minimize.mdp', 'r')
    linesin1 = inFile1.readlines()
    inFile2 = open('nvt.mdp', 'r')
    linesin2 = inFile2.readlines()
    inFile3 = open('simulation.mdp', 'r')
    linesin3 = inFile3.readlines()
    f.write('Simulation Conditions')
    f.write('\n')
    f.write('\n')
    f.write('Minimization Conditions: ')
    f.write('\n')
    for line in linesin1:
        f.write(str(line))
        f.write('\n')
    f.write('\n')
    f.write('NVT Ensemble Conditions: ')
    f.write('\n')
    for line in linesin2:
        f.write(str(line))
        f.write('\n')
    f.write('\n')
    f.write('Simulation Conditions: ')
    f.write('\n')
    for line in linesin3:
        f.write(str(line))
        f.write('\n')

    f.write('\n')
    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the subsample two dimensional matrix containing the box densities
    f.write('Sample Section Two Dimensional Complete Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(distribution00.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the subsample two dimensional matrix containing the box densities
    # filtered to remove insufficient density boxes
    f.write('Sample Section Two Dimensional Filtered Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(distribution0.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the subsample one dimensional array containing the box densities in the
    # horizontal direction
    f.write('Sample Section Horizontal One Dimensional Complete Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(horizDistribution00.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the cumulative two dimensional matrix containing the box densities
    f.write('Cumulative Two Dimensional Complete Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(distribution0.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the cumulative two dimensional matrix containing the box densities
    # filtered to remove insufficient density boxes
    f.write('Cumulative Two Dimensional Filtered Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(distribution.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the cumulative one dimensional array containing the box densities in
    # the horizontal direction
    f.write('Cumulative Horizontal One Dimensional Complete Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(horizDistribution0.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the cumulative one dimensional array containing the box densities
    # filtered to remove insufficient density boxes in the horizontal direction
    f.write('Cumulative Horizontal One Dimensional Filtered Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(horizDistribution.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the cumulative one dimensional array containing the box densities in
    # the vertical direction
    f.write('Cumulative Vertical One Dimensional Complete Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(vertDistribution0.tolist()))

    f.write('\n')
    f.write('________________________________________________________________________________')
    f.write('\n')
    f.write('\n')

    # Writes the cumulative one dimensional array containing the box densities
    # filtered to remove insufficient density boxes in the vertical direction
    f.write('Cumulative Vertical One Dimensional Filtered Distribution')
    f.write('\n')
    f.write('\n')
    numpy.set_printoptions(threshold='nan', linewidth='nan')
    f.write(str(vertDistribution.tolist()))

f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')

# Closes the 'Results.txt' file
f.flush()
f.close()
