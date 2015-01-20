
#===============================================================================
# Imports
#===============================================================================

import commands, math, os, pdb, scipy, shutil, sys
from random import randrange

#===============================================================================
# Parameters
#===============================================================================

# Sets the base path to be equal to the path to the 'superhydrophobic' directory
base_path = '/User/username/directory/superhydrophobic/'

# General Constants
# Avogadro's number (DO NOT CHANGE!)
NA = 6.0221413*10**23
# Nanometer distance at which inter-particle Lennard-Jones Potential is zero
# Note: Either specify sigma_wall or reducedSigmaWF and comment the other
sigma_fluid = 0.3045
# Comment out sigma_wall if a reduced unit is specified
# Note: even if chemical roughness is used, this sigma will be the value used
# to space all wall atoms, which will be one sigma from their closest neighbor
sigma_wall = 0.2705
# Comment out both of the following lines if a sigma_wall is explicitly given
# reducedSigmaWF = 0.5
# sigma_wall = ((2*sigma_fluid*reducedSigmaWF) - sigma_fluid)
# sigma(ij) = (sigma(i)+sigma(j))/2
# reduced sigma(wall-fluid) = (sigma(wall)+sigma(fluid))/(2*sigma_fluid)

# Fluid
# Spacing factor to set the space between particles, factor*sigma = spacing
# This value will be multiplied by sigma in the code so do not do it here
# Minimum energy spacing for staggered atoms is 2**(1/6)
# Typically need to multiply by a factor of 1.15 or so to prevent initial
# fluid configuration from splitting apart. A large value will cause the
# fluid to contract and separate from the walls
factor = 1.15
# Fluid epsilon: depth of the fluid potential wells (kJ/mol)
fluid_epsilon = 1.67*10**(-24)*NA
# Fraction of wall in the y dimension that will initially be filled with fluid
fluidFillFraction_y = 0.5
# Mass of fluid in kg
fluid_mass = 6.6*10**(-26)*1000*NA

# Walls
# Lattice parameter for face centered cubic wall configuration
latticeParameter = 0.39236
# X, Y, and Z dimensions of the system walls
wall_x = 35*latticeParameter
wall_y = 70*latticeParameter
# Spacing between the innermost roughness of the two walls
dz = 25*latticeParameter
# Wall epsilon: depth of the wall potential wells (kJ)
# Note: Either specify wall_epsilon or reducedEpsilonWF and comment the other
# Reduced epsilon wall-fluid
# epsilon(ij) = sqrt(epsilon(ii)*epsilon(jj))
# Reduced epsilon(wall-fluid)=sqrt(epsilon(fluid)*epsilon(wall))/epsilon(fluid)
# COMMENT OUT IF wall_epsilon IS EXPLICITLY SPECIFIED
reducedEpsilonWF = 0.4
# COMMENT OUT IF wall_epsilon IS EXPLICITLY SPECIFIED
wall_epsilon = ((reducedEpsilonWF)**2)*fluid_epsilon
# COMMENT OUT ALL wall_epsilon VALUES BELOW HERE IF A REDUCED VALUE IS SPECIFIED
# ONLY ONE wall_epsilon CAN BE SPECIFIED
# Hydrophobic potential well depth
# wall_epsilon = 0.135*10**(-24)*NA
# Hydrophilic potential well depth
# wall_epsilon = 0.835*10**(-24)*NA
# Mass of wall in kg
wall_mass = 3.24*10**(-26)*1000*NA

# Dynamic Contact Angles
# Results in dynamic fluid motion in the positive Y dimension if set to True
dynamic = False
# Specifies the fluid acceleration in nanometers per picosecond squared
acceleration = 0.001
# Specified the temperature in Kelvin
temperature = 85

# Roughness
# Boolean to turn the system roughness on or off
# Note: Any type of roughness requires this boolean to be true
# To generate a chemically but not physically rough the height to zero
roughness = True
# Number of solid wall layers to be generated regardless of roughness type
baseLayers = 5
# Type of surface feature used to generate roughness
# Options: Vertical Box, Cylinder, Cone, Horizontal Box X, Horizontal Box Y,
# Diagonal Box, Sphere, Ring, Prism and Grid. Any other value for roughnessType
# will result in walls height atoms thick and shaped as parallel planes.
roughnessType = 'Cylinder'
# Atom height of pillars
height = 0 #5
# Atom radius of pillars
radius = 2.5
# Atom spacing between pillars
spacing = 8

# Note: Both roughness and arbitrary physical roughness must be turned on for
# arbitrary physical roughness to be generated
# If arbitrary chemical or physical roughness is specified, the conditional
# corresponding to each will be executed and values inside each will be used
# Turns arbitrary physical roughness on or off
arbitraryPhysicalRoughness = False
# NOTE: THE RADIUS AND SPACING MUST BE SET SUCH THAT NO PILLARS OVERLAP! IF
# PILLARS OVERLAP, MULTIPLE ATOMS WILL BE GENERATED IN THE SAME LOCATION AND
# THE SIMULATION WILL FAIL.
# If there is to be arbitrary physical roughness
if (arbitraryPhysicalRoughness == True):
    # Sets the range of random pillar height values
    # Maximum possible height (in number of atoms)
    height = 8
    # Minimum possible height (in number of atoms)
    heightMin = 0
    # Sets the range of random pillar radius values
    # Maximum possible radius (in number of atoms)
    radiusMax = 3
    # Minimum possible radius (in number of atoms)
    radiusMin = 2
    # Sets the range of random pillar spacing values
    # Maximum possible spacing (in number of atoms)
    spacingMax = 10
    # Minimum possible spacing (in number of atoms)
    spacingMin = 7
else:
    heightMin = 'N/A'
    radiusMax = 'N/A'
    radiusMin = 'N/A'
    spacingMax = 'N/A'
    spacingMin = 'N/A'
    
# Turns chemical roughness on or off
chemicalRoughness = False
# Turns arbitrary chemical roughness on or off
arbitraryChemicalRoughness = False

roughAtoms = list()
# Enter atom information in the following format:
# [residue name, atom name, sigma, epsilon, mass, charge, minimum percentage, maximum percentage]
# A random number will be generated for the atoms of each pillar that is
# between 1 and 100. If the number generated is among the range specified,
# that atom type will be generated, otherwise, a default wall atom will be
# generated. Append each additional atom type to the roughAtoms list.
# Note: Number ranges of atoms CAN NOT OVERLAP
# Atomname must begin with 'WAL' and be no more than five characters long
roughAtoms.append(['WALa', 'WALa', 0.2000, 0.135*10**(-24)*NA, 3.24*10**(-26)*1000*NA, 0, 1, 33])
roughAtoms.append(['WALb', 'WALb', 0.2000, 0.070*10**(-24)*NA, 3.24*10**(-26)*1000*NA, 0, 33, 66])

"""
    mimimize.txt output options
    Note: File must end in 0 to exit the menu and continue the simulation
    Put every number on a separate line
    
    Select the terms you want from the following list by
    selecting either (part of) the name or the number or a combination.
    End your selection with an empty line or a zero.
    -------------------------------------------------------------------
    1  LJ-(SR)          2  Disper.-corr.    3  Coulomb-(SR)     4  Potential
    5  Pres.-DC         6  Pressure         7  Vir-XX           8  Vir-XY
    9  Vir-XZ          10  Vir-YX          11  Vir-YY          12  Vir-YZ
    13  Vir-ZX          14  Vir-ZY          15  Vir-ZZ          16  Pres-XX
    17  Pres-XY         18  Pres-XZ         19  Pres-YX         20  Pres-YY
    21  Pres-YZ         22  Pres-ZX         23  Pres-ZY         24  Pres-ZZ
    25  #Surf*SurfTen   26  T-rest
    
    
    
    nvt.txt output options
    Note: File must end in 0 to exit the menu and continue the simulation
    Put every number on a separate line
    
    Select the terms you want from the following list by
    selecting either (part of) the name or the number or a combination.
    End your selection with an empty line or a zero.
    -------------------------------------------------------------------
    1  LJ-(SR)          2  Disper.-corr.    3  Coulomb-(SR)     4  Potential
    5  Kinetic-En.      6  Total-Energy     7  Temperature      8  Pres.-DC
    9  Pressure        10  Vir-XX          11  Vir-XY          12  Vir-XZ
    13  Vir-YX          14  Vir-YY          15  Vir-YZ          16  Vir-ZX
    17  Vir-ZY          18  Vir-ZZ          19  Pres-XX         20  Pres-XY
    21  Pres-XZ         22  Pres-YX         23  Pres-YY         24  Pres-YZ
    25  Pres-ZX         26  Pres-ZY         27  Pres-ZZ         28  #Surf*SurfTen
    29  T-system        30  Lamb-system
    
    
    Note: simulation.txt should not be changed
"""

#===============================================================================
# Subroutines
#===============================================================================


def write_file(filename, lines):
    """
    Writes the desired file to the work path when given the name of the file to
    be written and the lines 
    """
    
    # Opens the desired file
    outfile = open(filename,'w+')
    # For every line to be written to that file
    for line in lines:
        # Writes the line to the file
        outfile.write(line)
    # Closes the file
    outfile.close()
    
    # Void return
    return


def generate_line(resnum, resname, atomname, atomnum, x, y, z):
    """
    Generates a line used to represent a single atom in a .gro file
    """
    
    # Generates the correctly formatted line
    line = '%(resnum)5d%(resname)-4s%(atomname)6s%(atomnum)5d%(x)8.3f%(y)8.3f%(z)8.3f\n'%vars()
    
    # Returns the line generated from the atom information
    return line


def chemically_rough_line(resname, atomname, pattern):
    """
    Randomly selects from the provided atoms with the probabilities provided
    in order to return a resname and atomname for the atom type. Lines required
    to add this type of atom to the Gromacs topology file are also generated
    """
    
    # If there is arbitrary chemical roughness in the system
    if (arbitraryChemicalRoughness):
        # Generates a random number that will be used to select an atom type
        genNum = randrange(1, 101)

        # Boolean to indicate that no non default atom has been selected
        selected = False
        # For every type of alternate atom specified by roughList
        for atom in roughAtoms:
            # If the random number is greater than or equal to the lower end of
            # the range and less than the higher end of the range
            if ((genNum >= atom[6]) and (genNum < atom[7])):
                # Return the residue name and atom name
                result = [atom[0], atom[1]]
                # Selected becomes true as a non default atom has been selected
                selected = True
        # If a non default atom has not been selected
        if (not selected):
            # Return the default atom type
            result = [resname, atomname]
    # Otherwise there is patterned chemical roughness in the system
    elif (chemicalRoughness):
        # Obtains the number of options available for different atoms
        options = 1 + len(roughAtoms)
        # Determines a number from which the atom will be selected
        selection = int(pattern%options)
        # If selection is among the index range of roughAtoms
        if (selection < len(roughAtoms)):
            # Select the atom type based on the selection number
            atomSelectionType = roughAtoms[selection]
            # Determine the resname and atomname based on the selected atom type
            resname = atomSelectionType[0]
            atomname = atomSelectionType[1]
            # Set the result equal to the obtained information
            result = [resname, atomname]
        # Otherwise the default wall atom type is returned
        else:
            # Return the default atom type
            result = [resname, atomname]
    # Otherwise there is no chemical roughness
    else:
        # Return the default atom type
        result = [resname, atomname]

    # String to return the lines needed for additional atoms to be added
    # to the topology file
    roughAtomsTop = ''
    # String to add the additional atom types to the topology file
    atomTypes = ''
    # Atom number begins at 2 as the default wall atoms are already present
    aNum = 2
    # For every atom type specified by the user
    for atom in roughAtoms:
        # Obtain required atom information and format it for the topology file
        atom_name = atom[1].ljust(12)
        atom_name2 = atom[1].rjust(5)
        atom_name3 = atom[1].ljust(4)
        atom_mass = atom[4]
        atom_sigma = atom[2]
        atom_epsilon = atom[3]
        atom_charge = atom[5]
        atom_number = str(1).rjust(6)
        # Generates a line for the atom type to be added to the topology file
        line1 = '%(atom_name)s35      %(atom_mass)1.2s    0.0000  A   %(atom_sigma)1.5e  %(atom_epsilon)1.5e'%vars()
        # Adds the atom type to the topology string and begins a new line
        roughAtomsTop = roughAtomsTop + str(line1)
        roughAtomsTop = roughAtomsTop + str('\n')
        # Generates a line for the atom to be added to the topology file
        line2 = """
[ moleculetype ]
; name  nrexcl
%(atom_name)s 1

[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass
%(atom_number)s  %(atom_name2)s   1   %(atom_name3)s%(atom_name2)s    1%(atom_charge)13.6f%(atom_mass)13.5f"""%vars()
        # Adds the atom type to the topology file string and begins a new line
        atomTypes = atomTypes + str(line2)
        atomTypes = atomTypes + str('\n')
                
    # Returns the residue name and atom name for generating the atom and the
    # atom type and atom lines to be added to the topology file
    return result, roughAtomsTop, atomTypes


def generate_pillar(x, y, pillars, level, radius):
    """
    Used by generate_wall() to determine where the wall should be more than one
    layer thick. This is used exclusively for generating surface roughness
    geometries on both walls
    """
    
    # By default, only the outermost layer of wall is generated
    generation = False

    # Variable to store the type of atom located at that pillar
    atomType = 'null'
    
    # If the x and y location meet the conditions for additional layers
    # of wall to be generated, the boolean will become true and is returned,
    # otherwise the boolean will remain false and is returned. Numerous
    # roughness patterns are listed below.

    # Height must always be more than one and problems may result from a low
    # height. The spacing between the walls (dz) will not change as the wall
    # positions will adjust to account for different thickness. A tested and
    # working value for height, radius, and spacing are recommended for each
    # geometry along with notes about some of the specific constraints in
    # the parameters for each geometry. 
    
    # Vertical cylinders
    # Recommended: Height = 6; Radius = 2.5; Spacing = 9
    # Note: Spacing should be greater than or equal to twice the radius
    if (roughnessType == 'Cylinder'):
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within radius of a pillar center
            if ((((x/latticeParameter) - pillar[0])**2 + (((y/latticeParameter) - pillar[1]))**2) < radius**2):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]

    # Diagonal Box
    # Recommended: Height = 6; Radius = 3; Spacing = 10
    # Note: Spacing should be greater than twice the radius
    elif (roughnessType == 'Diagonal Box'):
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within the desired range
            if ((((x/latticeParameter) - pillar[0]) + (((y/latticeParameter) - pillar[1])))**2 < radius**2):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]

    # Horizontal Box Y
    # Recommended: Height = 6; Radius = 3; Spacing = 5
    # Note: Spacing should be greater than or equal to the radius
    elif (roughnessType == 'Horizontal Box Y'):
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within a radius of the pillar x value
            if ((((x/latticeParameter) - pillar[0]) > (pillar[0] - radius)) and (((x/latticeParameter) - pillar[0]) < (pillar[0] + radius))):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]

    # Horizontal Box X
    # Recommended: Height = 6; Radius = 3; Spacing = 5
    # Note: Spacing should be greater than or equal to the radius
    elif (roughnessType == 'Horizontal Box X'):
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within a radius of the pillar y value
            if ((((y/latticeParameter) - pillar[1]) > (pillar[1] - radius)) and (((y/latticeParameter) - pillar[1]) < (pillar[1] + radius))):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]
                   
    # Grid
    # Recommended: Height = 6; Radius = 3; Spacing = 5
    # Note: Spacing should be greater than or equal to the radius
    elif (roughnessType == 'Grid'):
        # For every pillar center
        for pillar in pillars:
            # If the atom is within a radius of the pillar center in x or y
            if (((((x/latticeParameter) - pillar[0]) > (pillar[0] - radius)) and (((x/latticeParameter) - pillar[0]) < (pillar[0] + radius))) | ((((y/latticeParameter) - pillar[1]) > (pillar[1] - radius)) and (((y/latticeParameter) - pillar[1]) < (pillar[1] + radius)))):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]

    # Prism Y
    # Recommended: Height = 5; Radius = 5; Spacing = 3
    elif (roughnessType == 'Prism Y'):
        # Calculated the correct radius based on the fluid level
        radiusMod = radius - (radius*level/height)
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within a radius of the pillar y value
            if ((((x/latticeParameter) - pillar[0]) > (pillar[0] - radiusMod)) and (((x/latticeParameter) - pillar[0]) < (pillar[0] + radiusMod))):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]

    # Prism X
    # Recommended: Height = 5; Radius = 5; Spacing = 3
    elif (roughnessType == 'Prism X'):
        # Calculated the correct radius based on the fluid level
        radiusMod = radius - (radius*level/height)
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within a radius of the pillar y value
            if ((((y/latticeParameter) - pillar[1]) > (pillar[1] - radiusMod)) and (((y/latticeParameter) - pillar[1]) < (pillar[1] + radiusMod))):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]
    
    # Vertical Box
    # Recommended: Height = 5; Radius = 4; Spacing = 5
    # Note: Spacing should be greater than or equal to the radius
    elif (roughnessType == 'Vertical Box'):
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within a radius from the pillar center in
            # both x and y (Note: square will have a side length of 2*radius)
            if ((((x/latticeParameter) - pillar[0]) > (pillar[0] - radius)) and (((x/latticeParameter) - pillar[0]) < (pillar[0] + radius)) and (((y/latticeParameter) - pillar[1]) < (pillar[1] + radius)) and (((y/latticeParameter) - pillar[1]) > (pillar[1] - radius))):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]
    
    # Cone
    # Recommended: Height = 6; Radius = 6; Spacing = 7
    # Note: Spacing should be greater than or equal to the radius
    elif (roughnessType == 'Cone'):
        # Calculate the correct radius based on the level of the wall layer
        radiusMod = radius - ((level*radius)/height)
        # For every pillar center
        for pillar in pillars:
            # If the atom is within a modified radius of a pillar center
            if ((((x/latticeParameter) - pillar[0])**2 + (((y/latticeParameter) - pillar[1]))**2) < radiusMod**2):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]

    # Sphere
    # Recommended: Height = 6; Radius = 6; Spacing = 14
    # Note: Spacing should be greater than or equal to two times the radius
    elif (roughnessType == 'Sphere'):
        # Calculate the correct radius based on the level of the wall layer
        radiusMod = (radius**2 - level**2)**0.5
        # For every pillar center
        for pillar in pillars:
            # If the atom is within a modified radius of a pillar center
            if ((((x/latticeParameter) - pillar[0])**2 + (((y/latticeParameter) - pillar[1]))**2) < radiusMod**2):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]

    # Ring
    # Recommended: Height = 5; Radius = 10; Spacing = 20
    # Note: Spacing should be greater than or equal to two times the radius
    elif (roughnessType == 'Ring'):
        # Ratio of outer to inner radius
        # Note: Radius sets outer radius so this should be greater than one
        radiusRatio = 2
        # Inner radius of the ring
        inRadius = radius/radiusRatio
        # For every pillar center
        for pillar in pillars:
            # If the atom location is within the outer radius of a pillar center
            # and outside the inner radius
            if (((((x/latticeParameter) - pillar[0])**2 + (((y/latticeParameter) - pillar[1]))**2) < radius**2) and ((((x/latticeParameter) - pillar[0])**2 + (((y/latticeParameter) - pillar[1]))**2) > inRadius**2)):
                # Return true
                generation = True
                # Assigns atom type to be equal to the atom type for that pillar
                atomType = pillar[2]
    
    # Return a boolean to signify if a wall atom should be generated at a
    # particular x and y location based on the surface pattern used to generate
    # the desired level of roughness 
    return generation, atomType


def generate_wall(x0, y0, z0, natoms, location, pillarLocations, outsideWall, level, radius, resname, atomname):
    """
    Used by generate_wall_gro() to generate all wall atoms for the x and y
    locations of a specific z value 
    """

    # Creates a list to store the lines representing the wall atoms
    lines = list()
    # Integer to signal if the atom number has exceeded 99999 and has been reset
    counterReset = 0 
    # Integers to represent atoms of wall length in the x and y
    n = int(wall_x/latticeParameter)
    m = int(wall_y/latticeParameter)
    # Accumulator to represent a bottom wall (0) and top wall (1) in resnum
    l = 0
    # Sets initial z value to z0
    z = z0
    # If the wall is a top wall
    if (location == 'top'):
        # z is increased by dz so there is dz spacing between the walls
        z = (z + dz)
        # l is changed to 1 to reflect the fact that the wall is a top wall
        l = 1

    # If there is roughness in the wall
    if (roughness):
        # Iterates through all of the x values where wall will be generated 
        for i in range(n):
            # Iterates through all of the y values where wall will be generated
            for j in range(m):
                # Determines the x location of the current atom
                x = x0 + (i*latticeParameter)
                # Determines the y location of the current atom
                y = y0 + (j*latticeParameter)
                # Obtain results to see if a pillar should be generated
                generatePillar, typeA = generate_pillar(x, y, pillarLocations, level, radius)
                # If the wall layer is an outside one or the condition for
                # additional layers is true
                if (outsideWall):
                    # Assigns resnum to be equal to the current i + j + l
                    resnum = (i + j + 1)
                    # Increments natoms to reflect the addition of a wall atom
                    natoms = (natoms + 1)
                    # While natoms is greater than 5 digits it is reset
                    while (natoms > 99999):
                        # Resets natoms
                        natoms = (natoms - 100000)
                        # Increments counterReset
                        counterReset = (counterReset + 1)
                    # Generate the line to represent the current wall atom
                    line = generate_line(resnum, 'WAL', 'WAL', natoms, x, y, z)
                    # Adds the current line to a list containing all other lines
                    lines.append(line)
                # If this is a pillar location and not an outside wall location
                elif (generatePillar and not outsideWall):
                    # Assigns resnum to be equal to the current i + j + l
                    resnum = (i + j + 1)
                    # Increments natoms to reflect the addition of a wall atom
                    natoms = (natoms + 1)
                    # While natoms is greater than 5 digits it is reset
                    while (natoms > 99999):
                        # Resets natoms
                        natoms = (natoms - 100000)
                        # Increments counterReset
                        counterReset = (counterReset + 1)
                    # Generate the line to represent the current wall atom
                    line = generate_line(resnum, resname, atomname, natoms, x, y, z)
                    # Adds the current line to a list containing all other lines
                    lines.append(line)

    # If there is not roughness in the wall, generate the wall atoms
    else:
        for i in range(n):
            # Iterates through all of the y values where wall will be generated
            for j in range(m):
                # Determines the x location of the current atom
                x = x0 + (i*latticeParameter)
                # Determines the y location of the current atom
                y = y0 + (j*latticeParameter)
                # Assigns resnum to be equal to the current i + j + l
                resnum = (i + j + 1)
                # Increments natoms to reflect the addition of a wall atom
                natoms = (natoms + 1)
                # While natoms is greater than 5 digits it is reset
                while (natoms > 99999):
                    # Resets natoms
                    natoms = (natoms - 100000)
                    # Increments counterReset
                    counterReset = (counterReset + 1)
                # Generates the line to represent the current wall atom
                line = generate_line(resnum, resname, atomname, natoms, x, y, z)
                # Adds the current line to the list containing all other lines
                lines.append(line)

    # Returns the list of lines used to represent all atoms within the wall and
    # the number of times the counter has been reset if any
    return lines, counterReset


def generate_wall_gro(radius):
    """
    Generates the system wall(s) in the system.gro file by using generate_wall()
    at every desired z value
    """
    
    # Accumulator for the number of wall atoms added to the system
    natoms = 0
    # Starting coordinates for wall creation in the x, y, and z directions
    x0 = ((box_x - wall_x)/2)
    y0 = ((box_y - wall_y)/2)
    z0 = (box_z - (wall_z*nwall) - dz)/2
    # List to store lines for all wall atoms that are generated
    lines = list()
    # List to store location of all pillar centers
    pillarsTop = list()
    pillarsBottom = list()
    # Booleans for the top and bottom walls to signify outermost layers
    # These layers will be fully generated regardless pillar location
    bottomWallOutside = True
    topWallOutside = True
    # Sets the initial wall level from the outside wall to be used for roughness
    levelBottom = 1
    levelTop = 1
    # Counter to alternate among atoms for regular chemical roughness
    # Initially set to zero
    pattern = 0

    # Pillar generation
    # If there is to be roughness in the system
    if (roughness):
        # Top wall pillar generation
        # Roughness surface features preparation
        # Determines location of pillars
        # If there is arbitrary physical roughness
        if (arbitraryPhysicalRoughness):
            # Set x and y equal to their initial values
            xp = int(x0/latticeParameter)
            yp = int(y0/latticeParameter)
            # While x is less than the number of x wall atoms
            while (xp < int((x0 + wall_x)/latticeParameter)):
                # Increment x by a random number of atoms within the given range
                xp = xp + randrange(spacingMin, (spacingMax + 1))
                # While y is less than the number of y wall atoms
                while (yp < int((y0 + wall_y)/latticeParameter)):
                    # Assigns resname and atomname to be equal to 'WAL' to
                    # represent wall atoms. Default wall atoms are provided
                    # and remain unchanged unless another is selected
                    resname = 'WAL'
                    atomname = 'WAL'
                    # Increment y by a random spacing within the given range
                    yp = yp + randrange(spacingMin, (spacingMax + 1))
                    # If there is to be arbitrary chemical roughness
                    if (arbitraryChemicalRoughness | chemicalRoughness):
                        # Obtains a random atom type for the pillar
                        atomType, discard1, discard2 = chemically_rough_line(resname, atomname, pattern)
                        # Increments position by one so the next pillar will
                        # have a different atom type
                        pattern = pattern + 1
                        # Add the pillar to the list of all pillars
                        pillarsTop.append([xp, yp, atomType])
                    # Otherwise add a default wall atom to the wall
                    else:
                        # Add the pillar x and y values to a list of all pillars
                        pillarsTop.append([xp, yp, [resname, atomname]])
                # Resets y to its initial value for the next iteration in x
                yp = int(y0/latticeParameter)
        # If there is not arbitrary physical roughness
        else:
            # For every x location within the walls
            for pillarX in range(int(x0/latticeParameter), int((x0 + wall_x)/latticeParameter)):
                # For every y location within the walls at that x location
                for pillarY in range(int(y0/latticeParameter), int((y0 + wall_y)/latticeParameter)):
                    # Assigns resname and atomname to be equal to 'WAL' to
                    # represent wall atoms. Default wall atoms are provided
                    # and remain unchanged unless another is selected
                    resname = 'WAL'
                    atomname = 'WAL'
                    # If the location is spacing atoms from the initial walls or
                    # previous pillar in both x and y 
                    if ((pillarX%spacing == 0) and (pillarY%spacing == 0)):
                        # If there is to be chemical roughness of any type
                        if (arbitraryChemicalRoughness or chemicalRoughness):
                            # Obtains a random atom type for the pillar
                            atomType, discard1, discard2 = chemically_rough_line(resname, atomname, pattern)
                            # Increments position by one so the next pillar will
                            # have a different atom type
                            pattern = pattern + 1
                            # Add the pillar to the list of all pillars
                            pillarsTop.append([pillarX, pillarY, atomType])
                        # Otherwise add a default wall atom to the wall
                        else:
                            # Add the pillar to the list of all pillars
                            pillarsTop.append([pillarX, pillarY, [resname, atomname]])
        # Bottom wall pillar generation
        # Roughness surface features preparation
        # Determines location of pillars
        # If there is arbitrary physical roughness
        if (arbitraryPhysicalRoughness):
            # Set x and y equal to their initial values
            xp = int(x0/latticeParameter)
            yp = int(y0/latticeParameter)
            # While x is less than the number of x wall atoms
            while (xp < int((x0 + wall_x)/latticeParameter)):
                # Increment x by a random number of atoms within the given range
                xp = xp + randrange(spacingMin, (spacingMax + 1))
                # While y is less than the number of y wall atoms
                while (yp < int((y0 + wall_y)/latticeParameter)):
                    # Assigns resname and atomname to be equal to 'WAL' to
                    # represent wall atoms. Default wall atoms are provided
                    # and remain unchanged unless another is selected
                    resname = 'WAL'
                    atomname = 'WAL'
                    # Increment y by a random spacing within the given range
                    yp = yp + randrange(spacingMin, (spacingMax + 1))
                    # If there is to be arbitrary chemical roughness
                    if (arbitraryChemicalRoughness or chemicalRoughness):
                        # Obtains a random atom type for the pillar
                        atomType, discard1, discard2 = chemically_rough_line(resname, atomname, pattern)
                        # Increments position by one so the next pillar will
                        # have a different atom type
                        pattern = pattern + 1
                        # Add the pillar to the list of all pillars
                        pillarsBottom.append([xp, yp, atomType])
                    # Otherwise add a default wall atom to the wall
                    else:
                        # Add the pillar x and y values to a list of all pillars
                        pillarsBottom.append([xp, yp, [resname, atomname]])
                # Resets y to its initial value for the next iteration in x
                yp = int(y0/latticeParameter)
        # If there is not arbitrary physical roughness
        else:
            # For every x location within the walls
            for pillarX in range(int(x0/latticeParameter), int((x0 + wall_x)/latticeParameter)):
                # For every y location within the walls at that x location
                for pillarY in range(int(y0/latticeParameter), int((y0 + wall_y)/latticeParameter)):
                    # Assigns resname and atomname to be equal to 'WAL' to
                    # represent wall atoms. Default wall atoms are provided
                    # and remain unchanged unless another is selected
                    resname = 'WAL'
                    atomname = 'WAL'
                    # If the location is spacing atoms from the initial walls or
                    # previous pillar in both x and y 
                    if ((pillarX%spacing == 0) and (pillarY%spacing == 0)):
                        # If there is to be chemical roughness of any type
                        if (arbitraryChemicalRoughness or chemicalRoughness):
                            # Obtains a random atom type for the pillar
                            atomType, discard1, discard2 = chemically_rough_line(resname, atomname, pattern)
                            # Increments position by one so the next pillar will
                            # have a different atom type
                            pattern = pattern + 1
                            # Add the pillar to the list of all pillars
                            pillarsBottom.append([pillarX, pillarY, atomType])
                        # Otherwise add a default wall atom to the wall
                        else:
                            # Add the pillar to the list of all pillars
                            pillarsBottom.append([pillarX, pillarY, [resname, atomname]])
    
    # If there is any type of roughness, generate the base wall layers
    if(roughness):
        # Empty array so the code continues even without pillars
        pil = []
        # Obtains the residue name from the atom type
        resname = 'WAL'
        # Obtains the atom name from the atom type
        atomname = 'WAL'
        # Current layer in the wall base layers
        baseLayerHeight = 0
        # While the current layer is less then the base layer height
        while (baseLayerHeight < baseLayers):
            # Determine the appropriate z values for the wall layers
            baseLayerZBottom = z0 - ((baseLayerHeight + height)*latticeParameter)
            baseLayerZTop = z0 + ((baseLayerHeight + height)*latticeParameter)
            # Obtains the lines required to represent a wall and counter
            # information
            wallLines, counterReset = generate_wall(x0, y0, baseLayerZBottom, natoms, 'bottom', pil, True, levelBottom, radius, resname, atomname)
            # Add the wall atoms to the system for that layer
            lines.extend(wallLines)
            # Increase natoms to reflect the atoms that have been added
            natoms = natoms + len(wallLines)
            # If there are more than 100000 atoms, the counter was reset
            if (counterReset != 0):
                # Adjusts counter to reflect correct amount of atoms
                natoms = natoms - (counterReset*100000)
            # Obtains the lines required to represent a wall and counter
            # information
            wallLines, counterReset = generate_wall(x0, y0, baseLayerZTop, natoms, 'top', pil, True, levelBottom, radius, resname, atomname)
            # Add the wall atoms to the system for that layer
            lines.extend(wallLines)
            # Increase natoms to reflect the atoms that have been added
            natoms = natoms + len(wallLines)
            # If there are more than 100000 atoms, the counter was reset
            if (counterReset != 0):
                # Adjusts counter to reflect correct amount of atoms
                natoms = natoms - (counterReset*100000)
            # Increments the layer z value for the next iteration
            baseLayerHeight = baseLayerHeight + 1

        # The bottom wall z value must be decreased to make room for roughness
        z01 = z0 - (height*latticeParameter)
        # The upper wall z value must be increased to fit the roughness
        z02 = z0 + (height*latticeParameter)
        # Pillars must be given as a list but must also be analyzed individually
        # for arbitrary physical roughness so this list is generated to hold a
        # pillar and is cleared every time another pillar is to be analyzed
        pil = list()
        # If there is arbitrary physical roughness in the system
        if (arbitraryPhysicalRoughness):
            # For each of the individual pillars arbitrarily positioned
            for pillar in pillarsBottom:
                # Adds the current pillar to the list intended to hold it
                pil.append(pillar)
                # Determines a random radius for the current pillar that is
                # within the range specified
                radius = randrange(radiusMin, (radiusMax + 1))
                # Obtains the atom type from the pillar
                aType = pillar[2]
                # Obtains the residue name from the atom type
                resname = aType[0]
                # Obtains the atom name from the atom type
                atomname = aType[1]
                # Resets the wall layer to the outermost wall level
                levelBottom = 1
                # Layers of wall must be created to obtain the desired wall
                # thickness. The number of layers for each pillar is generated
                # randomly in the desired range of atoms provided
                # Wall height begins at one atom above the outermost wall layer
                r = 1
                # Randomly sets the height for the pillar in the provided range
                limit = randrange(heightMin, (height + 1))
                # While the current z value is less then the pillar height
                while (r <= limit):
                    # Determine the appropriate z value for the wall layer
                    z = z01 + (r*latticeParameter)
                    # Obtains the lines required to represent a wall and counter
                    # information
                    wallLines, counterReset = generate_wall(x0, y0, z, natoms, 'bottom', pil, False, levelBottom, radius, resname, atomname)
                    # After the outer wall layer bottomWallOutside becomes false
                    bottomWallOutside = False
                    # Increment level as another layer has been added
                    levelBottom = levelBottom + 1
                    # Add the wall atoms to the system for that layer
                    lines.extend(wallLines)
                    # Increase natoms to reflect the atoms that have been added
                    natoms = natoms + len(wallLines)
                    # If there are more than 100000 atoms, the counter was reset
                    if (counterReset != 0):
                        # Adjusts counter to reflect correct amount of atoms
                        natoms = natoms - (counterReset*100000)
                    # Increments the layer z value for the next iteration
                    r = r + 1
                # Clears the list holding a single pillar for the next iteration
                # While there is at least one element in the pil list
                while (len(pil) > 0):
                    # Remove the last element in the list
                    pil.pop()
            # If there are multiple walls in the system
            if (nwall != 1):
                # For each of the individual pillars arbitrarily positioned
                for pillar in pillarsTop:
                    # Adds the current pillar to the list intended to hold it
                    pil.append(pillar)
                    # Determines a random radius for the current pillar that is
                    # within the range specified
                    radius = randrange(radiusMin, (radiusMax + 1))
                    # Obtains the atom type from the pillar
                    aType = pillar[2]
                    # Obtains the residue name from the atom type
                    resname = aType[0]
                    # Obtains the atom name from the atom type
                    atomname = aType[1]
                    # Resets the wall layer from the outermost wall level
                    levelTop = 1
                    # Generate the top wall
                    # Wall Layers must be created to obtain the desired wall
                    # thickness. The number of layers for each pillar is
                    # generated randomly in the desired range of atoms provided
                    # Wall height begins at one atom above the outer wall layer
                    r = 1
                    # Randomly sets the pillar height within the provided range
                    limit = randrange(heightMin, (height + 1))
                    # While the current z value is less then the pillar height
                    while (r <= limit):
                        # Determine the appropriate z value for the wall layer
                        z = z02 - (r*latticeParameter)
                        # Obtains the lines to represent a wall and counter
                        # information
                        wallLines, counterReset = generate_wall(x0, y0, z, natoms, 'top', pil, False, levelTop, radius, resname, atomname)
                        # After the outer wall layer topWallOutside is false
                        topWallOutside = False
                        # Increment level as another layer has been added
                        levelTop = levelTop + 1
                        # Add the wall atoms to the system for that layer
                        lines.extend(wallLines)
                        # Increase natoms to reflect the added atoms
                        natoms = natoms + len(wallLines)
                        # If there are more than 100000 atoms, counter is reset
                        if (counterReset != 0):
                            # Adjusts counter to reflect correct amount of atoms
                            natoms = natoms - (counterReset*100000)
                        # Increments the layer z value for the next iteration
                        r = r + 1
                    # Clears the list holding one pillar for the next iteration
                    # While there is at least one element in the pil list
                    while (len(pil) > 0):
                        # Remove the last element in the list
                        pil.pop()
                            
        # If there is not arbitrary physical roughness in the system
        else:
            # For every pillar in the system
            for pillar in pillarsBottom:
                # Adds the current pillar to the list intended to hold it
                pil.append(pillar)
                # Obtains the atom type from the pillar
                aType = pillar[2]
                # Obtains the residue name from the atom type
                resname = aType[0]
                # Obtains the atom name from the atom type
                atomname = aType[1]
                # Resets the wall layer from the outermost wall level
                levelBottom = 1
                # Sets the initial wall level to one atom above the outer layer
                r = 1
                # While the level is less than the thickness to be generated
                while (r < (height + 1)):
                    # Determine the appropriate z value for the wall layer
                    z = z01 + (r*latticeParameter)
                    # Obtains the lines required to represent a wall and counter
                    # information
                    wallLines, counterReset = generate_wall(x0, y0, z, natoms, 'bottom', pil, False, levelBottom, radius, resname, atomname)
                    # After the outer wall layer bottomWallOutside becomes false
                    bottomWallOutside = False
                    # Increment level as another layer has been added
                    levelBottom = levelBottom + 1
                    # Add the wall atoms to the system for that layer
                    lines.extend(wallLines)
                    # Increase natoms to reflect the atoms that have been added
                    natoms = natoms + len(wallLines)
                    # If there are more than 100000 atoms, the counter was reset
                    if (counterReset != 0):
                        # Adjusts counter to reflect correct amount of atoms
                        natoms = natoms - (counterReset*100000)
                    # Increment the wall layer z value for the next iteration
                    r = r + 1
                    # Clears the list holding one pillar for the next iteration
                    # While there is at least one element in the pil list
                while (len(pil) > 0):
                    # Remove the last element in the list
                    pil.pop()
            # If there are multiple walls in the system
            if (nwall != 1):
                # Generate the top wall
                # For every pillar in the top wall
                for pillar in pillarsTop:
                    # Adds the current pillar to the list intended to hold it
                    pil.append(pillar)
                    # Obtains the atom type from the pillar
                    aType = pillar[2]
                    # Obtains the residue name from the atom type
                    resname = aType[0]
                    # Obtains the atom name from the atom type
                    atomname = aType[1]
                    # Resets the wall layer from the outermost wall level
                    levelTop = 1
                    # Sets the initial level to one atom above the outer layer
                    r = 1
                    # While the level is less than the thickness to be generated
                    while (r < (height + 1)):
                        # Determine the appropriate z value for the wall layer
                        z = z02 - (r*latticeParameter)
                        # Gets the lines to represent a wall and counter
                        # information
                        wallLines, counterReset = generate_wall(x0, y0, z, natoms, 'top', pil, False, levelTop, radius, resname, atomname)
                        # After the outer wall layer topWallOutside is false
                        topWallOutside = False
                        # Increment level as another layer has been added
                        levelTop = levelTop + 1
                        # Add the wall atoms to the system for that layer
                        lines.extend(wallLines)
                        # Increase natoms to reflect the added atoms
                        natoms = natoms + len(wallLines)
                        # If there are more than 100000 atoms, counter is reset
                        if (counterReset != 0):
                            # Adjusts counter to reflect correct amount of atoms
                            natoms = natoms - (counterReset*100000)
                        # Increment the layer z value for the next iteration
                        r = r + 1
                    # Clears the list holding one pillar for the next iteration
                    # While there is at least one element in the pil list
                    while (len(pil) > 0):
                        # Remove the last element in the list
                        pil.pop()
                        
    # If there is to be no roughness in the wall system                        
    else:
        # Generates the smooth bottom wall of the system
        # Sets z equal to the initial value to generate the bottom wall
        z = z0
        # Assigns resname and atomname to be equal to 'WAL' to represent wall
        # atoms. Default wall atoms are provided and remain unchanged unless
        # another is selected
        resname = 'WAL'
        atomname = 'WAL'
        # Obtains lines required to represent a wall and counter information
        wallLines, counterReset = generate_wall(x0, y0, z, natoms, 'bottom', pillarsBottom, True, levelBottom, radius, resname, atomname)
        # Add the wall atoms to the system for that layer
        lines.extend(wallLines)
        # Increase natoms to reflect the atoms that have been added
        natoms = natoms + len(wallLines)
        # If there are more than 100000 atoms, the counter is reset
        if (counterReset != 0):
            # Adjusts counter to reflect correct amount of atoms
            natoms = natoms - (counterReset*100000)
        # If there is not one wall, then the second smooth wall is generated
        if (nwall != 1):
            # Adds the top wall to lines if the number of walls is not one
            # Obtains the lines required to represent a wall and counter
            # information
            wallLines, counterReset = generate_wall(x0, y0, z, natoms, 'top', pillarsTop, True, levelTop, radius, resname, atomname)
            # Add the wall atoms to the system for that layer
            lines.extend(wallLines)
            # Increase natoms to reflect the atoms that have been added
            natoms = natoms + len(wallLines)
            # If there are more than 100000 atoms, the counter is reset
            if (counterReset != 0):
                # Adjusts counter to reflect correct amount of atoms
                natoms = natoms - (counterReset*100000)

    # Returns a list containing all the lines required to represent the wall(s)
    return lines


def generate_fluid(z, atomnum, znum):
    """
    Used by generate_fluid_gro() to generate the fluid atoms in the x and y
    directions at a particular z value
    """
    
    # Creates a list used to return the lines containing the fluid atoms
    lines = list()
    # Counter for atomnum to be reset with if it goes over five digits
    counterReset = 0
    # Y center of the box 
    yc = (box_y/2)
    # Assigns resname and atomname to be equal to 'FLD' to represent fluid
    resname = 'FLD'
    atomname = 'FLD'
    
    # Sets x to its initial value
    x = (box_x - wall_x)/2
    # Accumulator to track the number of atoms in the x direction
    i = 0
    # While loop iterates through all desired x values
    while (x < ((wall_x + (box_x - wall_x)/2) - (0.5*factor*sigma_wall))):
        # Sets y equal to its initial value for each x value
        y = (yc - ((fluidFillFraction_y/2)*wall_y))
        # Accumulator to track the number of atoms in the y direction
        j = 0
        # While loop iterates through all desired y values at each x value
        while (y < (yc + ((fluidFillFraction_y/2)*wall_y))):
            # Assigns resnum equal to i plus j
            resnum = (i + j)
            # Increments the atom number as one more atom of fluid
            # has been added to the system
            atomnum = (atomnum + 1)
            # While atomnum is greater than 5 digits it is reset
            while (atomnum > 99999):
                # Resets atomnum
                atomnum = (atomnum - 100000)
                # Increments counterReset
                counterReset = (counterReset + 1)
            # Increments j as one more atom in the y direction is added
            j = (j + 1)

            # Every other level of fluid is staggered in the x and y directions
            # If the level is to be staggered
            if (znum%2 == 0):
                # Increments x and y by 0.5 sigma
                y = (y + (0.5*factor*sigma_fluid))
                x = (x + (0.5*factor*sigma_fluid))
                # Generate a line for the atom
                line = generate_line(resnum, resname, atomname, atomnum, x, y, z)
                # Returns x and y to their not staggered values
                y = (y - (0.5*factor*sigma_fluid))
                x = (x - (0.5*factor*sigma_fluid))
            # If the line is one that will not be staggered
            else:
                # Writes the line that will be used to represent the new atom
                line = generate_line(resnum, resname, atomname, atomnum, x, y, z)

            # Adds the line to the list lines that will be returned
            lines.append(line)
            # Increments the value of y by sigma times the spacing factor
            y = (y + (sigma_fluid*factor))
        # Increments the value of x by sigma times the spacing factor
        x = (x + (sigma_fluid*factor))
        # Increments i by one as one more atom is added in the x direction
        i = (i + 1)
        
    # Returns a list containing the lines for all fluid atoms at a given z value
    return lines, counterReset


def generate_fluid_gro():
    """
    Generates the fluid within the system.gro file by using generate_fluid()
    to generate fluid at all x and y values for every desired z value
    """

    # Sets the initial z value to the z value immediately above the bottom wall
    z = ((factor*latticeParameter) + ((box_z - (nwall*wall_z) - dz)/2))
    # Tracks the number of atoms of fluid that have been added to the system
    atomnum = 0
    # Tracks the number of fluid atom levels added in the z direction
    znum = 1
    # Creates a list used to store a line representing every atom of fluid
    lines = list()
    
    # While loop to iterate through the space in the z direction
    while (z < ((dz + (box_z - dz - (nwall*wall_z))/2) - (sigma_wall*factor))):
        # Generates fluid at all of the desired x and y values
        line, counterReset = generate_fluid(z, atomnum, znum)
        # Adds the fluid at each z value to the list containing all lines
        lines.extend(line)
        # Increments the atom number by the number of atoms added
        atomnum = atomnum + len(line)
        # If there are more than 100000 atoms, the counter is reset
        if (counterReset != 0):
            # Adjusts counter to reflect correct amount of atoms
            atomnum = (atomnum - (counterReset*100000))
        # Increments z till the value immediately before the top wall is reached
        # Incrementing corresponds to a 60 degree angle on a 30-60-90 triangle
        z = z + (factor*sigma_fluid*((3**0.5)/2))
        # Increment znum by 1 as an additional level in the z direction is added
        znum = znum + 1
        
    # Returns the list lines containing the lines representing the fluid
    return lines


def generate_wall_top(gro_lines, nWall, nFluid, sigma_wall, wall_epsilon, wall_mass, sigma_fluid, fluid_epsilon, fluid_mass, roughAtoms):
    """
    Generates the wall topology file lines
    """

    # Obtains the additionalTopLines and atomTypes to be added to the topology
    # file
    discard, additionalTopLines, atomTypes = chemically_rough_line('Anything', 'Anything', 0 )

    # Charge for the default wall and fluid atoms
    charge = 0

    # List to store information about all other atom types besides WAL and FLD
    othAt = list()

    # For every non default atom type provided
    for atom in roughAtoms:
        # Add an array containing that atom and the number of those atoms in the
        # system (initially 0) to a list
        othAt.append([str(atom[1]), 0])
    # For each non default wall atom type
    for atom in othAt:
        # For every line that represents an atom in the system
        for line in gro_lines:
            # If the line represents that atom type
            if str(atom[0]) in str(line):
                # Increment the element of the array to represent an additional
                # atom of that type being present within the system
                atom[1] = atom[1] + 1
                # Decrement nWall as there is one fewer default wall atoms
                nWall = nWall - 1

    # Begins a string that will be used to represent non default wall atoms
    otherAtoms = ''
    # For every non default wall atom type
    for atom in othAt:
        # Obtain the atom name and the number of those atoms in the system
        atNam = str(atom[0])
        atNum = str(atom[1])
        # Add a line representing that atom type to otherAtoms
        otherAtoms = otherAtoms + '%(atNam)s%(atNum)17s\n'%vars() 

    # Generates the correctly formatted lines for the topology file
    lines = """\
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
  1             2               no              1.0     1.0

[ atomtypes ]
WAL         35      %(wall_mass)1.2s    0.0000  A   %(sigma_wall)1.5e  %(wall_epsilon)1.5e
FLD         35      %(fluid_mass)1.2s    0.0000  A   %(sigma_fluid)1.5e  %(fluid_epsilon)1.5e
%(additionalTopLines)s

[ moleculetype ]
; name  nrexcl
WAL          1

[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass
     1   WAL    1   WAL    WAL    1%(charge)13.6f%(wall_mass)13.5f
%(atomTypes)s
[ moleculetype ]
; name  nrexcl
FLD          1

[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass
     1   FLD    1   FLD    FLD    1%(charge)13.6f%(fluid_mass)13.5f

[ system ]
LJ_Fluid

[ molecules ]
WAL%(nWall)18s
%(otherAtoms)sFLD%(nFluid)18s
""" %vars()

    # Returns the list of lines required to generate the file and information
    # about the names and number of each non default wall atom
    return lines, othAt


#===============================================================================
# Main
#===============================================================================

# Other Constants. (DO NOT CHANGE!)
# Number of walls in the system. Note: should always equal two
nwall = 2
# The number of atoms of thickness for a flat wall. Note: should always be two
wall_z = 1*sigma_wall
# Box
# X, Y, and Z dimensions of the system box
# Note: box_x should equal wall_x and box_y should equal wall_y
# Note: box_z should exceed dz + nwall*wall_z
box_x = wall_x
box_y = wall_y
box_z = 2*((dz + (nwall*wall_z)) + 2*latticeParameter)
# Names of files that will be required in order to create the desired system
# (DO NOT CHANGE!)
minimize_mdp_file = 'minimize.mdp'
nvt_mdp_file = 'nvt.mdp'
gro_file = 'system.gro'
top_file = 'system.top'
ndx_file = 'index.ndx'
# Creates paths that will be used to generate various directories
work_path = base_path + 'simulation/'

# Creates the 'simulation' directory by making a new directory
os.mkdir(work_path)

# Changes directory to the work path
os.chdir(work_path)

# Generate Wall 
wall_lines = generate_wall_gro(radius)
nWall = len(wall_lines)

# Generate Fluid
fluid_lines = generate_fluid_gro()
nFluid = len(fluid_lines)
natoms = nWall + nFluid

#===============================================================================
# Files
#===============================================================================

# Arranges lines of .gro file in correct order according to the topology file
atomOrder = list()
atomOrder.append('WAL')
for atom in roughAtoms:
    atomOrder.append(str(atom[0]))
wall_lines_final = list()
for atomType in atomOrder:
    for line in wall_lines:
        if (str(line[5:10]).strip() == str(atomType)):
            wall_lines_final.append(line)

# Write to gro file
gro_lines = list()
gro_lines.append('system\n')
gro_lines.append('%(natoms)5d\n'%vars())
gro_lines = gro_lines + wall_lines_final + fluid_lines
gro_lines.append('%(box_x)10.5f%(box_y)10.5f%(box_z)10.5f\n'%vars())
write_file(os.path.join(work_path, gro_file), gro_lines)
          
# Write top file
top_lines, otherAtoms = generate_wall_top(gro_lines, nWall, nFluid, sigma_wall, wall_epsilon, wall_mass, sigma_fluid, fluid_epsilon, fluid_mass, roughAtoms)
write_file(os.path.join(work_path, top_file), top_lines)

# List to store all non default wall atoms that are present in the system
atomTypes = list()
# For every non default wall atom type
for atom in otherAtoms:
    # If there is at least one atom of this type in the system
    if (atom[1] != 0):
        # Add the atom type name to the list of non default wall atoms
        atomTypes.append(atom)
# Begins a string to represent all energy groups
energyGroups = 'WAL'
# Begins a string to represent frozen dimensions of energy groups
freezedim = 'Y Y Y'
# Begins a string to represent energy group bond exclusions
energyExcl = 'WAL WAL  '
# Accumulator for total non default atoms
totalAt = 0
# For every additional atom type with at least one occurrence
for atom in atomTypes:
    # Determine the atom name
    atoName = atom[0]
    # Add the information required to represent that atom in the lines for
    # energy groups, frozen dimensions, and energy group bond exclusions
    energyGroups = energyGroups + '  %(atoName)s'%vars()
    freezedim = freezedim + '  Y Y Y'
    energyExcl = energyExcl + 'WAL %(atoName)s  %(atoName)s %(atoName)s  '%vars()
    # Add number of atoms to the total number of atoms
    totalAt = totalAt + atom[1]
    # Note: Only the major and most computationally expensive wall bonds are
    # excluded, this reduces the time required to run the simulation as all
    # types of wall atoms are frozen anyways to any type of WAL-WAL interaction
    # can be neglected

# If dynamic fluid motion is specified generate the correct .mdp files
if (dynamic):
    # Generates the minimize.mdp file lines
    minmdp_lines ="""
; RUN CONTROL PARAMETERS =
integrator               = l-bfgs
tinit                    = 0
dt                       = 0.01
nsteps                   = 1000
nstcomm                  = 1
nsstcouple               = 1

; NEIGHBORSEARCHING PARAMETERS =
; nblist update frequency=
nstlist                  = 10
; ns algorithm (simple or grid) =
ns_type                  = grid
; Periodic boundary conditions: xyz or none =
pbc                      = xyz
; nblist cut-off         =
rlist                    = 1.1

; OPTIONS FOR ELECTROSTATICS AND VDW =
; Method for doing electrostatics =
coulombtype              = cutoff ; can be changed to cutoff -- there are no charges.  Should already turn off PME, but you can tell by looking at the logs if it does or not.
;rcoulomb-switch         = 0
rcoulomb                 = 1.1
; Dielectric constant (DC) for cut-off or DC of reaction field =
epsilon-r                = 1
; Method for doing Van der Waals =
vdw-type                 = switch
; cut-off lengths        =
rvdw-switch              = 0.9135 ;3 sigma
rvdw                     = 1.06575 ;3.5 sigma
; Apply long range dispersion corrections for Energy and Pressure =
DispCorr                 = EnerPres  
; Spacing for the PME/PPPM FFT grid =
fourierspacing           = 0.1
; FFT grid size, when a value is 0 fourierspacing will be used =
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters =
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no

; OUTPUT CONTROL OPTIONS =
; Output frequency for coords (x), velocities (v) and forces (f) =
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file =
nstlog                   = 500
nstenergy                = 100
; Output frequency and precision for xtc file =
nstxtcout                = 1000
xtc-precision            = 1000

; ENERGY MINIMIZATION OPTIONS = 
; Force tolerance and initial step-size = 
emtol                    = 10
emstep                   = 1.0e-4
; Max number of iterations in relax_shells = 
niter                    = 20
; Step size (1/ps^2) for minimization of flexible constraints = 
fcstep                   = 0
; Frequency of steepest descents steps when doing CG = 
nstcgsteep               = 1000

; Non-equilibrium MD stuff
acc-grps                 =
accelerate               =
freezegrps               = %(energyGroups)s  
freezedim                = %(freezedim)s 
energygrps		         = %(energyGroups)s
energygrp-excl		     = %(energyExcl)s
cos-acceleration         = 0
deform                   =
"""%vars()
    # Writes the minimize.mdp file
    write_file(os.path.join(work_path, 'minimize.mdp'), minmdp_lines)

    # Generates the nvt.mdp file lines
    nvtmdp_lines ="""
; RUN CONTROL PARAMETERS =
integrator               = sd
tinit                    = 0  ; start time in ps
dt                       = 0.01 ; timestep
nsteps                   = 100000
nstcomm                  = 0 ; com removal steps
nsstcouple               = 1

; NEIGHBORSEARCHING PARAMETERS =
; nblist update frequency=
nstlist                  = 10
; ns algorithm (simple or grid) =
ns_type                  = grid
; Periodic boundary conditions: xyz or none =
pbc                      = xyz
; nblist cut-off         =
rlist                    = 1.1

; OPTIONS FOR ELECTROSTATICS AND VDW =
; Method for doing electrostatics =
coulombtype              = cutoff
;rcoulomb-switch         = 0
rcoulomb                 = 1.1
; Dielectric constant (DC) for cut-off or DC of reaction field =
epsilon-r                = 1
; Method for doing Van der Waals =
vdw-type                 = switch
; cut-off lengths        =
rvdw-switch              = 0.9135 ;3 sigma
rvdw                     = 1.06575 ;3.5 sigma
; Spacing for the PME/PPPM FFT grid =
fourierspacing           = 0.1
; FFT grid size, when a value is 0 fourierspacing will be used =
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters =
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no

; OUTPUT CONTROL OPTIONS =
; Output frequency for coords (x), velocities (v) and forces (f) =
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file =
nstlog                   = 500
nstenergy                = 100
; Output frequency and precision for xtc file =
nstxtcout                = 20
xtc-precision            = 1000

;OPTIONS FOR TEMPERATURE COUPLING
tc_grps                  = system
tau_t                    = 10
ref_t                    = %(temperature)s
; Non-equilibrium MD 
acc-grps                 = FLD
accelerate               = 0 %(acceleration)s 0
freezegrps               = %(energyGroups)s  
freezedim                = %(freezedim)s 
energygrps               = %(energyGroups)s
energygrp-excl	         = %(energyExcl)s
cos-acceleration         = 0
deform                   =
"""%vars()
    # Writes the nvt.mdp file
    write_file(os.path.join(work_path, 'nvt.mdp'), nvtmdp_lines)

    # Generates the simulation.mdp file lines
    simmdp_lines ="""
; RUN CONTROL PARAMETERS =
integrator               = sd
tinit                    = 0 ; start time in ps
dt                       = 0.01 ; timestep
nsteps                   = 500000
nstcomm                  = 0 ; com removal steps
nsstcouple               = 1

; NEIGHBORSEARCHING PARAMETERS =
; nblist update frequency=
nstlist                  = 10
; ns algorithm (simple or grid) =
ns_type                  = grid
; Periodic boundary conditions: xyz or none =
pbc                      = xyz
; nblist cut-off         =
rlist                    = 1.1

; OPTIONS FOR ELECTROSTATICS AND VDW =
; Method for doing electrostatics =
coulombtype              = cutoff
;rcoulomb-switch         = 0
rcoulomb                 = 1.1
; Dielectric constant (DC) for cut-off or DC of reaction field =
epsilon-r                = 1
; Method for doing Van der Waals =
vdw-type                 = switch
; cut-off lengths        =
rvdw-switch              = 0.9135 ;3 sigma
rvdw                     = 1.06575 ;3.5 sigma
; Spacing for the PME/PPPM FFT grid =
fourierspacing           = 0.1
; FFT grid size, when a value is 0 fourierspacing will be used =
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters =
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no

; OUTPUT CONTROL OPTIONS =
; Output frequency for coords (x), velocities (v) and forces (f) =
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file =
nstlog                   = 500
nstenergy                = 100
; Output frequency and precision for xtc file =
nstxtcout                = 10
xtc-precision            = 1000

;OPTIONS FOR TEMPERATURE COUPLING
tc_grps                  = system
tau_t                    = 10
ref_t                    = %(temperature)s
; Non-equilibrium MD 
acc-grps                 = FLD
accelerate               = 0 %(acceleration)s 0
freezegrps               = %(energyGroups)s  
freezedim                = %(freezedim)s 
energygrps	         = %(energyGroups)s
energygrp-excl	         = %(energyExcl)s
cos-acceleration         = 0
deform                   =
"""%vars()
    # Writes the simulation.mdp file
    write_file(os.path.join(work_path, 'simulation.mdp'), simmdp_lines)

# Otherwise the system is static and the correct .mdp files are generated
else:
    # Generates the minimize.mdp file lines
    minmdp_lines ="""
; RUN CONTROL PARAMETERS =
integrator               = l-bfgs
tinit                    = 0
dt                       = 0.01
nsteps                   = 1000
nstcomm                  = 1
nsstcouple               = 1

; NEIGHBORSEARCHING PARAMETERS =
; nblist update frequency=
nstlist                  = 10
; ns algorithm (simple or grid) =
ns_type                  = grid
; Periodic boundary conditions: xyz or none =
pbc                      = xyz
; nblist cut-off         =
rlist                    = 1.1

; OPTIONS FOR ELECTROSTATICS AND VDW =
; Method for doing electrostatics =
coulombtype              = cutoff ; can be changed to cutoff -- there are no charges.  Should already turn off PME, but you can tell by looking at the logs if it does or not.
;rcoulomb-switch         = 0
rcoulomb                 = 1.1
; Dielectric constant (DC) for cut-off or DC of reaction field =
epsilon-r                = 1
; Method for doing Van der Waals =
vdw-type                 = switch
; cut-off lengths        =
rvdw-switch              = 0.9135 ;3 sigma
rvdw                     = 1.06575 ;3.5 sigma
; Apply long range dispersion corrections for Energy and Pressure =
DispCorr                 = EnerPres  
; Spacing for the PME/PPPM FFT grid =
fourierspacing           = 0.1
; FFT grid size, when a value is 0 fourierspacing will be used =
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters =
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no

; OUTPUT CONTROL OPTIONS =
; Output frequency for coords (x), velocities (v) and forces (f) =
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file =
nstlog                   = 500
nstenergy                = 100
; Output frequency and precision for xtc file =
nstxtcout                = 1000
xtc-precision            = 1000

; ENERGY MINIMIZATION OPTIONS = 
; Force tolerance and initial step-size = 
emtol                    = 10
emstep                   = 1.0e-4
; Max number of iterations in relax_shells = 
niter                    = 20
; Step size (1/ps^2) for minimization of flexible constraints = 
fcstep                   = 0
; Frequency of steepest descents steps when doing CG = 
nstcgsteep               = 1000

; Non-equilibrium MD stuff
acc-grps                 =
accelerate               =
freezegrps               = %(energyGroups)s  
freezedim                = %(freezedim)s 
energygrps	         = %(energyGroups)s
energygrp-excl	         = %(energyExcl)s
cos-acceleration         = 0
deform                   =
"""%vars()
    # Writes the minimize.mdp file
    write_file(os.path.join(work_path, 'minimize.mdp'), minmdp_lines)

    # Generates the nvt.mdp file lines
    nvtmdp_lines ="""
; RUN CONTROL PARAMETERS =
integrator               = md
tinit                    = 0  ; start time in ps
dt                       = 0.01 ; timestep
nsteps                   = 100000
nstcomm                  = 1 ; com removal steps
nsstcouple               = 1

; NEIGHBORSEARCHING PARAMETERS =
; nblist update frequency=
nstlist                  = 10
; ns algorithm (simple or grid) =
ns_type                  = grid
; Periodic boundary conditions: xyz or none =
pbc                      = xyz
; nblist cut-off         =
rlist                    = 1.1

; OPTIONS FOR ELECTROSTATICS AND VDW =
; Method for doing electrostatics =
coulombtype              = cutoff
;rcoulomb-switch         = 0
rcoulomb                 = 1.1
; Dielectric constant (DC) for cut-off or DC of reaction field =
epsilon-r                = 1
; Method for doing Van der Waals =
vdw-type                 = switch
; cut-off lengths        =
rvdw-switch              = 0.9135 ;3 sigma
rvdw                     = 1.06575 ;3.5 sigma
; Spacing for the PME/PPPM FFT grid =
fourierspacing           = 0.1
; FFT grid size, when a value is 0 fourierspacing will be used =
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters =
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no

; OUTPUT CONTROL OPTIONS =
; Output frequency for coords (x), velocities (v) and forces (f) =
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file =
nstlog                   = 500
nstenergy                = 100
; Output frequency and precision for xtc file =
nstxtcout                = 20
xtc-precision            = 1000

;OPTIONS FOR TEMPERATURE COUPLING
tc_grps                  = system
tcoupl                   = v-rescale
tau_t                    = 0.1 
ref_t                    = %(temperature)s
; Non-equilibrium MD 
acc-grps                 =
accelerate               =
freezegrps               = %(energyGroups)s  
freezedim                = %(freezedim)s 
energygrps	         = %(energyGroups)s
energygrp-excl	         = %(energyExcl)s
cos-acceleration         = 0
deform                   =
"""%vars()
    # Writes the nvt.mdp file
    write_file(os.path.join(work_path, 'nvt.mdp'), nvtmdp_lines)

    # Generates the simulation.mdp file lines
    simmdp_lines ="""
; RUN CONTROL PARAMETERS =
integrator               = md
tinit                    = 0 ; start time in ps
dt                       = 0.01 ; timestep
nsteps                   = 500000
nstcomm                  = 1 ; com removal steps
nsstcouple               = 1

; NEIGHBORSEARCHING PARAMETERS =
; nblist update frequency =
nstlist                  = 10
; ns algorithm (simple or grid) =
ns_type                  = grid
; Periodic boundary conditions: xyz or none =
pbc                      = xyz
; nblist cut-off         =
rlist                    = 1.1

; OPTIONS FOR ELECTROSTATICS AND VDW =
; Method for doing electrostatics =
coulombtype              = cutoff
;rcoulomb-switch          = 0
rcoulomb                 = 1.1
; Dielectric constant (DC) for cut-off or DC of reaction field =
epsilon-r                = 1
; Method for doing Van der Waals =
vdw-type                 = switch
; cut-off lengths        =
rvdw-switch              = 0.9135 ;3 sigma
rvdw                     = 1.06575 ;3.5 sigma
; Spacing for the PME/PPPM FFT grid =
fourierspacing           = 0.1
; FFT grid size, when a value is 0 fourierspacing will be used =
fourier_nx               = 0
fourier_ny               = 0
fourier_nz               = 0
; EWALD/PME/PPPM parameters =
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no

; OUTPUT CONTROL OPTIONS =
; Output frequency for coords (x), velocities (v) and forces (f) =
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file =
nstlog                   = 500
nstenergy                = 100
; Output frequency and precision for xtc file =
nstxtcout                = 10
xtc-precision            = 1000

;OPTIONS FOR TEMPERATURE COUPLING
tc_grps                  = system
tcoupl                   = v-rescale
tau_t                    = 0.1 
ref_t                    = %(temperature)s
; Non-equilibrium MD 
acc-grps                 =
accelerate               =
freezegrps               = %(energyGroups)s  
freezedim                = %(freezedim)s 
energygrps	         = %(energyGroups)s
energygrp-excl	         = %(energyExcl)s
cos-acceleration         = 0
deform                   =
"""%vars()
    # Writes the simulation.mdp file
    write_file(os.path.join(work_path, 'simulation.mdp'), simmdp_lines)

# Generates the minimize.txt file lines
mintxt_lines = """4
6
0
"""
# Writes the minimize.txt file
write_file(os.path.join(work_path, 'minimize.txt'), mintxt_lines)

# Generates the nvt.txt file lines
nvttxt_lines = """4
7
9
0
"""
# Writes the nvt.txt file
write_file(os.path.join(work_path, 'nvt.txt'), nvttxt_lines)

# Generates the simulation.txt file lines
simtxt_lines = """0
"""
# Writes the simulation.txt file
write_file(os.path.join(work_path, 'simulation.txt'), simtxt_lines)

# Generates the makeindex.txt file lines
makeindex_lines = """q
    """
# Writes the makeindex.txt file
write_file(os.path.join(work_path, 'makeindex.txt'), makeindex_lines)
        
# Change directories to the work path
os.chdir(work_path)
# Write parameters to output 'Results.txt' file
f = open('Results.txt','w+')
f.write('System Parameters')
f.write('\n')
f.write('Number of Fluid Atoms: ' + str(nFluid) + ' Atoms')
f.write('\n')
f.write('Number of Default Wall Atoms: ' + str(nWall) + ' Atoms')
f.write('\n')
f.write('Other Wall Atoms: ')
f.write('\n')
for atom in atomTypes:
    f.write(str(atom[0]) + ': ' + str(atom[1]))
    f.write('\n')
f.write('Total Number of Wall Atoms: ' + str(nWall + totalAt) + ' Atoms')
f.write('\n')
f.write('Total Number of Atoms: ' + str(nWall + nFluid + totalAt) + ' Atoms')
f.write('\n')
f.write('Default Atom Wall-Fluid Epsilon: ' + str((fluid_epsilon*wall_epsilon)**0.5) + ' kJ')
f.write('\n')
f.write('Reduced Default Atom Wall-Fluid Epsilon (wall-fluid epsilon divided by fluid epsilon): ' + str(((fluid_epsilon*wall_epsilon)**0.5)/fluid_epsilon))
f.write('\n')
f.write('Default Atom Wall-Fluid Sigma: ' + str((sigma_fluid + sigma_wall)/2) + ' nanometers')
f.write('\n')
f.write('Reduced Default Atom Wall-Fluid Sigma (wall-fluid sigma divided by fluid sigma): ' + str(((sigma_fluid + sigma_wall)/2)/sigma_fluid))
f.write('\n')
f.write('\n')
f.write('Fluid')
f.write('\n')
f.write('Fluid Sigma: ' + str(sigma_fluid) + ' nanometers')
f.write('\n')
f.write('Fluid Epsilon: ' + str(fluid_epsilon) + ' kJ')
f.write('\n')
f.write('Fluid Mass: ' + str(fluid_mass) + ' kg')
f.write('\n')
f.write('Fraction of Wall Y Dimension Initially Filled with Fluid: ' + str(fluidFillFraction_y))
f.write('\n')
f.write('Initial Spacing Between Fluid Atoms: ' + str(factor) + ' Sigma')
f.write('\n')
f.write('Note: Fluid atoms are initially arranged in a body-centered cubic configuration')
f.write('\n')
f.write('\n')
f.write('Walls')
f.write('\n')
f.write('Default Wall Sigma: ' + str(sigma_wall) + ' nanometers')
f.write('\n')
f.write('Default Wall Epsilon: ' + str(wall_epsilon) + ' kJ')
f.write('\n')
f.write('Wall Atom lattice Parameter: ' + str(latticeParameter) + ' nanometers')
f.write('\n')
f.write('Default Wall Mass: ' + str(wall_mass) + ' kg')
f.write('\n')
f.write('Note: Wall atom mass is not relevant as the atoms are fixed in position')
f.write('\n')
f.write('Atoms in the Wall X Dimension: ' + str((wall_x/latticeParameter)) + ' Atoms')
f.write('\n')
f.write('Atoms in the Wall Y Dimension: ' + str((wall_y/latticeParameter)) + ' Atoms')
f.write('\n')
f.write('Atoms of Space in Between Walls: ' + str((dz/latticeParameter)) + ' Atoms')
f.write('\n')
f.write('\n')
f.write('Dynamics')
f.write('\n')
f.write('Dynamic Fluid Motion: ' + str(dynamic))
if (dynamic):
    f.write('\n')
    f.write('Acceleration: ' + str(acceleration) + ' Nanometers per Picosecond Squared')
f.write('\n')
f.write('Temperature: ' + str(temperature) + ' Kelvin')
f.write('\n')
f.write('\n')
f.write('Roughness')
f.write('\n')
f.write('Roughness: ' + str(roughness))
f.write('\n')
f.write('Number of Base Wall Layers Generated Regardless of Roughness Characteristics: ' + str(baseLayers))
f.write('\n')
if (roughness and not arbitraryPhysicalRoughness):
    f.write('Patterned Roughness: ' + str(roughness))
    f.write('\n')
    f.write('Roughness Type: ' + str(roughnessType))
    f.write('\n')
    f.write('Atom Height of Pillars: '+ str(height) + ' Atoms')
    f.write('\n')
    f.write('Atom Radius of Pillars from Center: ' + str(spacing) + ' Atoms')
    f.write('\n')
    f.write('Atom Spacing Between Pillar Centers: ' + str(spacing) + ' Atoms')
    f.write('\n')
    f.write('Note: These parameters are only relevant if there is no arbitrary physical roughness')
    f.write('\n')
if (arbitraryPhysicalRoughness):
    f.write('Arbitrary Physical Roughness: ' + str(arbitraryPhysicalRoughness))
    f.write('\n')
    f.write('Arbitrary Physical Roughness Height Range: ' + str(heightMin) + '-' + str(height) + ' Atoms')
    f.write('\n')
    f.write('Arbitrary Physical Roughness Radius Range: ' + str(radiusMin) + '-' + str(radiusMax) + ' Atoms')
    f.write('\n')
    f.write('Arbitrary Physical Roughness Spacing Range: ' + str(spacingMin) + '-' + str(spacingMax) + ' Atoms')
    f.write('\n')
if (chemicalRoughness):
    f.write('Patterned Chemical Roughness: ' + str(chemicalRoughness))
    f.write('\n')
    f.write('Additional Patterned Chemical Roughness Atoms Besides Default Wall Atoms:')
    f.write('\n')
    f.write('Format: [residue name, atom name, sigma, epsilon, mass, charge, minimum percentage, maximum percentage]')
    f.write('\n')
    for atom in roughAtoms:
        writeLine = '[' + str(atom[0]) + ', ' + str(atom[1]) + ', ' + str(atom[2]) + ', ' + str(atom[3]) + ', ' + str(atom[4]) + ', ' + str(atom[5]) + ', ' + str(atom[6]) + ', ' + str(atom[7]) + ']'
        f.write(str(writeLine))
        f.write('\n')
    f.write('Pillars are generated according to a repeating and sequential pattern of atom assignment')
    f.write('\n')
if (arbitraryChemicalRoughness):
    f.write('Arbitrary Chemical Roughness: ' + str(arbitraryChemicalRoughness))
    f.write('\n')
    f.write('Additional Arbitrary Chemical Roughness Atoms Besides Default Wall Atoms:')
    f.write('\n')
    f.write('Format: [residue name, atom name, sigma, epsilon, mass, charge, minimum percentage, maximum percentage]')
    f.write('\n')
    for atom in roughAtoms:
        writeLine = '[' + str(atom[0]) + ', ' + str(atom[1]) + ', ' + str(atom[2]) + ', ' + str(atom[3]) + ', ' + str(atom[4]) + ', ' + str(atom[5]) + ', ' + str(atom[6]) + ', ' + str(atom[7]) + ']'
        f.write(str(writeLine))
        f.write('\n')
    f.write('Atom is generated for a pillar if a random number between 1 and 100 is equal to or greater than minimum percentage but less than maximum percentage. If no atom is selected, the default wall atom is used.')
f.write('\n')
f.write('\n')
f.write('________________________________________________________________________________')
f.flush()
f.close()
# Change directories back to the base path
os.chdir(base_path)
