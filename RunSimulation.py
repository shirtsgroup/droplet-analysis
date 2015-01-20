
#===============================================================================
# Imports
#===============================================================================

import commands, math, os, pdb, scipy, shutil, sys, time

#===============================================================================
# Parameters
#===============================================================================

# Path to the directory containing the 'simulation' directory
base_path = '/User/username/directory/droplet-analysis/'
# Path to the Gromacs installation directory with no '/' at the end
gmx_path = '/User/username/directory/gromacs'
# Sets the maximum number of cores that will be utilized by Gromacs
cores = 8

#===============================================================================
# Subroutines
#===============================================================================

def run_Gmx(cmd):
    """
        Runs a command from command line that is provided as a string.
    """
 
    # Executes the command
    output = commands.getoutput(cmd)
    
    # If there is a 'Fatal error', show which command gave the error
    if ((output.count('Fatal error') + output.count('Segmentation fault')) > 0):
        print ('There was a fatal error in the following command:')
        print (cmd)
        print (os.getcwd())
        print ('Exiting...')
        sys.exit(1)
    pass
    

#===============================================================================
# Main
#===============================================================================

# Sets start time
start_time = time.time()
# Additional parameters (DO NOT CHANGE!)
# The work path for running the simulation
work_path = base_path + 'simulation/'
# Files names of important input files to be generated
minimize_mdp_file = 'minimize.mdp'
nvt_mdp_file = 'nvt.mdp'
gro_file = 'system.gro'
top_file = 'system.top'
ndx_file = 'index.ndx'

# Changes directory to the base path
os.chdir(base_path)

# Runs a Python module to generate the initial system
print ('Generating Initial System...')
os.system('python ' + base_path + 'SetupSystem.py')

# Changes directory to the work path
os.chdir(work_path)

# Makes the index.ndx file
print ('Making Index File...')

# Runs 'make_ndx -f system.gro -o index.ndx < makeindex.txt'
mkindex = '%(gmx_path)sbin/make_ndx -f %(gro_file)s -o %(ndx_file)s < makeindex.txt'%vars()
output = run_Gmx(mkindex)

# Determines generation time
generation_time = time.time()

# Runs the energy minimization component of the simulation
print ('Energy Minimization...')

# Runs 'grompp -f minimize.mdp -c system.gro -n index.ndx -p system.top -o minimize.tpr -maxwarn 100'
grompp = '%(gmx_path)sbin/grompp -f %(minimize_mdp_file)s -c %(gro_file)s -n %(ndx_file)s -p %(top_file)s -o minimize.tpr -maxwarn 100'%vars()
output = run_Gmx(grompp)

# Runs 'mdrun -deffnm minimize -nt (Number of Cores)'
mdrun = '%(gmx_path)sbin/mdrun -deffnm minimize -nt %(cores)s'%vars()
output = run_Gmx(mdrun)
# Runs 'g_energy -f minimize.edr -o minimize.xvg'
g_energy = '%(gmx_path)sbin/g_energy -f minimize.edr -o minimize.xvg < minimize.txt'%vars()
output = run_Gmx(g_energy)

# Determines minimization time
minimization_time = time.time()

# Runs the NVT ensemble
print ('NVT Ensemble...')

# Runs 'grompp -f nvt.mdp -c minimize.gro -n index.ndx -p system.top -o nvt.tpr -maxwarn 100'
grompp = '%(gmx_path)sbin/grompp -f %(nvt_mdp_file)s -c minimize.gro -n %(ndx_file)s -p %(top_file)s -o nvt.tpr -maxwarn 100'%vars()
output = run_Gmx(grompp)

# Runs 'mdrun -deffnm nvt -nt (Number of Cores)'
mdrun = '%(gmx_path)sbin/mdrun -deffnm nvt -nt %(cores)s'%vars()
output = run_Gmx(mdrun)

# Runs 'g_energy -f nvt.edr -o nvt.xvg'
g_energy = '%(gmx_path)sbin/g_energy -f nvt.edr -o nvt.xvg < nvt.txt'%vars()
output = run_Gmx(g_energy)

# Determines nvt time
nvt_time = time.time()

# Runs the production simulation
print ('Running Simulation...')

# Runs 'grompp -c nvt.gro -p system.top -f simulation.mdp -o simulation.tpr -maxwarn 100'
grompp = '%(gmx_path)sbin/grompp -c nvt.gro -p %(top_file)s -f simulation.mdp -o simulation.tpr -maxwarn 100'%vars()
output = run_Gmx(grompp)

# Runs 'mdrun -deffnm simulation -nt (Number of Cores)'
mdrun = '%(gmx_path)sbin/mdrun -deffnm simulation -nt %(cores)s'%vars()
output = run_Gmx(mdrun)

# Determines simulation time
simulation_time = time.time()

# Processes the simulation output and enerates a .gro file for every frame so contact angles can be determined
print ('Processing Output...')

# Runs 'trjconv -s simulation.tpr -f simulation.xtc -o analyze.gro -ndec 3'
trjconv = '%(gmx_path)sbin/trjconv -s simulation.tpr -f simulation.xtc -o analyze.gro -ndec 3 < simulation.txt'%vars()
output = run_Gmx(trjconv)

# Change directory to the base path so that 'ComputeContactAngle.py' can be run
os.chdir(base_path)

# Determines contact angle calculation time
conversion_time = time.time()

# Runs a Python module to generate contact angle data from the simulation output
print ('Calculating Contact Angles...')
os.system('python ' + base_path + 'ProcessResults.py')
os.system('python ' + base_path + 'AnalyzeResults.py')

# Determines contact angle calculation time
analysis_time = time.time()

# Simulation steps are now complete, data has been generated, and this module can be closed
print ('Simulation Complete!')

# Determines total time
end_time = time.time()

os.chdir(base_path + 'simulation/')
f = open('Results.txt','a')
f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')
f.write('Simulation Time Analysis')
f.write('\n')
f.write('System Generation Time: ' + str((generation_time - start_time)/(60*60)) + ' Hours')
f.write('\n')
f.write('Energy Minimization Time: ' + str((minimization_time - generation_time)/(60*60)) + ' Hours')
f.write('\n')
f.write('NVT Ensemble Time: ' + str((nvt_time - minimization_time)/(60*60)) + ' Hours')
f.write('\n')
f.write('Production Simulation Time: ' + str((simulation_time - nvt_time)/(60*60)) + ' Hours')
f.write('\n')
f.write('Trajectory Conversion Time: ' + str((conversion_time - simulation_time)/(60*60)) + ' Hours')
f.write('\n')
f.write('Analysis Time: ' + str((analysis_time - conversion_time)/(60*60)) + ' Hours')
f.write('\n')
f.write('Total Time: ' + str((end_time - start_time)/(60*60)) + ' Hours')
f.write('\n')
f.write('________________________________________________________________________________')
f.write('\n')
f.write('\n')
f.flush()
f.close()

os.chdir(base_path)

# Moves desired output files to the base directory
shutil.move(base_path + 'simulation/Results.txt', base_path + 'Results.txt')
shutil.move(base_path + 'simulation/Data.pkl', base_path + 'Data.pkl')
shutil.move(base_path + 'simulation/minimize.gro', base_path + 'minimize.gro')
shutil.move(base_path + 'simulation/minimize.log', base_path + 'minimize.log')
shutil.move(base_path + 'simulation/minimize.xvg', base_path + 'minimize.xvg')
shutil.move(base_path + 'simulation/nvt.gro', base_path + 'nvt.gro')
shutil.move(base_path + 'simulation/nvt.log', base_path + 'nvt.log')
shutil.move(base_path + 'simulation/nvt.xvg', base_path + 'nvt.xvg')
shutil.move(base_path + 'simulation/simulation.gro', base_path + 'simulaton.gro')
shutil.move(base_path + 'simulation/simulation.log', base_path + 'simulaton.log')
shutil.move(base_path + 'simulation/system.gro', base_path + 'system.gro')

# Removes the simulation directory
os.system('rm -r ' + base_path + 'simulation/')
