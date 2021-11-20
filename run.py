#!/usr/bin/env

import os

##START SETUP
def get_mut_script(mutation):
    '''
    Generate a script to invoke pmx performing a requested SNP and make an associated gro file
    :param mutation: string naming the mutation eg: s450l
    :return: script as a sting
    '''

    mutate_script = '''#! /bin/bash
    
source {}

sim={}
timestep={}         # e.g. 10ns
mutation={}         # e.g. s428c

{}/pmx/scripts/mutate.py -f $sim/rpob-5uh6-3-$sim-$timestep-chainC.pdb -o $mutation/rpob-5uh6-3-$sim-$timestep-chainC-$mutation-0 -ff amber99sb-star-ildn-mut.ff -script mutations.txt

perl -pi -e 's/GLU C1126/CGLUC1126/g' $mutation/rpob-5uh6-3-$sim-$timestep-chainC-$mutation-0.pdb

grep CRYST $sim/rpob-5uh6-3-$sim-$timestep-chainAB.pdb > $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.pdb
grep ATOM $sim/rpob-5uh6-3-$sim-$timestep-chainAB.pdb >> $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.pdb
grep ATOM $mutation/rpob-5uh6-3-$sim-$timestep-chainC-$mutation-0.pdb >> $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.pdb
grep ATOM $sim/rpob-5uh6-3-$sim-$timestep-chainDEF.pdb >> $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.pdb
grep ATOM $sim/rpob-5uh6-3-$sim-$timestep-not-protein.pdb >> $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.pdb
echo "END" >> $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.pdb

gmx editconf -f $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.pdb -o $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.gro

    '''.format(GMXRC, sim, timestep, mutation, pmx_root)

    return mutate_script

def get_num_atoms(gro_file):
    '''
    Function to extract the number of atoms from a gro file. It reads the second line.
    :param gro_file:
    :return: int for the number of atoms.
    '''
    with open(gro_file) as f:
        lines = f.readlines()
        num_atoms = int(lines[1])
    return num_atoms

def get_hack_apo_script(mutation, num_atoms):
    '''
    Function genreates a script, the script removes the drug from a gro file to hack the apo structure.
    :param mutation: str , name of the mutation
    :param num_atoms: int, the number of atoms in the og gro file we want to make apo.
    :return: script as string
    '''
    hack = '''#! /bin/bash
    
sim={}
timestep={}         # e.g. 10ns
mutation={}         # e.g. s428c

# hack an apo structure
grep -v RFP $mutation/rpob-5uh6-3-$sim-$timestep-$mutation.gro > $mutation/rpob-5uh6-3-$sim-$timestep-$mutation-apo.gro

# correct the number of atoms (subtract 119 for RFP)
perl -pi -e 's/{}/{}/g' $mutation/rpob-5uh6-3-$sim-$timestep-$mutation-apo.gro
    '''.format(sim, timestep, mutation, num_atoms, num_atoms-119)

    return hack

def get_gen_topo_script(filestem, mutation):
    '''
    Function to generate topology of mutant chain C
    :param filestem: str, identifying string that is prefix to most files.
    :param mutation: str, name of the mutation.
    :return:
    '''
    topo = '''#! /bin/bash

source {}

filestem={}
mutation={}

echo "1
1" | gmx pdb2gmx -f ../seeds/$mutation/$filestem-chainC-$mutation-0.pdb -o $filestem-chainC-$mutation-1.gro -p $filestem-chainC-$mutation.top -i $filestem-chainC-$mutation\_posre.itp

{}/pmx/scripts/generate_hybrid_topology.py -p $filestem-chainC-$mutation.top -o $filestem-chainC-$mutation-1 -ff amber99sb-star-ildn-mut.ff

    '''.format(GMXRC, filestem, mutation, pmx_root)

    return topo

def get_itp_script(mutation):
    '''
    Function that generates script to convert top files to itp. We can then bundle the itps later into the full topo.
    :param mutation: str, name of the mutation
    :return: script as a sting
    '''

    #end of itps currently bad

    itp = '''#! /bin/bash

mutation={}
time={}

awk 'NR>1 {{print $0}}' rpob-5uh6-3-md-1-$time-chainC-$mutation-1.top > rpob-5uh6-rnap_beta_$mutation.itp
perl -pi -e 's/Protein_chain_C/RNAP_BETA/g' rpob-5uh6-rnap_beta_$mutation.itp
perl -pi -e 's/3-md-1-'$time'-chainC-/rnap_beta_/g' rpob-5uh6-rnap_beta_$mutation.itp
cp rpob-5uh6-3-md-1-$time-chainC-$mutation\_posre.itp rpob-5uh6-rnap_beta_$mutation\_posre.itp

    '''.format(mutation, timestep)

    return itp

def gen_topo(mutation):
    '''
    Wrapper function around all the script that must be executed to generate the topology
    :param mutation: str, name of the mutation
    :return: None
    '''
    os.chdir('../generate-topology')
    filestem = 'rpob-5uh6-3-{}-{}'.format(sim, timestep)
    topo_script = get_gen_topo_script(filestem, mutation.name)
    os.system(topo_script)
    itp_script = get_itp_script(mutation.name)
    os.system(itp_script)


def mutate(mutation):
    '''
    Wrapper function around all the script that must be excuted to build the mutant system
    :param mutation: str, name of the mutation
    :return: None
    '''
    os.mkdir(os.path.join(cwd, mutation.name))
    os.mkdir(os.path.join(cwd, mutation.name, sim))
    os.mkdir(os.path.join(cwd, mutation.name, sim, timestep))
    mutation.write_txt()
    mut_script = get_mut_script(mutation.name)
    os.system(mut_script)
    gro_file = os.path.join(cwd, mutation.name, 'rpob-5uh6-3-{}-{}-{}.gro'.format(sim, timestep, mutation.name))
    num_atoms = get_num_atoms(gro_file)
    hack_apo = get_hack_apo_script(mutation.name, num_atoms)
    os.system(hack_apo)

def clean_topo(mutation):
    '''
    Function to remove some extra lines that are not needed at the end of the topology file
    If these lines are left in there is some conflict about what force fields to use.
    :param mutation:
    :return:
    '''

    top_file = os.path.join(cwd, '../generate-topology/rpob-5uh6-rnap_beta_{}.itp'.format(mutation.name))
    with open(top_file) as f:
        to_write = f.readlines()[:-10]
        to_write = ''.join(to_write)

    with open(top_file, 'w') as f:
        f.write(to_write)

def setup(mutation):
    '''
    Wrapper arrond the function required to setup the mutant systems.
    :param mutation: str, name of the mutant.
    :return:
    '''
    mutate(mutation)
    gen_topo(mutation)
    clean_topo(mutation)

def write_full_tops(mutation, leg):
    '''
    Function to write top file combining all itps
    :param mutation: str, name of the nutation
    :param leg: str, name of simulation leg (apo/rfp)
    :return: None
    '''

    if leg == 'apo':
        include = '; NO DRUG'
        sys_name = 'un-bound'
        mol = include
    else:
        include = '#include \"rpob-5uh6-rfp-1.itp\"'
        sys_name = 'bound'
        mol = 'RFP                  1'


    template = '''; Include forcefield parameters
#include "amber99sb-star-ildn-mut.ff/forcefield.itp"

#include "../../build/atomtypes.itp"

; Include water topology
#include "amber99sb-star-ildn-mut.ff/tip3p.itp"

; Include topology for ions
#include "amber99sb-star-ildn-mut.ff/ions.itp"

#include "../../build/5uh6-rnap_alpha_1.itp"
#include "../../build/5uh6-rnap_alpha_2.itp"
#include "../../generate-topology/rpob-5uh6-rnap_beta_s450l.itp"
#include "../../build/5uh6-rnap_beta_prime.itp"
#include "../../build/5uh6-rnap_omega.itp"
#include "../../build/5uh6-rnap_sigma.itp"
#include "../../build/5uh6-dna_1.itp"
#include "../../build/5uh6-dna_2.itp"
#include "../../build/5uh6-rna_2nt.itp"
{}

[ system ]
; Name
{} rpoB complex with rifampicin {} (PDB:5UH6)

[ molecules ]
; Compound         #mols
RNAP_ALPHA_1         1
RNAP_ALPHA_2         1
RNAP_BETA            1
RNAP_BETA_PRIME      1
RNAP_OMEGA           1
RNAP_SIGMA           1
DNA_1                1
DNA_2                1
RNA_2NT              1
{}
ZN                   2
MG                   1
SOL             114838
NA                 127
        '''.format(include, mutation, sys_name, mol)


    file_name = 'rpob-5uh6-{}-{}.top'.format(mutation, leg)
    with open(os.path.join(cwd, mutation, file_name), 'w') as f:
        f.write(template)
##END SETUP

class Mutation():
    def __init__(self, og, target, res_id, shift):
        '''
        Class to hold infomation about the mutation
        :param og:str, denoting the letter of the original amino acid
        :param target:str, denoting the letter of the amino acid we want to mutate to
        :param res_id:int, denoting the res id in the full stucture
        :param shift:int, denoting the differnce between the res id in the full structure and the res id in just Chain C
        '''
        self.og = og
        self.target = target
        self.res_id = res_id
        self.shift = shift
        self.true_id = res_id-shift
        self.name = '{}{}{}'.format(og, res_id, target)

    def write_txt(self):
        '''
        Write a file containing id of residue to mutate and letter of amino acid to mutate in to.
        :return: None
        '''
        path = os.path.join(cwd, 'mutations.txt')
        f = open(path, "w")
        f.write("{} {}".format(self.true_id, self.target))
        f.close()

def main():
    '''
    Main function to run all components of the pipeline (setup and run)
    :return: None
    '''

    #Iterate over requested muataions
    #this preps both apo and drug legs
    for mut in mutations:

        #from seed files build mutant inputs
        setup(mut)
        for leg in legs:
            write_full_tops(mut.name, leg)

    #leg_lines={{} for leg in legs}
    #for mut in mutations:
    #    for leg in legs:
    #        #Collect run line for all mutants
    #        run_lines = main(mut)
    #        leg_lines[leg][mut] = run_lines


    #write run_lines to file

if __name__ == "__main__":
    ###GLOBAL OPTIONS

    #must run this is user env
    #export GMXLIB="/gpfs/alpine/scratch/adw62/chm155/Phil_work/pmx/pmx/data/mutff45"
    #must load gromacs module

    #What mutation do we want to look at
    mut1 = Mutation('s', 'l', 450, 21)

    mutations = [mut1]

    legs = ['apo', 'rfp']

    #these two variable define the seed strutures used. Provided by Phil Fowler
    sim = 'md-1'
    timestep = '0ns'

    #Where is GMX installed
    GMXRC = '/sw/summit/spack-envs/base/opt/linux-rhel8-ppc64le/' \
            'gcc-8.3.1/gromacs-2021.2-ijh7cfdbr6ckwxu5jruhe5am2efmogm3/bin/GMXRC'

    #Where is PMX installed
    pmx_root = '/gpfs/alpine/scratch/adw62/chm155/Phil_work/pmx'
    ###########

    cwd = os.getcwd()
    main()




