import os

#helper functions
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


def make_all_dirs(muts_to_build):
    for mut in muts_to_build:
        #os.mkdir(os.path.join(seeds_dir, mut))
        for leg in legs:
            os.mkdir(os.path.join(seeds_dir, mut, leg))
            for state in lambdas:
                os.mkdir(os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num)))
                for rep in reps:
                    os.mkdir(os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep)))
                    os.mkdir(os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep' + str(rep), 'equilibration'))
                    os.mkdir(os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep' + str(rep), 'simulation'))


def main():
    make_all_dirs(all_muts)
    if os.path.isdir('./bash_scripts'):
        pass
    else:
        os.makedirs('./bash_scripts')

    for i, group in enumerate(chunker(all_muts, 5)):
        #make sub script
        with open('./sub{}.sh'.format(i), 'w') as f:

            header_ = header.format(i, i, i)
            f.write(header_)

            # gromp0
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            if state.win_type == 'q':
                                continue
                            else:
                                mdp = os.path.join(ref_dir, 'EM-{}.mdp'.format(state.num))
                                gro = os.path.join(seeds_dir, mut, 'rpob-5uh6-3-md-1-0ns-{}-{}.gro'.format(mut, leg))
                                top = os.path.join(seeds_dir, mut,
                                                   'rpob-5uh6-{}-{}-{}.top'.format(mut, leg, state.win_id))
                                out = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration', 'EM')
                                run_line = ' '.join([jsrun, gromp, '-f', mdp, '-c', gro, '-p', top,
                                                     '-o', out, '-po', out + '.mdp', '-maxwarn 2 &\n'])
                                f.write(run_line)
            f.write('wait\n')


            # minimization of mutant
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            if state.win_type == 'q':
                                continue
                            else:
                                out = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration', 'EM')
                                tpr = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration', 'EM.tpr')
                                run_line = ' '.join([jsrun, gmx_in_cpu.format('false', tpr), gmx_out.format(out), '&\n'])
                                f.write(run_line)
            f.write('wait\n')


            #gromp1
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            if state.win_type == 'q':
                                continue
                            else:
                                root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                                mdp = os.path.join(ref_dir, 'alchembed.mdp')
                                gro = os.path.join(root, 'EM.gro')
                                top = os.path.join(seeds_dir, mut, 'rpob-5uh6-{}-{}-{}.top'.format(mut, leg, state.win_id))
                                out = os.path.join(root, 'alchembed')
                                run_line = ' '.join([jsrun, gromp, '-f', mdp, '-c', gro, '-p', top,
                                                    '-o', out, '-po', out+'.mdp', '-maxwarn 2 &\n'])
                                f.write(run_line)
            f.write('wait\n')


            #alchembed
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            if state.win_type == 'q':
                                continue
                            else:
                                root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                                out = os.path.join(root, 'alchembed')
                                tpr = os.path.join(root, 'alchembed.tpr')
                                run_line = ' '.join([jsrun, gmx_in.format('false', tpr), gmx_out.format(out), '&\n'])
                                f.write(run_line)
            f.write('wait\n')


            #extract gros
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            if state.win_type == 'q':
                                continue
                            else:
                                root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                                in_ = os.path.join(root, 'alchembed')
                                xtc = in_+'.xtc'
                                gro = os.path.join(root, 'bed_out.gro')
                                tpr = os.path.join(root, 'alchembed.tpr')
                                script = ' '.join([trajconv.format(xtc, str(state.bed_val), gro, tpr)])
                                with open(os.path.join('./bash_scripts', '{}{}{}{}'.format(mut, leg, state.num, rep)), 'w') as s:
                                    s.write(script)
                                os.system('chmod +x ./bash_scripts/{}{}{}{}'.format(mut, leg, state.num, rep))
                                run_line = ' '.join([jsrun, './bash_scripts/{}{}{}{}'.format(mut, leg, state.num, rep), '&\n'])
                                f.write(run_line)
            f.write('wait\n')


            #copy bed grow files from vdw runs to q runs
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            if state.win_type == 'q':
                                if int(state.num) < 3:
                                    src = os.path.join(seeds_dir, mut, leg, 'LAMBDA_'+str(3), 'rep'+str(rep), 'equilibration', 'bed_out.gro')
                                    dest = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                                    run_line = ' '.join([jsrun, 'cp {} {}'.format(src, dest), '&\n'])
                                    f.write(run_line)
                                else:
                                    src = os.path.join(seeds_dir, mut, leg, 'LAMBDA_'+str(9), 'rep'+str(rep), 'equilibration', 'bed_out.gro')
                                    dest = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                                    run_line = ' '.join([jsrun, 'cp {} {}'.format(src, dest), '&\n'])
                                    f.write(run_line)
                            else:
                                continue
            f.write('wait\n')

            #gromp warming sims
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                            mdp = os.path.join(ref_dir, 'warm-short-{}.mdp'.format(state.num))
                            gro = os.path.join(root, 'bed_out.gro')
                            top = os.path.join(seeds_dir, mut, 'rpob-5uh6-{}-{}-{}.top'.format(mut, leg, state.win_id))
                            out = os.path.join(root, 'warm')
                            run_line = ' '.join([jsrun, gromp, '-f', mdp, '-c', gro, '-p', top,
                                                 '-o', out, '-po', out + '.mdp', '-maxwarn 2 &\n'])
                            f.write(run_line)
            f.write('wait\n')

            #warming
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                            out = os.path.join(root, 'warm')
                            tpr = os.path.join(root, 'warm.tpr')
                            run_line = ' '.join([jsrun, gmx_in.format('false', tpr), gmx_out.format(out), '&\n'])
                            f.write(run_line)
            f.write('wait\n')

            # gromp NVT eq
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                            mdp = os.path.join(ref_dir, 'NVT-{}'.format(state.num))
                            gro = os.path.join(root, 'warm.gro')
                            top = os.path.join(seeds_dir, mut, 'rpob-5uh6-{}-{}-{}.top'.format(mut, leg, state.win_id))
                            cpt = os.path.join(root, 'warm.cpt')
                            out = os.path.join(root, 'NVT')
                            run_line = ' '.join([jsrun, gromp, '-f', mdp, '-c', gro, '-p', top, '-t', cpt,
                                                 '-o', out, '-po', out + '.mdp', '-maxwarn 2 &\n'])
                            f.write(run_line)
            f.write('wait\n')

            # run NVT eq
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                            out = os.path.join(root, 'NVT')
                            tpr = os.path.join(root, 'NVT.tpr')
                            run_line = ' '.join([jsrun, gmx_in.format('false', tpr), gmx_out.format(out), '&\n'])
                            f.write(run_line)
            f.write('wait\n')

            # gromp NPT eq
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                            mdp = os.path.join(ref_dir, 'NPT-{}'.format(state.num))
                            gro = os.path.join(root, 'NVT.gro')
                            top = os.path.join(seeds_dir, mut, 'rpob-5uh6-{}-{}-{}.top'.format(mut, leg, state.win_id))
                            cpt = os.path.join(root, 'NVT.cpt')
                            out = os.path.join(root, 'NPT')
                            run_line = ' '.join([jsrun, gromp, '-f', mdp, '-c', gro, '-p', top, '-t', cpt,
                                                 '-o', out, '-po', out + '.mdp', '-maxwarn 2 &\n'])
                            f.write(run_line)
            f.write('wait\n')

            # run NPT eq
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'equilibration')
                            out = os.path.join(root, 'NPT')
                            tpr = os.path.join(root, 'NPT.tpr')
                            run_line = ' '.join([jsrun, gmx_in.format('false', tpr), gmx_out.format(out), '&\n'])
                            f.write(run_line)
            f.write('wait\n')

            # gromp prod
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep))
                            mdp = os.path.join(ref_dir, 'prod-{}'.format(state.num))
                            gro = os.path.join(root, 'equilibration', 'NPT.gro')
                            top = os.path.join(seeds_dir, mut, 'rpob-5uh6-{}-{}-{}.top'.format(mut, leg, state.win_id))
                            cpt = os.path.join(root, 'equilibration', 'NPT.cpt')
                            out = os.path.join(root, 'simulation', 'prod')
                            run_line = ' '.join([jsrun, gromp, '-f', mdp, '-c', gro, '-p', top, '-t', cpt,
                                                 '-o', out, '-po', out + '.mdp', '-maxwarn 2 &\n'])
                            f.write(run_line)
            f.write('wait\n')

            # run production
            for mut in group:
                for leg in legs:
                    for state in lambdas:
                        for rep in reps:
                            root = os.path.join(seeds_dir, mut, leg, 'LAMBDA_{}'.format(state.num), 'rep'+str(rep), 'simulation')
                            out = os.path.join(root, 'prod')
                            tpr = os.path.join(root, 'prod.tpr')
                            run_line = ' '.join([jsrun, gmx_in.format('false', tpr), gmx_out.format(out), '&\n'])
                            f.write(run_line)
            f.write('wait\n')

class Window():
    def __init__(self, num, lam_val, win_type, win_id, bed_val):
        self.num = str(num)
        self.lam_val = str(lam_val)
        self.win_type = win_type
        self.win_id = win_id
        self.bed_val = str(bed_val)

class Mutation():
    def __init__(self, shift_, og, target, res_id):
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
        self.shift = shift_
        self.true_id = res_id - shift_
        self.name = '{}{}{}'.format(og, res_id, target)

def get_mutations(mut_file):
    all_muts = []
    with open(mut_file) as f:
        for line in f:
            data = line.strip('\n').split(' ')
            data = [data[0].lower(), data[2].lower(), int(data[1])]
            all_muts.append(Mutation(shift, *data))

    return all_muts


if __name__ == '__main__':
    # define mutants
    # Build class for each mutation to carry its info
    shift = 21
    mutations = get_mutations('./all_muts.dat')

    # filer the mutations down, useful to test one at a time
    filter_ = range(0, 5)
    all_muts = [mutations[x].name for x in filter_]

    print(all_muts)

    # build lambdas
    state_id = [11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13]
    state_type = ['q', 'q', 'q', 'vdw', 'vdw', 'vdw', 'vdw', 'vdw', 'vdw', 'vdw', 'q', 'q', 'q']
    q_lamb = [0.0, 0.5, 1.0]
    vdw_lamb = [0.00, 0.125, 0.25, 0.50, 0.75, 0.875, 1.0]
    lambda_ = q_lamb + vdw_lamb + q_lamb
    bed = [0.7 * x for x in lambda_]
    assert len(state_id) == len(state_type)
    assert len(state_type) == len(lambda_)
    lambdas = [Window(i, x, y, z, k) for i, (x, y, z, k) in enumerate(zip(lambda_, state_type, state_id, bed))]

    # define replicas
    reps = [0, 1, 2, 3, 4]


    legs = ['apo', 'rfp']

    # dir locations
    seeds_dir = '/gpfs/alpine/scratch/adw62/chm155/Phil_work/test_auto/seeds'
    ref_dir = '/gpfs/alpine/scratch/adw62/chm155/Phil_work/test_auto/reference'

    # run commands
    jsrun = 'jsrun -n 1 -a 1 -c 7 -g 1 -bpacked:7'
    gromp = 'gmx_mpi grompp'
    gmx_in = 'gmx_mpi mdrun -ntomp 28 -nb gpu -bonded gpu -pme gpu -nstlist 200 -pmefft gpu -dlb no' \
             ' -gpu_id 0 -pin on -pinoffset 0 -pinstride 1 -tunepme {} -s {}'
    gmx_in_cpu = 'gmx_mpi mdrun -ntomp 28'
    gmx_out = '-deffnm {}'
    trajconv = 'echo \"0\" | gmx_mpi trjconv -f {} -dump {} -o {} -s {}'


    header = '''#!/bin/bash
#BSUB -P CHM155_001
#BSUB -W 720
#BSUB -alloc_flags "gpudefault"
#BSUB -q batch
#BSUB -nnodes 109
#BSUB -J {}
#BSUB -o {}_%J.out
#BSUB -e {}_test_%J.out

module load gcc/9.3.0
module load cuda/11.0.3
module load fftw/3.3.9

export GROMACS_DIR=/gpfs/alpine/proj-shared/chm155/benjha/gromacs/2021.2
export PATH=$GROMACS_DIR/bin:$PATH
export LD_LIBRARY_PATH=$GROMACS_DIR/lib64:$LD_LIBRARY_PATH
export GMXLIB="/gpfs/alpine/scratch/adw62/chm155/Phil_work/pmx/pmx/data/mutff45"

export PATH="/gpfs/alpine/scratch/adw62/chm155/Phil_work/py/bin:$PATH"

export OMP_NUM_THREADS=28
'''

    main()

