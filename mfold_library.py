import os
import re
import string
import subprocess


class Strand:
    allowed_bases = set('ATCG')
    allowed_constraints = set(string.ascii_letters + string.digits)
    base_pair = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def __init__(self, bases, constraints):
        """
        Args:
            bases: A string representing the bases in a strand
            constraints: A list of Regions representing the structure of a strand
        """
        self.bases = bases.upper()
        self.constraints = constraints
        if set(self.bases) > Strand.allowed_bases:
            raise TypeError('The selected bases contain letters '
                          + 'other than A, T, C, and G: ' + bases)
        if set(self.constraints) > Strand.allowed_constraints:
            raise TypeError('The selected constraints contain '
                          + 'non-alphanumeric characters: ' + constraints)

    def complement(bases):
        return "".join([Strand.base_pair[base] for base in bases])

class Region:
    def __init__(self, name, length):
        # the name that represents the region, e.g. 'A3' -> 'A'
        self.name = name
        # the length of the region
        self.length = length

class Mfold:
    energy_string = ' Initial dG = '
    linker_sequence = 'LLL'
    output_suffixes = ['.aux', '.cmd', '.con', '.log', '.pnt', '.sav', '.seq'
            '-local.pnt', '-local.seq', '.ann', '.ct', '.ps', '.det', '.out',
            '.h-num', '.plot', '.pdf', '.ss-count', '-temp.det', '-temp.out']


    def __init__(self, output_folder='', mfold_command=''):
        self.folder = output_folder
        self.command = mfold_command


    def run(self, strand1, strand2, sequence_file='a.seq', settings_file='a.aux'):
        seq_path = os.path.join(self.folder, sequence_file)
        set_path = os.path.join(self.folder, settings_file)

        with open(seq_path, 'w') as seqfile:
            seqfile.write(strand1.bases + Mfold.linker_sequence + strand2.bases)
        with open(set_path, 'w') as setfile:
            for constraint in Mfold.get_constraints(strand1, strand2):
                setfile.write(constraint)

        #subprocess.run([self.command, f'SEQ={seq_path}', f'AUX={set_path}'],
        #        cwd=self.folder)


    def clean(self, file_prefix):
        for suffix in Mfold.output_suffixes:
            file_path = os.path.join(self.folder, f'{file_prefix}{suffix}')
            if os.path.exists(file_path):
                os.remove(file_path)


    def get_energy(self, details_file='a.det'):
        details_path = os.path.join(self.folder, details_file)
        if os.path.exists(details_path):
            with open(details_path, 'r') as detfile:
                for line in detfile:
                    if line.startswith(Mfold.energy_string):
                        return float(line[len(Mfold.energy_string):])
        return None
            

    def get_constraints(strand1, strand2):
        constraints = []
        all_regions = {}

        curr_index = 0
        for region in strand1.constraints:
            all_regions[region.name] = (curr_index, curr_index + region.length)
            curr_index += region.length

        curr_index += 3
        for region in strand2.constraints:
            all_regions[region.name] = (curr_index, curr_index + region.length)
            curr_index += region.length


        for region in all_regions:
            if region.isupper():
                constraints.append(
                        f'P {all_regions[region.lower()][0]} {all_regions[region.lower()][1]} '
                        + f'{all_regions[region][0]} {all_regions[region][1]}')
        return constraints


class EnergyMatrix:
    def __init__(self, mfold, strands):
        self.mfold = mfold
        self.strands = strands
        self.matrix = [[None for strand1 in strands] for strand2 in strands]

    def create(self):
        for i, strand1 in enumerate(self.strands):
            for j, strand2 in enumerate(self.strands):
                self.mfold.clean(f'{i}_{j}')
                self.mfold.run(strand1, strand2, f'{i}_{j}.seq', f'{i}_{j}.aux')
        self.matrix = [[self.mfold.get_energy(f'{i}_{j}.det')
                for i in range(len(self.strands))] 
                    for j in range(len(self.strands))]
