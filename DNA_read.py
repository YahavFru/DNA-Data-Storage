import re

is_print = True

def pam_finder(dna_seq, pam):
    pam_indices = re.finditer('.GG', dna_seq)
    return [pam.end() for pam in pam_indices] # Returns the start position of each pam sequence
    
# TODO: a function that finds the PAMs
# input: sequence + pam
# output: a list of pam indices 

def main(dna_seq, pam):
    pam_indices = pam_finder(dna_seq, pam)
    if is_print:
        print(f'The PAM indeces are: ${pam_indices}, The PAM sequences are: ${[dna_seq[pam_index] for pam_index in pam_indices]}')
    dna_data_storage_proccess.__init__(dna_data_storage_proccess, dna_seq, pam_indices)    
    dna_data_storage_proccess.write(dna_data_storage_proccess)


# TODO: joint class for encoder+writer+reader+decoder
# Shared metadata: DNA seq, PAM indices (or at least PAM sequence + logic to generate a list of indices)
# functions:
# encoder: binary message -> "ideal" edited sequence/list of PAM indices to edit (V)
# writer: simulate the enzymatic process. DNA seq + list of PAM indices to edit -> list of *many* edited DNA sequence
#         parameters: edit probability, number of copies, other noise
# reader: simulate NGS. list of *many* edited DNA sequence + DNA seq (unedited) -> list of edited PAM indices
#         parameters: logic for identifying editing from a noisy signal
# decoder: "ideal" edited sequence/list of PAM indices to edit -> binary message

class dna_data_storage_proccess:
    bit_list = [1] #array of booleans / 0s or 1s
    dna_seq = '' #original DNA sequence, string
    tagged_ideal_seq = ''
    grna_indexes_for_edit = []
    grna_indexes = [] #array of indexes for each gRNA block on DNA
    tagged_dna_seqs = [] #array of dna seqs, each edited.

    def __init__(self, origin_dna, pam_end_indices): #Defines all inputs for the class
        self.dna_seq = origin_dna
        for index in pam_end_indices:
            self.grna_indexes.append(index)
        
    def encode(self):
        if len(self.bit_list) != len(self.grna_indexes):
            print(f'Error: Bit list length ({len(self.bit_list)}) does not match the number of pam sequences present ({len(self.grna_indexes)})')
        else:
            for index in self.grna_indexes:
                if self.bit_list[index]: # is this bit 1? i.e, does it require editing?
                    self.grna_indexes_for_edit.append(index)
                    for base_jump in range(20):
                        if self.dna_seq[index + base_jump] == 'C':
                            self.tagged_ideal_seq += 'T'
                        else: 
                            self.tagged_ideal_seq += self.dna_seq[index + base_jump]


main('AGGCCCCCCCCCCCCCCCCCCCC', 'NGG')
