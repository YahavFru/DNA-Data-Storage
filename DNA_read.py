import re, random

is_print = True

def pam_finder(dna_seq, pam):
    pam_indices = re.finditer(pam.replace('N', '.'), dna_seq)
    return [pam.end() for pam in pam_indices] # Returns the start position of each pam sequence

def randomizer(success_chance): # recieves odds as a double between 0-1, returns wether the odds were met
    if success_chance >= random.random():
        return True
    else:
        return False
    
def has_mutated(original_ratio, new_ratio, required_ratio_growth):
    if original_ratio * required_ratio_growth <= new_ratio:
        return True
    else:
        return False
    
# TODO: a function that finds the PAMs
# input: sequence + pam
# output: a list of pam indices 

def main(dna_seq, pam):
    pam_indices = pam_finder(dna_seq, pam)
    if is_print:
        print(f'The PAM indices are: {pam_indices}')
    dna_data_storage_proccess.__init__(dna_data_storage_proccess, dna_seq, pam_indices)    
    dna_data_storage_proccess.encode(dna_data_storage_proccess)
    dna_data_storage_proccess.write(dna_data_storage_proccess)
    if is_print:
        print(f'Tagged sequence is: {dna_data_storage_proccess.tagged_dna_seqs[0]}')
    mutated_pams = dna_data_storage_proccess.read(dna_data_storage_proccess)
    if is_print:
        print(f'The lists of mutated pams are {mutated_pams}')
    bit_result = dna_data_storage_proccess.decode(dna_data_storage_proccess, mutated_pams)
    print(f'The resulting bit is: {bit_result}') # not in if is_print because it is the end result


# TODO: joint class for encoder+writer+reader+decoder
# Shared metadata: DNA seq, PAM indices (or at least PAM sequence + logic to generate a list of indices)
# functions:
# encoder: binary message -> "ideal" edited sequence/list of PAM indices to edit (V)
# writer: simulate the enzymatic process. DNA seq + list of PAM indices to edit -> list of *many* edited DNA sequence
#         parameters: edit probability, number of copies, other noise (V)
# reader: simulate NGS. list of *many* edited DNA sequence + DNA seq (unedited) -> list of edited PAM indices
#         parameters: logic for identifying editing from a noisy signal
# decoder: "ideal" edited sequence/list of PAM indices to edit -> binary message

class dna_data_storage_proccess:
    bit_list = [False, True] #array of booleans / 0s or 1s
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
            grna_index_counter = 0
            for index in self.grna_indexes:
                if self.bit_list[grna_index_counter]: # is this bit 1? i.e, does it require editing?
                    self.grna_indexes_for_edit.append(index)
                    for base in range(20):
                        if self.dna_seq[index + base] == 'C':
                            self.tagged_ideal_seq += 'T'
                        else: 
                            self.tagged_ideal_seq += self.dna_seq[index + base]
                grna_index_counter += 1

    def write(self, edit_probability = 1, dna_copy_num = 10):
        is_grna = False
        grna_counter = 0
        for dna_copy in range(dna_copy_num):
            self.tagged_dna_seqs.append('')
            for index, base in enumerate(self.dna_seq):
                if index == self.grna_indexes_for_edit[grna_counter]:
                    is_grna = True
                if is_grna:
                    if base == 'C' and randomizer(edit_probability):
                        self.tagged_dna_seqs[dna_copy] += 'T'
                    else:
                        self.tagged_dna_seqs[dna_copy] += base

                    if index == self.grna_indexes_for_edit[grna_counter] + 19:
                        is_grna = False
                        if not grna_counter == len(self.grna_indexes_for_edit) -1:
                            grna_counter += 1
                else:
                    self.tagged_dna_seqs[dna_copy] += base

    def read(self, read_accuracy = 1, required_ratio_growth = 1): #Doesnt work with only c in original seq 
        tc_ratio = [] # list that contains, for each grna sequence, the UNEDITED Thymine to Cytosine ratio
        for index in self.grna_indexes:
            t_amount = 0
            c_amount = 0
            for base in range(20):
                if self.dna_seq[index + base] == 'C': 
                    c_amount += 1
                elif self.dna_seq[index + base] == 'T': 
                    t_amount += 1
            tc_ratio.append(t_amount/c_amount) #What if no Cs? 
        

        mutated_pams = [] # list of lists, each containing wether each pam site was edited

        for tagged_seq in self.tagged_dna_seqs:
            edited_tc_ratio = [] # list that contains, for each grna sequence, the UNEDITED Thymine to Cytosine ratio
            for index in self.grna_indexes:
                edited_t_amount = 0
                edited_c_amount = 0
                for base in range(20):
                    if self.dna_seq[index + base] == 'C' and randomizer(read_accuracy): #Consult if this is the correct way to randomize reading (needs improvement to determine wether faulty reading results in what base)
                        edited_c_amount += 1
                    elif self.dna_seq[index + base] == 'T' and randomizer(read_accuracy): 
                        edited_t_amount += 1
                edited_tc_ratio.append(t_amount/c_amount)
            mutated_pams.append([])
            for original_ratio, edited_ratio in zip(tc_ratio, edited_tc_ratio):
                mutated_pams[-1].append(has_mutated(original_ratio, edited_ratio, required_ratio_growth)) # Appends True if mutated, False if not

        return mutated_pams

    def decode(self, mutated_pams, mutation_ratio = 0.5):
        pam_mutation = []
        for pam_index in range(len(mutated_pams[0])): # which pam are we going through?
            is_mutated_average = 0
            for pams in mutated_pams: #which seq's pam are we going through?
                is_mutated_average += pams[pam_index]
            is_mutated_average /= len(mutated_pams[0])
            pam_mutation.append(is_mutated_average)
        
        result_bit_list = []
        for mutation_average in pam_mutation:
            result_bit_list.append(mutation_average > mutation_ratio) #If the average mutated pams is larger than the set threshold, append 1, else append 0
        return result_bit_list
        


main('ACCTTGCAGGAATTGAGGCCGTCCGTTAATTTCCCTTGCATACATATTGCGTTTTTTTGTCCTTTTATCCGCTCACTTAGAAAAGAGACAGATAGCTTCT', 'NGG')
