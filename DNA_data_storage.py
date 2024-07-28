import re, random, config, numpy as np
from scipy.spatial.distance import pdist, squareform


def pam_finder(dna_seq, pam): # Dna + pam -> pam end indices
    pam_indices = re.finditer(pam.replace('N', '.'), dna_seq) #RegEx ignores overlap giving only the first one
    return [pam.end() for pam in pam_indices] # Returns the start position of each pam sequence


def randomizer(success_chance): # recieves odds as a double between 0-1, returns wether the odds were met
    return success_chance >= random.random()

def has_mutated(original_ratios, new_ratios): # List of each protospacer's t/t+c original ratios + list of each tagged seq's protospacer's t/t+c edited ratios -> Bit list
    min_confidence_exponent = config.parameters.confidence_exponent
    sum_mutations = []
    for ratio in original_ratios:
        sum_mutations.append(0)
    for ratios in new_ratios:
        for ratio_index, edited_ratio in enumerate(ratios):
            sum_mutations[ratio_index] += edited_ratio > original_ratios[ratio_index] * min_confidence_exponent

    mutated_bits = []
    for sum_mut in sum_mutations:
        mutated_bits.append(sum_mut / len(new_ratios) >= 0.5)

    return mutated_bits

def hamming_distance_matrix(dna_seq, prtspcr_index):
    base_dict ={
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }
    protospacers = [[base_dict[base] for base in dna_seq[index:index+20]] for index in prtspcr_index]
    arr = np.array([list(s) for s in protospacers])
    distances = pdist(arr, metric='hamming')
    distance_matrix = squareform(distances)
    return distance_matrix

def protospacer_editor(prtspcr_start_index, all_sequences, sequence_edit_index): #protospacer start index + all sequences + sequences to edit -> ALL edited protospacers.
    sequence_edit_index.sort()
    index_counter = 0
    for seq_index, sequence in enumerate(all_sequences):
        if index_counter >= len(sequence_edit_index):
            break
        if seq_index == sequence_edit_index[index_counter]:
            index_counter += 1
            all_sequences[seq_index] = all_sequences[seq_index][:prtspcr_start_index] + ''.join(["T" if base == "C" and randomizer (config.parameters.edit_probability) else base for base in all_sequences[seq_index][prtspcr_start_index:prtspcr_start_index+20]]) + all_sequences[seq_index][prtspcr_start_index+20:]
    return all_sequences


class dna_data_storage_process:

    def __init__(self, origin_dna, pam_end_indices, bit_list): #Defines all inputs for the class + Resets for multiple back to back runs
        #Settings for class:
        self.dna_seq = origin_dna #original DNA sequence, string
        self.prtspcr_indices = pam_end_indices #array of indexes for each protospacer block on DNA
        self.bit_list = bit_list #array of booleans
        
    def encode(self): #Bit list, Protospacer indices, DNA sequence, -> Ideal edited seq

        if len(self.bit_list) != len(self.prtspcr_indices): #Bits and bit storage spaces dont match
            print(f'Error: Bit list length ({len(self.bit_list)}) does not match the number of pam sequences present ({len(self.prtspcr_indices)})')
            exit() 
        
    
        protspcr_indices_for_edit = [prtspcr_index for list_index, prtspcr_index in enumerate(self.prtspcr_indices) if self.bit_list[list_index]] #For each protospacer, check wether it is a 1, ie does it need editing.
        self.tagged_ideal_seq = ''
        
        prtspcr_counter = 0
        for index, base in enumerate(self.dna_seq):
            if index in range(protspcr_indices_for_edit[prtspcr_counter], protspcr_indices_for_edit[prtspcr_counter] + 20) and base == 'C':
                self.tagged_ideal_seq += 'T'
            else:
                self.tagged_ideal_seq += base
            # If we passsed the prtspcr area and we arent on the last pam
            if index > protspcr_indices_for_edit[prtspcr_counter] + 19 and not prtspcr_counter == len(protspcr_indices_for_edit) - 1: 
                prtspcr_counter += 1

    

    def channel(self): #Ideal seq + original seq + parameters for read/write chances -> all edited seqs
        
        edit_probability = config.parameters.edit_probability
        read_accuracy = config.parameters.read_accuracy 
        dna_copy_num = config.parameters.copy_nums
        self.tagged_dna_seqs = dna_copy_num * [self.dna_seq] #Creates n copies of DNA

        dist_matrix = hamming_distance_matrix(self.dna_seq, self.prtspcr_indices) #n x n matrix of similarity between protospacers
        print(dist_matrix)
        for prtspcr_num in range(len(self.prtspcr_indices)): #For each COLUMN in dist matrix
            for second_prtspcr_num in range(len(self.prtspcr_indices)): #For each ROW in dist matrix
                if self.bit_list[prtspcr_num]:
                    copy_edit_amount = np.random.binomial(dna_copy_num, (1-dist_matrix[prtspcr_num][second_prtspcr_num]) * config.parameters.offtarget_exponent) # How many protospacers to edit
                    seq_edit_indices = random.sample(range(dna_copy_num), copy_edit_amount) #Which ones are edited (Randomly selected)
                    protospacer_editor(self.prtspcr_indices[second_prtspcr_num], self.tagged_dna_seqs, seq_edit_indices)

        

    def decode(self): # Edited seqs, original seq, prtspcr_indices (Can be replaced) -> bit list
                 
        tc_ratio = [] # list that contains, for each prtspcr sequence, the UNEDITED Thymine to ALL T and C ratio
        for index in self.prtspcr_indices:
            end_index = index + 20 if index + 20 <= len(self.dna_seq)-1 else len(self.dna_seq) -1
            t_amount = self.dna_seq[index:end_index].count('T')
            c_amount = self.dna_seq[index:end_index].count('C')
            tc_ratio.append(t_amount/(c_amount + t_amount)) # uses t/t+c to avoid division by 0
        
        edited_tc_ratios = [] # list that contains, for each prtspcr sequence in each dna molecule, the EDITED Thymine to ALL T and C ratio
        for tagged_seq in self.tagged_dna_seqs:
            edited_tc_ratios.append([])
            for index in self.prtspcr_indices:
                end_index = index + 20 if index + 20 <= len(tagged_seq)-1 else len(tagged_seq) -1
                edited_t_amount = tagged_seq[index:end_index].count('T')
                edited_c_amount = tagged_seq[index:end_index].count('C')
                edited_tc_ratios[-1].append(edited_t_amount/(edited_c_amount + edited_t_amount)) # uses t/t+c to avoid division by 0
        return has_mutated(tc_ratio, edited_tc_ratios)
        
def main(dna_seq, pam, bit_list):

    pam_indices = pam_finder(dna_seq, pam) 
    dna_data_storage_process.__init__(dna_data_storage_process, dna_seq, pam_indices, bit_list)    
    dna_data_storage_process.encode(dna_data_storage_process)
    dna_data_storage_process.channel(dna_data_storage_process)
    bit_result = dna_data_storage_process.decode(dna_data_storage_process)
    
    if config.parameters.is_print:
        print(f'The resulting bit is: {bit_result}')
        print(f'Tagged {len(dna_data_storage_process.tagged_dna_seqs)} DNA sequences')

    return bit_result
