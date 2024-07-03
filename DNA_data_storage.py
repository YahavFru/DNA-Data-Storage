import re, random, config


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

class dna_data_storage_process:
    bit_list = [] #array of booleans / 0s or 1s
    dna_seq = '' #original DNA sequence, string
    tagged_ideal_seq = ''
    grna_indices_for_edit = []
    grna_indices = [] #array of indexes for each gRNA block on DNA
    tagged_dna_seqs = [] #array of dna seqs, each edited.
    mutated_pams = [] # list of lists, each containing wether each pam site was edited


    def __init__(self, origin_dna, pam_end_indices, bit_list): #Defines all inputs for the class
        self.bit_list = [] #array of booleans / 0s or 1s
        self.dna_seq = '' #original DNA sequence, string
        self.tagged_ideal_seq = ''
        self.grna_indices_for_edit = []
        self.grna_indices = [] #array of indexes for each gRNA block on DNA
        self.tagged_dna_seqs = [] #array of dna seqs, each edited.
        self.mutated_pams = [] # list of lists, each containing wether each pam site was edited

        self.dna_seq = origin_dna
        self.grna_indices = pam_end_indices
        self.bit_list = bit_list
        
    def encode(self):

        if len(self.bit_list) != len(self.grna_indices): #Bits and bit storage spaces dont match
            print(f'Error: Bit list length ({len(self.bit_list)}) does not match the number of pam sequences present ({len(self.grna_indices)})')
            exit() 
        
        else:
            for list_index, grna_index in enumerate(self.grna_indices):
                if self.bit_list[list_index]:
                    self.grna_indices_for_edit.append(grna_index)
            grna_counter = 0
            for index, base in enumerate(self.dna_seq):
                if index in range(self.grna_indices_for_edit[grna_counter], self.grna_indices_for_edit[grna_counter] + 20) and base == 'C':
                    self.tagged_ideal_seq += 'T'
                else:
                    self.tagged_ideal_seq += base
                # If we passsed the gRNA area and we arent on the last pam
                if index > self.grna_indices_for_edit[grna_counter] + 19 and not grna_counter == len(self.grna_indices_for_edit) - 1: 
                    grna_counter += 1

    def channel(self): #Ideal seq + original seq + parameters for read/write chances -> all edited seqs
        
        edit_probability = config.parameters.edit_probability
        read_accuracy = config.parameters.read_accuracy 
        dna_copy_num = config.parameters.copy_nums

        for tagged_seq in range(dna_copy_num):
            self.tagged_dna_seqs.append('')
            for base, edited_base in zip(self.dna_seq, self.tagged_ideal_seq):
                if base == edited_base: # Same base -> no edit needs to happen, just randomizing sequencing errors
                    if randomizer(read_accuracy):
                        self.tagged_dna_seqs[tagged_seq] += edited_base
                    else:
                        self.tagged_dna_seqs[tagged_seq] += random.choice('ACTG')
                else: # Different bases -> edit needs to happen, randomizes sequencing and edit errors
                    if randomizer(read_accuracy):
                        if randomizer(edit_probability):
                            self.tagged_dna_seqs[tagged_seq] += edited_base
                        else:
                            self.tagged_dna_seqs[tagged_seq] += base
                    else:
                        self.tagged_dna_seqs[tagged_seq] += random.choice('ACTG')


    def decode(self): # Edited seqs, original seq -> bit list

        tc_ratio = [] # list that contains, for each grna sequence, the UNEDITED Thymine to ALL T and C ratio
        for index in self.grna_indices:
            t_amount = self.dna_seq[index:index + 20 if index + 20 <= len(self.dna_seq)-1 else len(self.dna_seq)- 1 - index].count('T')
            c_amount = self.dna_seq[index:index + 20 if index + 20 <= len(self.dna_seq)-1 else len(self.dna_seq)- 1 - index].count('C')
            tc_ratio.append(t_amount/(c_amount + t_amount)) # uses t/t+c to avoid division by 0
        
        edited_tc_ratios = [] # list that contains, for each grna sequence in each dna molecule, the EDITED Thymine to ALL T and C ratio
        for tagged_seq in self.tagged_dna_seqs:
            edited_tc_ratios.append([])
            for index in self.grna_indices:
                edited_t_amount = tagged_seq[index:index + 20 if index + 20 <= len(self.dna_seq)-1 else len(self.dna_seq)- 1 - index].count('T')
                edited_c_amount = tagged_seq[index:index + 20 if index + 20 <= len(self.dna_seq)-1 else len(self.dna_seq)- 1 - index].count('C')
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

        
# Code limitations: RegEx doesnt find overlapping PAMS (AGGG -> AGG (not GGG))