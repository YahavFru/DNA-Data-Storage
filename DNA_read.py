# TODO: a function that finds the PAMs
# input: sequence + pam
# output: a list of pam indices 

def main(dna_seq, pam):
    enumerated_dna_seq = enumerate(dna_seq)
    pam_end_indexes = []
    # TODO: Replace with regex
    pam_match_counter = 0
    for base in enumerated_dna_seq:
        is_pam, pam_match_counter = pam_checker(base[1], pam, pam_match_counter)
        if is_pam:
            pam_match_counter += 1
            if pam_match_counter == len(pam):
                pam_end_indexes.append(base[0] +1)
                pam_match_counter = 0
        else:
            pam_match_counter = 0

    print(f'Pam indexes in DNA: {pam_end_indexes}')
    encoder.__init__(encoder, dna_seq, pam_end_indexes)    
    encoder.write(encoder)

# TODO: Replace with regex
def pam_checker(base, pam, pam_index): #Returns: Is part of PAM? + what part of PAM?
    if pam[pam_index] == 'N' or base == pam[pam_index]:
        return True, pam_index
    
    elif pam_index != 0:
        return pam_checker(base, pam, 0)
    
    return False, -1

# TODO: joint class for encoder+writer+reader+decoder
# Shared metadata: DNA seq, PAM indices (or at least PAM sequence + logic to generate a list of indices)
# functions:
# encoder: binary message -> "ideal" edited sequence/list of PAM indices to edit
# writer: simulate the enzymatic process. DNA seq + list of PAM indices to edit -> list of *many* edited DNA sequence
#         parameters: edit probability, number of copies, other noise
# reader: simulate NGS. list of *many* edited DNA sequence + DNA seq (unedited) -> list of edited PAM indices
#         parameters: logic for identifying editing from a noisy signal
# decoder: "ideal" edited sequence/list of PAM indices to edit -> binary message

class encoder:
    bit_list = [1] #array of booleans / 0s or 1s
    dna_seq = [] #original DNA sequence, array of single char string
    grna_indexes = [] #array of indexes for each gRNA block on DNA
    tagged_dna_seqs = [] #array of dna seqs, each edited.

    def __init__(self, origin_dna, pam_end_indexes): #Defines all inputs for the class
        for i in range(10):
            self.tagged_dna_seqs.append(origin_dna)
        for index in pam_end_indexes:
            self.grna_indexes.append(index)
        if len(self.bit_list) != len(self.grna_indexes):
            print(f'Bit list length ({len(self.bit_list)}) does not match the number of pam sequences present ({len(self.grna_indexes)})')

    def write(self):
        grna_counter = 0
        for dna_seq in self.tagged_dna_seqs:
            for grna in self.grna_indexes:
                if self.bit_list[grna_counter]:
                    for base in range(20):
                        if dna_seq[grna + base] == 'C': 
                            del dna_seq[grna + base]
                            dna_seq.insert(grna + base, 'T') # This is the part where the DNA gets "Tagged", assuming 100% success rate
        print(self.tagged_dna_seqs[9])
main(['A','G','G','T','C','C','C','C','C','C','C','C','G','A','C','C','C','C','C','C','C','C','G'], 'NGG')
