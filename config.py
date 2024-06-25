class required_inputs:
    dna_sequence = 'CTCAGCTCTATTTTAGTGGTCATGGGTTTTGGTCCGCCCGAGCGGTGCAACCGATTAGGACCATGTAAAACATTTGTTACAAGTCTTCTTTTAAACACAATCTTCCTGCTCAGTGGCGCATGATTATCGTTGTTGCTAGCCAGCGTGGTAAGTAACAGCACCACTGCGAGCCTAATGTGCCCTTTCCACGAACACAGGGCTGTCCGATCCTATATTAGGACTCCGCAATGGGGTTAGCAAGTCGCACCCTAAACGATGTTGAAGACTCGCGATGTACATGCTCTGATACAATACATACGT'
    pam = 'NGG'
    bit_list = [True, True, True, True, False, False, True, False, True, False]

class parameters:
    is_print = True
    edit_probability = 1
    read_accuracy = 1
    copy_nums = 10
    confidence_exponent = 1.1 # The larger the number, the more the ratio of t/tc sum has to grow. 1 would mean the same ratio is marked as mutated.