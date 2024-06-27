class required_inputs:
    dna_sequence = 'GTTCCGGCTGTTATCCTGCATCGGAACCTCAATCATGCATCGCACCAGCGTATTCGTGTCATCTAGGAGGGGCGCGTAGGATAAATAATTCAATTAAGATGTCGTTATGCTAGTATACGCCTACCCGTCACCGGCCATCTGTGTGCAGATGGGGCGACGAGTTACTGGCCCTGATTTCTCCGCTTCTAATACCACACACTGGGCAATACGAGCTCAAGCCAGTCTCGCAGTAACGCTCATCAGCTAACGAAAGAGTTAGAGGCTCGCTAATTCGCACTGTCGGGGTCCCTTGGGTGTTTT'
    pam = 'NGG'
    bit_list =  [True, False, True, False, True]

class parameters:
    is_print = False
    edit_probability = 0.8
    read_accuracy = 0.99 #According to oxford's webpage
    copy_nums = 10
    confidence_exponent = 1.1 # The larger the number, the more the ratio of t/tc sum has to grow. 1 would mean the same ratio is marked as mutated.

