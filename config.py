class required_inputs:
    dna_sequence = 'GGGACTGACTGACTGACTGACTGAGGGTCAGTCAGTCAGTCAGTCA'
    pam = 'NGG'
    bit_list =  1 * [True, False, True, False, True]

class parameters:
    is_print = False
    edit_probability = 1
    read_accuracy = 1 #According to oxford's webpage
    copy_nums = 100
    confidence_exponent = 1.1 # The larger the number, the more the ratio of t/tc sum has to grow. 1 would mean the same ratio is marked as mutated
    offtarget_exponent = 0.8