def main(dna_seq, pam):
    pam = list(pam)
    print(pam)
    dna_seq = enumerate(dna_seq)
    pam_start_indexes = []
    pam_match_counter = 0
    for base in dna_seq:
        is_pam, pam_match_counter = pam_checker(base[1], pam, pam_match_counter)
        if is_pam:
            pam_match_counter += 1
            if pam_match_counter == len(pam):
                pam_start_indexes.append(base[0] - len(pam) + 1)
                pam_match_counter = 0

        else:
            pam_match_counter = 0
    print(pam_start_indexes)

def pam_checker(base, pam, pam_index): #Returns: Is part of PAM? + what part of PAM?
    if pam[pam_index] == 'N' or base == pam[pam_index]:
        return True, pam_index
    
    elif pam_index != 0:
        return pam_checker(base, pam, 0)
    
    return False, -1
    

main(['G', 'C', 'G', 'A', 'A', 'A', 'A', 'T', 'G', 'G', 'T', 'C', 'C', 'T', 'T', 'A', 'A', 'G', 'A', 'A', 'T', 'C', 'A', 'T', 'C', 'A', 'T', 'C', 'G', 'G', 'G', 'T', 'T', 'C', 'T', 'T', 'T', 'T', 'A', 'A', 'G', 'A', 'C', 'G', 'A', 'G', 'G', 'C', 'C', 'A'], 'NGG')
