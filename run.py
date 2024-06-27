import DNA_data_storage as main
import config
import pandas as pd
import random
import plotly.express as px
import numpy as np

def dna_generator(count):
    pam = config.required_inputs.pam
    dna_checkpoint = ''
    counter = 1
    while True:
        dna_seq = dna_checkpoint

        for base in pam:
            if base == 'N':
                dna_seq += random.choice('ACTG')
            else:
                dna_seq += base
            
        dna_seq += ''.join(random.choices('ACTG', k= random.randint(20, 20))) #20, 39
        
        pam_indices = main.pam_finder(dna_seq, pam)
        if len(pam_indices) == count and len(dna_seq) >= pam_indices[-1] + 19:
            return dna_seq
        elif len(pam_indices) == counter:
            counter += 1
            dna_checkpoint = dna_seq
        

confidence_exponents = []
avg_results = []
single_seq_results = []
accuracy_std = [] #Standard deviation

for dna_seq_num in range(10):
    config.required_inputs.dna_sequence = dna_generator(5)
    single_seq_results.append([])
    confidence_exponents = [] #Same for each seq
    config.parameters.confidence_exponent = 1.1

    for i in range(10):
        config.parameters.confidence_exponent += i * 0.05
        print(config.parameters.confidence_exponent)
        result_bits = main.main(config.required_inputs.dna_sequence,
                                config.required_inputs.pam,
                                config.required_inputs.bit_list)
        
        correct_bits = sum(result_bit == bit for result_bit, bit in zip(result_bits, config.required_inputs.bit_list))
        accuracy = correct_bits / len(config.required_inputs.bit_list)
        
        confidence_exponents.append(config.parameters.confidence_exponent)
        single_seq_results[-1].append(accuracy)

for seq_i_num in range(len(single_seq_results[0])):
    std = []
    sum = 0
    for each in single_seq_results:
        std.append(each[seq_i_num])
        sum += each[seq_i_num]
    avg_results.append(sum/len(single_seq_results))
    accuracy_std.append(np.std(std)) #Compute and append std




df = pd.DataFrame({
    'Confidence Exponent': confidence_exponents,
    'Average Result': avg_results,
    'Accuracy Std': accuracy_std
})

# Create the line plot with error bars
fig = px.line(df, x='Confidence Exponent', y='Average Result', title='Confidence Exponents vs Average Result with Std Dev')

# Add error bars
fig.update_traces(error_y=dict(type='data', array=df['Accuracy Std']))

# Show the plot
fig.show()