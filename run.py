import DNA_data_storage as main
import config
import pandas as pd
import matplotlib.pyplot as plt

df = pd.DataFrame(columns=['Confidence Exponent', 'Result'])
for i in range(10):
    config.parameters.confidence_exponent += i * 0.05
    result_bits = main.main(config.required_inputs.dna_sequence, config.required_inputs.pam, config.required_inputs.bit_list)
    correct_bits = 0
    for result_bit, bit in zip(result_bits, config.required_inputs.bit_list):
        correct_bits += result_bit == bit
    correct_bits /= len(config.required_inputs.bit_list)
    new_row = {'Confidence Exponent': config.parameters.confidence_exponent, 'Result': correct_bits}
    df = df._append(new_row, ignore_index=True)

fig, ax = plt.subplots()
ax.plot(df['Confidence Exponent'], df['Result'], marker='o', linestyle='-', color='blue')
ax.set_xlabel('Confidence Exponent')
ax.set_ylabel('Result')
ax.set_title('Input vs Output')
print(df)
plt.show()

