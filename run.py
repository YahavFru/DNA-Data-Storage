import DNA_data_storage as main
import config
import pandas as pd
import random
import plotly.express as px
import numpy as np
import plotly.graph_objects as go

def clean_dna_generator(count):
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
            
        dna_seq += ''.join(random.choices('ACTG', k= random.randint(20, 39)))
        
        pam_indices = main.pam_finder(dna_seq, pam)
        if len(pam_indices) == count and len(dna_seq) >= pam_indices[-1] + 19:
            return dna_seq
        elif len(pam_indices) == counter:
            counter += 1
            dna_checkpoint = dna_seq
        
##############################################################################################

def dna_generator(count):
    pam = config.required_inputs.pam
    dna_checkpoint = ''
    while True:
        dna_seq = dna_checkpoint

        for base in pam:
            if base == 'N':
                dna_seq += random.choice('ACTG')
            else:
                dna_seq += base
            
        dna_seq += ''.join(random.choices('ACTG', k= random.randint(20, 39)))
        
        pam_indices = main.pam_finder(dna_seq, pam)
        if len(pam_indices) == count and len(dna_seq) >= pam_indices[-1] + 19:
            return dna_seq
        elif len(pam_indices) > count:
            dna_checkpoint = ''
        elif len(dna_seq) >= pam_indices[-1] + 19:
            dna_checkpoint = dna_seq        

########################################################################################

def conf_exponent_graph():
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
        avg_sum = 0
        for each in single_seq_results:
            std.append(each[seq_i_num])
            avg_sum += each[seq_i_num]
        avg_results.append(avg_sum/len(single_seq_results))
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


#########################################################################################

def edit_probability_graph():
    # Initialize lists for the sequences
    copy_nums_list = [100, 25, 1]  # Different copy numbers to plot
    edit_rates = np.arange(1, 0.0, -0.05).tolist()  # Edit rates from 0.8 to 0.05
    num_repeats = 100  # Number of repetitions for averaging

    avg_results_dict = {}
    std_results_dict = {}

    for copy_nums in copy_nums_list:
        all_avg_results = []
        all_std_results = []

        for edit_rate in edit_rates:
            all_results = []

            for _ in range(num_repeats):
                config.parameters.copy_nums = copy_nums
                config.parameters.edit_probability = edit_rate
                config.required_inputs.dna_sequence = dna_generator(5)
                result_bits = main.main(config.required_inputs.dna_sequence,
                                        config.required_inputs.pam,
                                        config.required_inputs.bit_list)
                correct_bits = sum(result_bit == bit for result_bit, bit in zip(result_bits, config.required_inputs.bit_list))
                accuracy = correct_bits / len(config.required_inputs.bit_list)
                all_results.append(accuracy)

            avg_result = np.mean(all_results)
            std_result = np.std(all_results)

            all_avg_results.append(avg_result)
            all_std_results.append(std_result)

        avg_results_dict[copy_nums] = all_avg_results
        std_results_dict[copy_nums] = all_std_results

    # Create dataframes for plotting
    data_frames = []
    for copy_nums in copy_nums_list:
        df = pd.DataFrame({
            'Edit Rate': edit_rates,
            'Average Result': avg_results_dict[copy_nums],
            'Accuracy Std': std_results_dict[copy_nums]
        })
        data_frames.append(df)

    # Create the line plot with error bars
    fig = go.Figure()

    colors = ['#00bcd4', '#ffc107', '#757575']
    for i, copy_nums in enumerate(copy_nums_list):
        fig.add_trace(go.Scatter(
            x=data_frames[i]['Edit Rate'],
            y=data_frames[i]['Average Result'],
            mode='lines+markers',
            name=f'{copy_nums} Copies',
            line=dict(color=colors[i]),
            error_y=dict(type='data', array=data_frames[i]['Accuracy Std'], color=colors[i])
        ))

    fig.update_layout(
        title='Edit Probability vs Accuracy for Different Copy Numbers',
        xaxis_title='Edit Probability',
        yaxis_title='Accuracy',
        # template='plotly_dark'
    )

    fig.show()

# Example usage
# edit_probability_graph()

#############################################################################

def copy_num_graph(num_seq_1, num_seq_2, num_seq_3):
    # Initialize lists for the sequences
    copy_nums = []
    avg_results_1 = []
    avg_results_2 = []
    avg_results_3 = []
    single_seq_results_1 = []
    single_seq_results_2 = []
    single_seq_results_3 = []
    accuracy_std_1 = []  # Standard deviation for the first number of sequences
    accuracy_std_2 = []  # Standard deviation for the second number of sequences
    accuracy_std_3 = []  # Standard deviation for the third number of sequences

    # Process sequences for num_seq_1
    for dna_seq_num in range(num_seq_1):
        config.required_inputs.dna_sequence = clean_dna_generator(5)
        single_seq_results_1.append([])
        copy_nums = []  # Same for each seq
        config.parameters.copy_nums = 1

        for i in range(20):
            config.parameters.copy_nums = 1
            config.parameters.copy_nums += 2 * i
            result_bits = main.main(config.required_inputs.dna_sequence,
                                    config.required_inputs.pam,
                                    config.required_inputs.bit_list)
            correct_bits = sum(result_bit == bit for result_bit, bit in zip(result_bits, config.required_inputs.bit_list))
            accuracy = correct_bits / len(config.required_inputs.bit_list)

            copy_nums.append(config.parameters.copy_nums)
            single_seq_results_1[-1].append(accuracy)

    # Process sequences for num_seq_2
    for dna_seq_num in range(num_seq_2):
        config.required_inputs.dna_sequence = clean_dna_generator(5)
        single_seq_results_2.append([])
        copy_nums = []  # Same for each seq
        config.parameters.copy_nums = 1

        for i in range(20):
            config.parameters.copy_nums = 1
            config.parameters.copy_nums += 2 * i
            result_bits = main.main(config.required_inputs.dna_sequence,
                                    config.required_inputs.pam,
                                    config.required_inputs.bit_list)
            correct_bits = sum(result_bit == bit for result_bit, bit in zip(result_bits, config.required_inputs.bit_list))
            accuracy = correct_bits / len(config.required_inputs.bit_list)

            copy_nums.append(config.parameters.copy_nums)
            single_seq_results_2[-1].append(accuracy)
        
    copy_nums = []
    # Process sequences for num_seq_3
    for dna_seq_num in range(num_seq_3):
        config.required_inputs.dna_sequence = clean_dna_generator(5)
        single_seq_results_3.append([])
        copy_nums = []  # Same for each seq
        config.parameters.copy_nums = 1

        for i in range(20):
            config.parameters.copy_nums = 1
            config.parameters.copy_nums += 2 * i
            result_bits = main.main(config.required_inputs.dna_sequence,
                                    config.required_inputs.pam,
                                    config.required_inputs.bit_list)
            correct_bits = sum(result_bit == bit for result_bit, bit in zip(result_bits, config.required_inputs.bit_list))
            accuracy = correct_bits / len(config.required_inputs.bit_list)

            copy_nums.append(config.parameters.copy_nums)
            single_seq_results_3[-1].append(accuracy)

    def calculate_avg_and_std(single_seq_results):
        avg_results = []
        accuracy_std = []
        for seq_i_num in range(len(single_seq_results[0])):
            std = []
            avg_sum = 0
            for each in single_seq_results:
                std.append(each[seq_i_num])
                avg_sum += each[seq_i_num]
            avg_results.append(avg_sum / len(single_seq_results))
            accuracy_std.append(np.std(std))
        return avg_results, accuracy_std

    avg_results_1, accuracy_std_1 = calculate_avg_and_std(single_seq_results_1)
    avg_results_2, accuracy_std_2 = calculate_avg_and_std(single_seq_results_2)
    avg_results_3, accuracy_std_3 = calculate_avg_and_std(single_seq_results_3)
    df_1 = pd.DataFrame({
        'Number of copies': copy_nums,
        'Average Result': avg_results_1,
        'Accuracy Std': accuracy_std_1
    })

    df_2 = pd.DataFrame({
        'Number of copies': copy_nums,
        'Average Result': avg_results_2,
        'Accuracy Std': accuracy_std_2
    })

    df_3 = pd.DataFrame({
        'Number of copies': copy_nums,
        'Average Result': avg_results_3,
        'Accuracy Std': accuracy_std_3
    })

    # Create the line plot with error bars
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df_1['Number of copies'],
        y=df_1['Average Result'],
        mode='lines+markers',
        name=f'{num_seq_1} Sequences',
        line=dict(color = '#00bcd4'),
        error_y=dict(type='data', array=df_1['Accuracy Std'],color = '#00bcd4')
    ))

    fig.add_trace(go.Scatter(
        x=df_2['Number of copies'],
        y=df_2['Average Result'],
        mode='lines+markers',
        name=f'{num_seq_2} Sequences',
        line=dict(color = '#ffc107'),
        error_y=dict(type='data', array=df_2['Accuracy Std'], color = '#ffc107')
    ))

    fig.add_trace(go.Scatter(
        x=df_3['Number of copies'],
        y=df_3['Average Result'],
        mode='lines+markers',
        name=f'{num_seq_3} Sequences',
        line=dict(color = '#757575'),
        error_y=dict(type='data', array=df_3['Accuracy Std'], color = '#757575')   
    ))

    fig.update_layout(
        title='Accuracy as a function of DNA copy number',
        xaxis_title='Copy number',
        yaxis_title='Accuracy',
        # template='plotly_dark'
    )

    fig.show()
# copy_num_graph(1000, 100, 10)

###################################################

def copy_num_graph():
    # Initialize lists for the sequences
    copy_nums = [0] + np.arange(1, 20, 1).tolist()  # Start with 0
    edit_rates = [0.8, 0.4, 0.1]
    num_repeats = 150
    show_std = False
    # Function to calculate average and standard deviation
    def calculate_avg_and_std(single_seq_results):
        avg_results = np.mean(single_seq_results, axis=0)
        accuracy_std = np.std(single_seq_results, axis=0)
        return avg_results, accuracy_std

    # Collect results for each edit rate
    avg_results_dict = {}
    std_results_dict = {}

    for edit_rate in edit_rates:
        config.parameters.edit_probability = edit_rate
        all_results = []

        for _ in range(num_repeats):
            single_seq_results = [0]  # Start with accuracy 0 for 0 copy nums
            for copy_num in copy_nums[1:]:  # Skip the first element (0)
                config.required_inputs.dna_sequence = clean_dna_generator(5)
                config.parameters.copy_nums = copy_num
                result_bits = main.main(config.required_inputs.dna_sequence,
                                        config.required_inputs.pam,
                                        config.required_inputs.bit_list)
                correct_bits = sum(result_bit == bit for result_bit, bit in zip(result_bits, config.required_inputs.bit_list))
                accuracy = correct_bits / len(config.required_inputs.bit_list)
                single_seq_results.append(accuracy)
            all_results.append(single_seq_results)

        avg_results, accuracy_std = calculate_avg_and_std(np.array(all_results))
        avg_results_dict[edit_rate] = avg_results
        std_results_dict[edit_rate] = accuracy_std

    # Create dataframes for plotting
    df_0_8 = pd.DataFrame({
        'Number of copies': copy_nums,
        'Average Result': avg_results_dict[edit_rates[0]],
        'Accuracy Std': std_results_dict[edit_rates[0]]
    })

    df_0_4 = pd.DataFrame({
        'Number of copies': copy_nums,
        'Average Result': avg_results_dict[edit_rates[1]],
        'Accuracy Std': std_results_dict[edit_rates[1]]
    })

    df_0_1 = pd.DataFrame({
        'Number of copies': copy_nums,
        'Average Result': avg_results_dict[edit_rates[2]],
        'Accuracy Std': std_results_dict[edit_rates[2]]
    })

    # Create the line plot with error bars
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df_0_8['Number of copies'],
        y=df_0_8['Average Result'],
        mode='lines+markers',
        name='Edit rate 0.8',
        line=dict(color='#00bcd4'),
        # error_y=dict(type='data', array=df_0_8['Accuracy Std'], color='#00bcd4')
    ))

    fig.add_trace(go.Scatter(
        x=df_0_4['Number of copies'],
        y=df_0_4['Average Result'],
        mode='lines+markers',
        name='Edit rate 0.4',
        line=dict(color='#ffc107'),
        # error_y=dict(type='data', array=df_0_4['Accuracy Std'], color='#ffc107')
    ))

    fig.add_trace(go.Scatter(
        x=df_0_1['Number of copies'],
        y=df_0_1['Average Result'],
        mode='lines+markers',
        name='Edit rate 0.1',
        line=dict(color='#757575'),
        # error_y=dict(type='data', array=df_0_1['Accuracy Std'], color='#757575')
    ))

    fig.update_layout(
        title='Accuracy as a function of DNA copy number',
        xaxis_title='Copy number',
        yaxis_title='Accuracy',
        # template='plotly_dark'
    )

    fig.show()

# copy_num_graph()

##################################################

print(main.main(clean_dna_generator(5), config.required_inputs.pam, config.required_inputs.bit_list))