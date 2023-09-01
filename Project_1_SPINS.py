#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import glob
import pandas as pd
import os
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt


# In[ ]:


#NOTE: this has the files removed from the inclusion criteria, the ones from the missing scores, and due to NA values
#NOTE: THIS SETS THE DATA UP FOR FURTHER ANALYSIS 

'''SETTING THINGS UP'''

'''Reads the files'''
#reads the files imported and sorts them numerically 
files_EA_resid_ts =[] # lowkey dont know if this is needed, cause idk if glob automatically puts it into a file or something
files_EA_resid_ts = sorted(glob.glob("/scratch/a/arisvoin/<user>/project_1/data_extracted/NEW_DATA_SPINS/SPINS_Removed_Files/output_exclusion_criteria/output/**/sub-*_RS_2mm_GSR_glasser_tian_meants.csv", recursive = True)) #for SPINS

#getting the data for the first columm here for region names
data = pd.read_csv("G_atlas_info_reorder_correct.csv")
list_of_region_names = data['atlas_roi'].tolist() #takes column that has the brain regions here and puts it within a list

'''Filters the files - to make sure there is no empty csv'''
#filters the files to make sure that there is no empty csv file - when print(filtered_files), you can see that there are no files that are 0 bytes
filtered_files = [file for file in files_EA_resid_ts if os.path.getsize(file) == 0]

'''Creates list of IDs'''
#from Teo's code: gets substrings from the files_EA_resid_ts file and attaches it to ptlist, takes characters 5 to 7, and then 8 to 11
#in python, it takes the starts at the first character mentioned and ends at one before the end of the list
ptlist = []
ptlist = ["SPN01_" + file[126:129] + "_" + file[129:133] for file in files_EA_resid_ts] #85 is the character for the "-" after second "sub" mentioned in csv title 
                                                                                        #file[86:89] + "_" + file[89:93] from using the one without the removed files - for SPINS
          
'''Reads the time series files'''
#creating the read in the time series files
resid_ts = [pd.read_csv(file, header=None) for file in files_EA_resid_ts] #checks out with Teo's code 

'''Names the dictionary with the Participant IDs'''    
#Names dataframe with participant IDs within a dictionary now (participant ID -> time ser|ies data for participant ID) so hopefully that is fine I hope
resid_ts_named = {name: df for name, df in zip(ptlist, resid_ts)}

'''Changes the display limit - to see final matrix'''
pd.options.display.max_rows = 999




'''Testing Puporses'''
#use the first, one dude to test the code below
#test_dude = resid_ts_named['SPN01_CMH_0001'] - for SPINS
#test_dude = resid_ts_named['SPA01_CMP_0001'] - for SPASD
#print(resid_ts_named)

'''Optional: To Transpose Dataframe'''
#transposing dataframes - CAN do this if needed
#transposed_resid_ts = [df.transpose() for df in resid_ts] #in Teo's this is equivalent to "resid_ts" - CHANGE after you check with the varaible equivalents here 
                                                          #but this does transpose the above read_
'''STEPS TO FOLLOW SHOWN BELOW'''
#Z-score here (fisher z score) - "how to fisher z score in ptyhon" in google
#Develop timeserise variablity score here - calculation MSSD thing from gulia?
#Create dataframe of participants by variability score (348 x 392)
#need all of the partipants - need like the dataframe going from 1-209 participants 

#### Running neuroComBat to account for scanner effects - couldnt figure it out - i yeeted it into R and then put the csv back into here after calculating the mssd put path to where the neurocombat code is


# In[ ]:


'''TESTING MENU'''

'''FUNCTIONS'''
def mean_centering(participant):
    center_function = lambda x: x - x.mean()
    centered_data = center_function(participant)
    return centered_data

def z_score(participant):
    data = np.array(participant)
    z_scores = stats.zscore(data, axis=1) #axis = 1 is of each row 
    return z_scores

def mssd_value(participant_scores):
    f_diff = np.diff(participant_scores) 
    f_sc_diff = (f_diff)**2
    mean = np.nanmean(f_sc_diff, axis=1) #axis=1 gets it from the rows
    sq_mean = np.sqrt(mean)
    mssd_values = sq_mean
    return mssd_values

'''MAIN LOOP'''
#initializes variables used later
z_scores_list = []
centered_data_list = []
mssd_values_list = [] #for z-scoring
mssd_second_values = [] #for mean-centering

#creates a numpy array filled with 0's at the start, for row length of "len(resid_ts_named)" and then column length based on the first participant + 1 for the labels
mssd_matrix = np.zeros((len(resid_ts_named), len(resid_ts_named['SPN01_CMH_0001'][0]) + 1), dtype=object)  #have to set dtype because going to attach string to integer - for SPINS
mssd_matrix_mean = np.zeros((len(resid_ts_named), len(resid_ts_named['SPN01_CMH_0001'][0]) + 1), dtype=object) 

while True:
    try:
        #creating user input to see if the data should be mean-centered or z_scored:
        user = int(input("Would you like the data to be z_scored or mean-centered? Enter 1 for z_score and 2 for mean-centered: "))
        if user == 1:
            # this will be part of the main loop
            for i, (key, participant) in enumerate(resid_ts_named.items()): #iterates through all of the participants in resid_ts_named while keepying both (key, value - participant) AND index (i)
                z_scores = (z_score(participant))
                z_scores_list.append(z_scores)
                mssd_values_list.append(mssd_value(z_scores)) #calls functions for calculations

                mssd_values = np.array(mssd_values_list)
                mssd_values_rounded = np.round(mssd_values[-1], decimals=10)  # Round to 10 decimal places
                mssd_matrix[i, 0] = key
                mssd_matrix[i, 1:] = mssd_values_rounded
            
            mssd_matrix = np.vstack(([''] + list_of_region_names, mssd_matrix)) #makes [0][0] an empty space here to accomodate the region name labels, and then vertically stacks the region names atop the matrix
            
            '''MAKES IT LOOK BETTER HERE - bless pandas'''
            df = pd.DataFrame(mssd_matrix)
            break;
            
        if user == 2:
             # this will be part of the main loop
            for i, (key, participant) in enumerate(resid_ts_named.items()): #iterates through all of the participants in resid_ts_named while keepying both (key, value - participant) AND index (i)
                centered_data = (mean_centering(participant))
                centered_data_list.append(centered_data)
                mssd_second_values.append(mssd_value(centered_data)) #calls functions for calculations

                mssd_values_mean = np.array(mssd_second_values)
                mssd_values2_rounded = np.round(mssd_values_mean[-1], decimals=10)  # Round to 10 decimal places
                mssd_matrix_mean[i, 0] = key
                mssd_matrix_mean[i, 1:] = mssd_values2_rounded
            
            
            mssd_matrix_mean = np.vstack(([''] + list_of_region_names, mssd_matrix_mean))
            '''MAKES IT LOOK BETTER HERE - bless pandas'''
            df = pd.DataFrame(mssd_matrix_mean)
            break;
            
        else:
            print("Invalid input. Please enter either 1 for z_score or 2 for mean-centered (integer only): ")
    except ValueError:
        print("Invalid input. Please enter either 1 for z_score or 2 for mean-centered (integer only): ")





print(df)


# In[ ]:


'''NEEDED THIS FOR THE CORRELATION GRAPH AND THE CSV - update: not anymore, cause figured it out'''
#puts the matrix into the shape without the labels
#mssd_matrix_without_labels = np.reshape(mssd_values, (346, 392))
#print(mssd_matrix_without_labels)

'''GETTING CSVs - make sure to check which matrix variable is produced'''
pd.DataFrame(mssd_matrix_mean).to_csv("old_mssd_matrix_mean.csv") #gets a csv file here
#pd.DataFrame(mssd_matrix_without_labels).to_csv("without_labels.csv") 


# In[ ]:


''' THIS SHOWS FOR ALL THE REGIONS AGAINST THE MOTION DATA'''
# Load the region data from CSV file - without the labels
region_data = pd.read_csv('mssd_matrix_mean.csv', header = 1, index_col = 1)

# Load the motion data from CSV file
motion_data = pd.read_csv('NEW_SPINS_final_regressiondata_UPDATED_r_RS_reordered.csv')

# Iterate over each region column
correlations_list = []
for column in region_data.columns:
    region_values = region_data[column].values #this is taken from each column for the regions
    motion_values = motion_data['avg_fd_RS'].values #motion_values is taken from the csv
    
    # Calculate correlation between region values and motion data values
    correlation = np.corrcoef(region_values, motion_values)[0,1]
    correlations_list.append(correlation) #stores the correlations within the list

# Plot a histogram of the correlations
plt.hist(correlations_list, bins=10)
plt.xlabel('Correlation')
plt.ylabel('Frequency')
plt.title('Correlation Histogram')
plt.show()


'''THIS SHOWS FOR EACH REGION SPECIFICALLY - w/o names because its easier that way :p'''
'''
# Load the region data from CSV file
region_data = pd.read_csv('testing.csv', header = 1, index_col = 1)

# Load the motion data from CSV file
motion_data = pd.read_csv('SPINS_final_regressiondata_UPDATED_r_RS_2.csv')

# Get the values of the selected region column
region_values = region_data.iloc[:, 4].values #if you change the value here - you can get each region 
motion_values = motion_data['avg_fd_RS'].values  # Replace 'motion_column' with the actual column name of motion data

# Calculate correlation between region values and motion data values
correlation = np.corrcoef(region_values, motion_values)[0, 1]

# Plot a histogram of the correlations
plt.hist(correlation, bins=10)
plt.xlabel('Correlation')
plt.ylabel('Frequency')
plt.title('Correlation Histogram for Region')
plt.show()
'''


# In[ ]:


'''SENDING MSSD VALUES (NOT THE FIRST 32) INTO TXT FILES FOR EACH PARTICIPANT'''

#for creating some txt files for p-scaler step later 
import csv

def write_rows_to_text(csv_file, folder_name):
    
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)  # Read and store the header row
        
        for i, row in enumerate(reader):
            if i == 0:
                continue  # Skip the first row
            
            # Generate a unique file name for each row
            file_name = f'{row[1]}.txt'
            file_path = os.path.join(folder_name, file_name)  # Create the file path inside the folder
            
            # Write the desired portion of the row to a text file with the corresponding label
            with open(file_name, 'w') as text_file:
                #text_file.write(f'{header[0]}: {row[0]}\n')  # Customize the label format as needed
                text_file.write('\n'.join(row[34:]))  # Skips the first 32 rows
            print(f'Row {i} written to {file_path}')
            
# Usage example
csv_file = 'mssd_matrix_mean.csv'
folder_name = '/project_1/data_extracted/NEW_DATA_SPINS/text_files'
write_rows_to_text(csv_file, "/project_1/data_extracted/NEW_DATA_SPINS/resting_mean_centered_mssd_each_participant")


# In[ ]:


'''CREATING AND SEPARATING THE AVERAGE MSSD CSV's HERE - REARRANGING ROWS BASED ON THE CASE/CONTROL CSV - lowkey has to follow a format look at mssd_matrix_mean.csv'''
# Read the first CSV file for the mssd values
df = pd.read_csv('neuroCombat_mssd.csv')

# Read the second CSV file in cases/controls
df2 = pd.read_csv('NEW_SPINS_final_regressiondata_UPDATED_r_RS_cases_controls.csv')

# Create a dictionary mapping the values in the first column of the second CSV file to their corresponding row index
index_map = {value: index for index, value in enumerate(df2.iloc[:,0])}

# Create a copy of df1 to store the swapped rows
df_swapped = df.copy()

# Iterate through the index_map dictionary
for key, index in index_map.items():
    # Find the row in df1 that matches the key
    matching_row = df[df.iloc[:, 1] == key].values[0]

    # Replace the row at the specified index in df1_swapped
    df_swapped.iloc[index] = matching_row

# Assign the swapped rows back to df1
df1 = df_swapped
# Drop the first column so it leaves out the Unnamed column -> need this for R code
df1_drop = df1.drop(df1.columns[0], axis=1)

# Split the DataFrame into two parts based on row indices
df1_part1 = df1.iloc[:189]  # First n rows for cases
df1_part2 = df1.iloc[189:]  # Remaining rows for controls

'''FOR PART ONE'''
# Convert the columns in df1_part1 to numeric
df1_part1 = df1_part1.apply(pd.to_numeric, errors='coerce')

# Calculate the average for each column in df1_part1
df1_part1_avg = df1_part1.mean()

'''FOR PART TWO'''
# Convert the columns in df1_part1 to numeric
df1_part2 = df1_part2.apply(pd.to_numeric, errors='coerce')

# Calculate the average for each column in df1_part1
df1_part2_avg = df1_part2.mean()

'''average values'''
# Create a new DataFrame for average values
avg_df = pd.DataFrame({'Case Average': df1_part1_avg, 'Control Average': df1_part2_avg})
avg_df_removed = avg_df.iloc[2:] #removes the first two averages

avg_df_trans = avg_df_removed.transpose()


# In[ ]:


pd.DataFrame(df1_drop).to_csv("OLD_neuroCombat_mssd_mean_center_ordered_by_case_control.csv")
#pd.DataFrame(avg_df_trans).to_csv("neuroCombat_average_mean_centered_mssd_values.csv")
#pd.DataFrame(df1_part1).to_csv("controls.csv")


# In[ ]:


'''CREATED CSV AND THEN TXT FILE'''


#for creating some txt files for p-scaler step later 
import csv

def write_rows_to_text(csv_file, folder_name):
    
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)  # Read and store the header row
        
        for i, row in enumerate(reader):
            # Generate a unique file name for each row
            file_name = f'{row[0]}.txt'
            file_path = os.path.join(folder_name, file_name)  # Create the file path inside the folder
            
            # Write the desired portion of the row to a text file with the corresponding label
            with open(file_name, 'w') as text_file:
                #text_file.write(f'{header[0]}: {row[0]}\n')  # Customize the label format as needed
                text_file.write('\n'.join(row[33:]))  # Skips the first 32 rows
            print(f'Row {i} written to {file_path}')
            
# Usage example
csv_file = 'average_mean_centered_mssd_values.csv'
folder_name = '/project_1/data_extracted/NEW_DATA_SPINS/average_mssd'
write_rows_to_text(csv_file, "/project_1/data_extracted/NEW_DATA_SPINS/average_mssd")


# In[ ]:


'''CALCULATING THE T-TEST OF AVERAGE CASES VS CONTROLS - GIVES ONE OVERALL'''

import pandas as pd
from scipy.stats import ttest_ind

# Replace 'your_file.csv' with the path to your actual CSV file
data1 = pd.read_csv('neuroCombat_average_mean_centered_mssd_values.csv', header=None)
data_dropped1 = data1.drop(data1.columns[0], axis=1) #drops the first column
data_dropped2 = data_dropped1.drop(data_dropped1.index[0], axis=0) #drops the first row

# Transpose the DataFrame so that controls and cases are in separate columns
data = data_dropped2.T

# Extract controls and cases data into separate Series
controls = data.iloc[:,1]
cases = data.iloc[:, 0]

# Perform the t-test
t_stat, p_value = ttest_ind(cases, controls)

# Output the results
print("T-Statistic:", t_stat)
print("P-Value:", p_value)


# In[ ]:


'''For T-TEST, FDRCORRECTION AND LOG10 ON THE MSSD MATRIX CSV - GIVES 392 VALUES'''
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection

# Assuming you have two CSV files: cases.csv and controls.csv
# Replace 'cases.csv' and 'controls.csv' with the paths to your actual CSV files

n = 189 #the number of cases in top half of csv according to PYTHON row index
start_row = 191 #the starting point for the bottom half oc csv + 1 (idk python is weird)

# Read the CSV files into DataFrames
controls_data = pd.read_csv('neuroCombat_mssd_mean_center_ordered_by_case_control.csv', header = None, skiprows = start_row -1)
con_data_drop = controls_data.drop(controls_data.columns[0], axis=1)
con_data_drop2 = con_data_drop.drop(con_data_drop.columns[0], axis=1)


cases_data = pd.read_csv('neuroCombat_mssd_mean_center_ordered_by_case_control.csv', nrows = n)
cases_data_drop = cases_data.drop(cases_data.columns[0], axis=1)
cases_data_drop2 = cases_data_drop.drop(cases_data_drop.columns[0], axis=1)


# Assuming the CSV files have 392 columns for each group
num_columns = 392

# Create a list to store the p-values for each column
p_values = []
t_stats = []

# Perform the t-test for each column independently
for column in range(num_columns):
    t_stat, p_value = ttest_ind(cases_data_drop2.iloc[:, column], con_data_drop2.iloc[:, column])
    p_values.append(p_value)
    t_stats.append(t_stat)

# p_values now contains the p-values for each column, representing the significance of the differences
# between cases and controls for each corresponding value within the columns.

# Apply FDR correction to the p-values
#rejected, corrected_p_values, _, _ = multipletests(p_values, method='fdr_tsbh')
rejected, corrected_p_values = fdrcorrection(p_values, alpha=0.05, method='indep')

# Extract statistically significant results (p < 0.05) after FDR correction with their indexes - accounts for the relabelling of the first column as "2"
sig_results = [(i, p_value) for i, p_value in enumerate(corrected_p_values) if p_value < 0.05]

# Apply log-10 transformation
log10_p_values = [(index, -np.log10(p_value)) for index, p_value in sig_results]

index = []
#gets the indexes out of the list to be applied later for grabbing indexes for other variables   
index = [index for index, _ in log10_p_values]

#get the t_stats for the indexes that are significiant
sig_t_stats = [(num, t_stats[num]) for num in index]
#get the t_stats that are negative - THIS ISNT RIGHT
neg_t_stats = [(i, neg_t) for i, neg_t in enumerate(t_stats) if neg_t > 0]


#gets the region names from the indexes of the list
regions = [list_of_region_names[int(num)] for num in index]


'''YEETING THINGS INTO A TXT FILE but gotta 0 regions beforehand'''
list = [0]*392 #intiializing a list with 392 places but filled with 0 

for index, p_value in log10_p_values:
    list[index] = p_value

#removes the first 32 items
final_list = list[32:]


# In[ ]:


'''SENDING LOG10_P_VALUES(NOT THE FIRST 32) INTO TXT FILES FOR EACH PARTICIPANT'''

#for creating some txt files for p-scaler step later 
import csv

# Assuming you have a list named 'my_list' with 360 items
# Replace 'my_list' with the actual name of your list

# Specify the file path for the new text file
file_path = 'log10_p_values_list.txt'

# Open the file in write mode and write the elements of the list to the file
with open(file_path, 'w') as file:
    for item in final_list:
        file.write(str(item) + '\n')


# In[ ]:


''' COPIED HERE FROM PREVIOUS CODE IF NEEDED - THIS SHOWS FOR ALL THE REGIONS AGAINST THE MOTION DATA'''
'''
# Load the region data from CSV file - without the labels
region_data = pd.read_csv('testing.csv', header = 1, index_col = 1)

# Load the motion data from CSV file
motion_data = pd.read_csv('SPINS_final_regressiondata_UPDATED_r_RS_reordered.csv')

# Iterate over each region column
correlations_list = []
for column in region_data.columns:
    region_values = region_data[column].values #this is taken from each column for the regions
    motion_values = motion_data['avg_fd_RS'].values #motion_values is taken from the csv
    
    # Calculate correlation between region values and motion data values
    correlation = np.corrcoef(region_values, motion_values)[0,1]
    correlations_list.append(correlation) #stores the correlations within the list

# Plot a histogram of the correlations
plt.hist(correlations_list, bins=10)
plt.xlabel('Correlation')
plt.ylabel('Frequency')
plt.title('Correlation Histogram')
plt.show()
'''

'''THIS SHOWS FOR EACH REGION SPECIFICALLY - w/o names because its easier that way :p'''
'''
# Load the region data from CSV file
region_data = pd.read_csv('testing.csv', header = 1, index_col = 1)

# Load the motion data from CSV file
motion_data = pd.read_csv('SPINS_final_regressiondata_UPDATED_r_RS_2.csv')

# Get the values of the selected region column
region_values = region_data.iloc[:, 4].values #if you change the value here - you can get each region 
motion_values = motion_data['avg_fd_RS'].values  # Replace 'motion_column' with the actual column name of motion data

# Calculate correlation between region values and motion data values
correlation = np.corrcoef(region_values, motion_values)[0, 1]

# Plot a histogram of the correlations
plt.hist(correlation, bins=10)
plt.xlabel('Correlation')
plt.ylabel('Frequency')
plt.title('Correlation Histogram for Region')
plt.show()
'''

