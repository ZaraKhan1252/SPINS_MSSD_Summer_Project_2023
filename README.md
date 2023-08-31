# Summer2023
Summer Project 2023 within the KIMEL Lab under the co-supervision of Teo Secara and Dr. Colin Hawco. This lab was absolutely epic and I learned a TON of stuff about neurosci, coding, statistics and applications of research throughout my first research placement! :)

---
title: "README"
output: html_document

DISCLAIMER:
This will follow the SPINS dataset, but the code and steps are equivalent for the SPASD Dataset. 

Who are you?
# Hello future gamers and gamees, my name is Zara Khan and I was a summer student working in the KIMEL Lab under the co-supervision of Teo Secara and Dr.Colin Hawco. I was part of "The Piposaurs" which included Clara Sun and Fariah Sandhu, since we started in May and attened Brainhack School, instead of starting in June like "Group 2". Yes, we called them Group 2. I will be going into my third year of Electrical and Biomedical Engineering at McMaster University coming this fall but idk where I'll be when you are reading this hehe. Hopefully somewhere nice (please) mentally and physically (hopefully not in delulu land). Anyway, over the four months of my summer in 2023, this README shall document the more important things I did - instead of the other things as in exploring CAMH, Toronto, and making some epic summer friends from the other summer students researching at the same time as me hehe. Hoping that this isn't confusing because man, oh man, did I have to yeet around all my files for this (it was so disorganized, idk how I found anything in my files), so I hope, wish and pray that you'll be able to find everything that you're ever looking for within my projects folder - in a reasonable context, that is.
#
# Oh! And I also QC'd for the SPINS and SPASD database that was used here, along with Fariah Sandhu - and then also QC'd for OPTIMUM - but like I have no idea what's going on with that dataset, I just QC'd. It was the death of me, until I realized that my code was the death of me when it didn't choose to work - AND OMG CITRIX PLEASE THAT GODFORSAKEN PROGRAM YEETED MY LAPTOP INTO BLUE SCREEN TERRITORRY A HANDFUL OF TIMES- anyway

What did you do
#SO basically my entire code is lowkey on SciNet's Jupyter Notebook, I'm a python gworlie. But I also used R for PLSC, getting loadings, etc. shall be explained later. This shall be as chronoligcal as the homosapiens timeline.

#EXTRACTING TIMESERIES DATA

# Code can found in /scratch/zkhan/SPINS_Timeseries_Extraction/extract_timr_series_RS_2mm_GSR_glasser_tian_SPINS.sh

# 1) Descending into thy folder listed above - this code when run in terminator, bash, control panel, etc. allows you to get the timeseries data from preprocessed data (bless Ju-Chi). There are comments that explain what each line means, but basically it pulls each subject's data from sub_dir and then either makes a folder if there isn't a folder made for that subject already and then computes the timeseries extraction through "ciftify_meants" (a function that Erin coded, very big brain). The directory that you put after "--outputcsv" is where the extracted timeseries csv files will be placed into. 

# Mine is found at /scratch/zkhan/SPINS_Timeseries_Extraction/output. 

#To run in terminal: "sbatch [name of entire path - mine would be /scratch/zkhan/SPINS_Timeseries_Extraction/extract_timr_series_RS_2mm_GSR_glasser_tian_SPINS.sh ] if not in the actual directory or just sbatch[name of the .sh file] if in the actual directory

#COPYING TO JUPYTER NOTEBOOK

# 2) You'll know this command by the end of this, but you have to copy the entire folder to jupyter notebook in order to access it. All files that are used with the code in Jupyter Notebook are on Jupyter Notebook. But in order to do this, you need to make sure that the desktop that you are using is SSH'd with your SciNet account (ask Dawn or Kevin about this). If it is not, then the command below will not work.

# Code to be run in terminal: 
# scp -r [the file path of the files you want to copy over] [scinetusername]@niagara.scinet.utoronto.ca:[file path of where you want the copied files to be stored into on Jupyter Notebook]
# An example:
# scp - r /scratch/zkhan/SPINS_Timeseries_Extraction/extract_timr_series_RS_2mm_GSR_glasser_tian_SPINS.sh zara1252@niagara.scinet.utoronto.ca:/scratch/project_1/data_extracted/NEW_DATA_SPINS/SPINS_Removed_Files/output

#ONCE IN JUPYTER NOTEBOOK

#Code can be found in: projects/zkhan/proj1/JupyterNotebook_Code/project_1/data_extracted/NEW_DATA_SPINS/Project_1_SPINS.ipynb
#NOTE: Everytime you see a "In[ ]" within the python file, it is another executable cell within Jupyter Notebook. You need to run the cells chronologically from top to bottom so that all variables are executed and can be run within the other cells. 

# 3) I am running this within the environment [conda env:.conda-3.6_nilearn_01], but like I never used nilearn within the code. SO. Don't ask - okay fine, it was for BrainHacks School and I never switched it and here we are. I think you should still be able to run the code within the standard environment. 

#Don't mind the red lines - they just show here in R, but is fine in Jupyter Notebook - and hopefully will help with understanding what code block I am talking about

# The first cell block is 
```{r setup, include=FALSE}
import glob
import pandas as pd
import os
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
```

# which include all the libraries used within the code. This should just download the libraries and be ready to use as you go on. You would need to run this cell block everytime you restart the kernel.


# The second cell block is 
```{r setup, include=FALSE}
#NOTE: this has the files removed from the inclusion criteria, the ones from the missing scores, and due to NA values
#NOTE: THIS SETS THE DATA UP FOR FURTHER ANALYSIS 

'''SETTING THINGS UP'''

'''Reads the files'''
#reads the files imported and sorts them numerically 
files_EA_resid_ts =[] # lowkey dont know if this is needed, cause idk if glob automatically puts it into a file or something - has 2 participants taken out MRC_166 because
                      # they didnt complete the study and then MRC_0007 because they were an outlier in the data
files_EA_resid_ts = sorted(glob.glob("/scratch/a/arisvoin/zara1252/project_1/data_extracted/NEW_DATA_SPINS/SPINS_Removed_Files/output_exclusion_criteria/output/**/sub-*_RS_2mm_GSR_glasser_tian_meants.csv", recursive = True)) #for SPINS

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

```

#The code block above was actually a python rendition of Teo's R code that can be found here: /projects/zkhan/proj1/R_code/SPINS_OLD_CODE/timeseries_variability.R

#What this code block does is formulates the timeseries data that was extracted into the csv's from Step 1 and places them into a python dictionary called "resid_ts_named" where the key is the participant ID and the value is the timeseries data for each participant. Simply, "resid_ts_named" is a massive dictionary variable which stores all timeseries data for each participant (in this case 344) and names them. 

#if you want to change how the files are named - check out "ptlist"


#The third block of code is:
```{r setup, include=FALSE}
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

#### Running neuroComBat to account for scanner effects - couldnt figure it out - i yeeted it into R and then put the csv back into here after calculating the mssd put path to where the neurocombat code is
```

#the code above is kinda the main part of my entire summer project - my baby. Handle with care please. To kinda sum this up, once the cell block is run, a menu appears where you either enter "1" or "2" depending on how you want the timeseries data to be analyzed. From there, the timeseries data is then z-scored or mean-centered and then we get the mssd values after it runs through the mssd_value function. It spits out a dataframe with the participants lining the first column and then their respective mssd values for each of the 392 regions (or however much there is from the parcellated atlas that was used). I won't go into code detail here - because the code is commented - and chat gpt exists.

# The fourth code block is
```{r setup, include=FALSE}
'''GETTING CSVs - make sure to check which matrix variable is produced'''
pd.DataFrame(mssd_matrix_mean).to_csv("old_mssd_matrix_mean.csv") #gets a csv file here
```

#the above code block was to get a csv of the mssd values

#BUT moving along - I didn't get the chance to code neuroCombat - a function used to mitigate any scanner errors - into this code here. So basically, I had to take the csv file of the mssd values from Jupyter Notebook and then put them into the code shown below in R (a snippet taken out of Julia's code :) ). This was done by, again, using the secure copy command through terminal, but now set-up the other way since I was copying from SciNet's Jupyter Notebook to the CAMH Desktops (example shown below): 

# scp - r zara1252@niagara.scinet.utoronto.ca:/scratch/project_1/data_extracted/NEW_DATA_SPINS/mssd_csv_without_neurocombat/mssd_matrix_mean.csv /projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/

#the code for neuroCombat is shown below - the full version can also be accessed through /projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/NEW_plsc.Rmd
#NOTE: all variables that have a file path lead to the correct csv used
  
Load in data for analysis and run neurcoCombat 
```{r}
behavioural <- read.csv('/projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/NEW_SPINS_final_regressiondata.csv')
mssd <- read.csv('/projects/zkhan/proj1/New_Data/SPINS/old_mssd_matrix_mean.csv')
X <- mssd
Y <- behavioural

#non0_behavioural = Y[(!rowSums(Y[,2:9])==0),]
#non0_brain = X[(!rowSums(Y[,2:9])==0),]

X_non01 <- X[,-1] #deletes the first column here
Y_non01 <- Y 

#Creating rownames and colnames based on the first rows/columns for formatting 
X_non0 <- X_non01[-1,]
colnames(X_non0) <-X_non01[1,]
X_non0f <- X_non0[,-1]
rownames(X_non0f) <-X_non0[,1]

Y_non0f <- Y_non01[,-1]
rownames(Y_non0f) <-Y_non01[,1]

#NeuroCombat step 
#Must transpose so that the participants are in the columns and the MSSD values are the rows 
master_2 <- as.matrix(X_non0f) #first convert to a matrix, and then transpose 
master_1 <- t(master_2) #I did this step but idk R and it's weird row/column thing it has going on 

# dat is a data matrix of the data to harmonize - rows are features (connections) and columns are participants
resting_mssd_values <- master_1
numeric_test <- apply(resting_mssd_values, c(1,2), as.numeric) #changes the matrix values into numeric so it doesnt run into issues for the neruocombat step
#as.numeric(resting_mssd_values)

class###Import data to be used as design matrix 
add_data <- read.csv("/projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/SPINS_finalregression_RS_reorderd_excluded.csv")
class(add_data$diagnostic_group)
add_data$diagnostic_group <- as.factor(add_data$diagnostic_group)
class(add_data$diagnostic_group)
class(add_data$scanner)
add_data$scanner <- as.factor(add_data$scanner)

# mod is a design matrix specifying biological covariates that should be protected - here diagnosis, age, sex, and cog variables
modcombat <- model.matrix(~diagnostic_group +sex + age +
                              avg_fd_RS + scog_mean_ea + simulation + mentalizing + wtar_std_score + scog_iri_total +  np_domain_tscore_process_speed +
                            np_domain_tscore_att_vigilance + np_domain_tscore_work_mem +np_domain_tscore_verbal_learning + np_domain_tscore_visual_learning
                          + np_domain_tscore_reasoning_ps + np_domain_tscore_social_cog + bsfs_sec1_total + bsfs_sec2_total + bsfs_sec3_total +
                            bsfs_sec4_total + bsfs_sec5_total + bsfs_sec6_total, data=add_data)

# R run ComBat
# batch is a vector (length should be equal to the number of columns in the data matrix) that specifies the id for the batch, site, or scanner to correct for
rest_combat <- neuroCombat(dat=numeric_test, batch=c(add_data$scanner), mod=modcombat) #try and remove the first column where all the regions are shown

# transpose the harmonized data matrix, give harmonized matrix of MSSD values in the end 
rest_combat_data <- t(rest_combat$dat.combat)

```

Confounds of CPZ equivalence, motion (Mean FD), sex and age were regressed out of both brain and clinical data
```{r warning=FALSE}
#Define confounds.
testing <- read.csv('/projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/SPINS_regressors_excluded.csv') 
column_names <- colnames(testing) #testing to see the column names
con_mat= data.frame(matrix(ncol=4)) #added the rows here
colnames(con_mat) <- c("record_id","sex","age","avg_fd_RS")
non0sub = rownames(X_non0f)


#comment out
for (subject in non0sub){
  con_mat= rbind(con_mat, testing[testing$record_id == subject, c(1,2,3,4)])
}

#removes the first column
con_mat = con_mat %>%
  filter(!row_number() %in% 1) 

#keep this for now nad see what this does - update: helps with formatting
con_mat <- con_mat[,-1]
rownames(con_mat) <- non0sub 

#behavioural matrix
behavmat <- Y_non0f

#Define values.
val_pred <- colnames(behavmat)

#Join confounds and values.
resmat <- cbind(behavmat,con_mat)

#For each value.
for (val_lab in val_pred) {
  
  #Produce residuals.
  inform <- as.formula(paste0(val_lab,'~',con_mat))
  respred <- resid(lm(inform,data=resmat))
  
  #Replace.
  #behavmat[,val_lab] <- respred
  behavmat[, val_lab][!is.nan(behavmat[, val_lab])] <- respred[!is.nan(respred)] #added to skip over any NaN values, COME BACK TO THIS 
}


#brain matrix
brainmat <- rest_combat_data 
colnames(brainmat) <- gsub('-','.',colnames(brainmat))
#Define values.
val_pred <- colnames(brainmat)
#Join confounds and values.
resmat <- cbind(brainmat,con_mat)

#For each value.
for (val_lab in val_pred) {
  
  #Produce residuals.
  inform <- as.formula(paste0(val_lab,'~',con_mat))
  respred <- resid(lm(inform,data=resmat))
  
  #Replace.
  brainmat[,val_lab] <- respred
}
```

#at the very end I had taken the "brainmat" variable and made it into a csv to send back to Jupyter Notebook. 
#the line of code for the "brainmat" variable to become a csv: write.csv(brainmat, file = "neuroCombat_mssd.csv", row.names = TRUE)
#this csv was then "scp" over to Jupyter Notebook


#The fifth code block:
```{r setup, include=FALSE}
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
```

#all this does is write .txt files of the mssd values from a certain csv file for each participant separately that you can update on the bottom of the cell block. No matter what it tells you, it will not dump the .txt files where the path is listed on the bottom. It will dump them in whatever directory the ipynb file is in. Trust me. Don't trust it. 


#The sixth block is:
```{r setup, include=FALSE}
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
```

#basically all this does is order a csv file based on controls vs cases or however the csv is ordered for "df2" and then gets the average value for both cases/controls for each region. You're left with a csv of 2 rows by 392 columns - in this case. 

#The next couple of cell blocks are creating csv's and then creating .txt files for further analysis and plotting.


#The ninth cell block is: 
```{r setup, include=FALSE}
'''CALCULATING THE T-TEST OF AVERAGE CASES VS CONTROLS - GIVES ONE OVERALL'''

import pandas as pd
from scipy.stats import ttest_ind

# Replace 'your_file.csv' with the path to your actual CSV file
data1 = pd.read_csv('neuroCombat_average_mean_centered_mssd_values.csv', header=None)
data_dropped1 = data1.drop(data1.columns[0], axis=1) #drops the first column
data_dropped2 = data_dropped1.drop(data_dropped1.index[0], axis=0) #drops the first row

# Transpose the DataFrame so that controls and cases are in separate columns
data = data_dropped2.T

# Extract controls and cases data into separate Series by columns
controls = data.iloc[:,1]
cases = data.iloc[:, 0]

# Perform the t-test
t_stat, p_value = ttest_ind(cases, controls)

# Output the results
print("T-Statistic:", t_stat)
print("P-Value:", p_value)
```

#this calculates the t_test on the average cases/controls - yay! first statistic done! Shall give you a number and a corresponding p-value. Make sure the inputted csv was the neuroCombatted version! 

#The tenth (we hit double digits O_O) cell block is:
```{r setup, include=FALSE}
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
```

#This basically now does a t-test on each region specfically between the cases and control group. To set the number of columns that a t-test should be done for, change the number for the "num_columns" variable, since the "for" loop runs for that many times. In the end, it gives you a list where it inserts the log_-10_p_values, where the indexed region would be - so it places a numerical value for a significant region, leaving the rest of the regions with a numerical value of 0. 


#The 11th block of code is:
```{r setup, include=FALSE}
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
```

#another .txt file creator. To be fair, I could've just kept one that couldve been used over and over again - but this kinda kept my thoughts in order for when I created csv files and in what order. But this just creates one .txt file where all the log_-10_p_values are stored in their respective region index.

#The 12th block of code is just to check the motion data with a histogram - yolo. Not what's important in thy project tho. Hence why she is commented out.  


#OKAY SO - IT'S PLSC TIME

# 4) So basically my PLSC code is found here: /projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/NEW_plsc.Rmd

#I had coded and used versions: R/4.1.0 and rstudio/1.4.1717

# This was taken from Julia's code and adapted wtih csv's needed for PLSC. Nothing was drastically changed and the lines of code have been commented 

#For this project, we had found Latent Variable 1 to be stable and significant to be included within the analysis. All code should work up until you reach "Important Brain Connection contributors for LV3 /n original connections valence weighted by contributions" . I only focused on LV1, so if you scroll down to "Important Brain Connection contributors for LV1 /n original connections valence weighted by contributions" - some lines of code work here as I transitioned from R to Jupyter Notebook again. 

# 5) THE MAJOR SWITCHAROO
# To find the top ten positive and negative loadings, I needed to send a csv of the variable "Cx_conn" to Jupyter Notebook. This csv can be found here:  /projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/NEW_SPINS_output_file.csv

#After sending this file onto Jupyter Notebook, I used the code called "fixing_file_for_brain_plot_in_r_backup.ipynb" found within this folder: /scratch/project_1/random_backup/NEW_SPINS_output_file.csv.

#Make sure that the Cx_conn csv is within the same directory as the python code. 

#The code block below is basically what was used to amplify the values, you need a csv with ROIS, area, network all filled in and then an "r" column filled with 0's so that it can be replaced by the first column of the Cx_conn csv in the "file" variable below.
```{r setup, include=FALSE}
file = pd.read_csv('NEW_SPINS_output_file.csv')
plsc_brain = pd.read_csv('NEW_rs_conn_LV1_correct_with_zeros.csv')

#extracting the first column because that was used for LV1
column = file.iloc[:, 1]

#multiplying values by 10000
big_column = column*10000

#changes the values in a column named 'r' all to 0
plsc_brain['r'] = 0 
#attach the amplified column values to the end of the other csv
plsc_brain['r'] += big_column

#delete the first two columns:
plsc_brain_drop = plsc_brain.drop(plsc_brain.columns[0], axis=1)
#plsc_brain_drop2 = plsc_brain_drop.drop(plsc_brain_drop.columns[0], axis=1) #weird cause it doesnt consider the first unnamed column in the csv to be anything
```

#this was then yeeted into a csv through: pd.DataFrame(plsc_brain_drop).to_csv("NEW_rs_conn_LV1.csv"). This csv was then brought back into R using the "scp" command in terminal.

#From there, a couple lines of code were used here:
```{r setup, include=FALSE}
rs_conn_LV1 <- read.csv('/projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/NEW_rs_conn_LV1.csv') #changed from jupyter notebook
lv1_pos_loadings = rs_conn_LV1[order(-rs_conn_LV1$r),] #put the csv back in after i add in the other column with it matched up the regions
top_10_lv1 = lv1_pos_loadings[1:10,]
```

#this gave the top 10 positive loadings, and to get the top 10 negative loadings it would be:
```{r setup, include=FALSE}
lv1_neg_loadings = rs_conn_LV1[order(rs_conn_LV1$r),]
top_10_neg_lv1 = lv1_neg_loadings[1:10,]
```

#From here, the top 10 region values were taken and manually placed into a previously filled "0" value csv respectfully.


#6) NOW HERE COMES THE BRAIN PLOTS
# so basically I adapted Julia's code for this in R. 
# My code can be found at /projects/zkhan/proj1/R_code/SPINS_NEW_R_CODE/testing_brain_plots. 
# Minor things to note: you need to install "ggsegGlasser" by enabling the universe instead of the common command of install.packages("ggsegGlasser"), and you have to manually input the loading values into the functions. 
#Another thing is that the region names are the letters/numbers in the middle of a region name

# But that's about it! Thanks for reading :p

#The Jupyter Notebook code can be found at: /projects/zkhan/proj1/JupyterNotebook_Code/project_1/data_extracted/NEW_DATA_SPINS/Project_1_SPINS.ipynb

