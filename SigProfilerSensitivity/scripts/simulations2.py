import pandas as pd 
import numpy as np
from numpy import random

def inject_signature(input_matrix, percentage, sd, signature_weights=[], samples=[]):

    xs = [] #of injected mutations
    sds = []

    
    df = pd.read_csv(input_matrix, sep="\t")

    #either operate on a specific set of samples or on all samples
    if len(samples) == 0:
        samples = list(df.columns)
    else:
        samples = [df.columns[0]] + samples
    df = df[samples]
    
    #print(len(df.columns))

    assert(percentage >= 0 and percentage <= 100)
    assert(len(signature_weights == df.shape[0]))
    

    matrix = [list(df[df.columns[0]])]
    for i, col in enumerate(df.columns[1:]):
        counts = np.array(list(df[col]))
        if len(samples) != 0 and col not in samples: #ignore if sample not of interest
            matrix.append(counts)
        else:
            print("Injecting signature...")

            #determine mean and variance of signature to be injected as Gaussian noise
            total = np.sum(counts)
            mean = (percentage / 100) * total #mean TMB 
            stddev = (sd / 100) * mean #standard deviation as a percentage of the mean
            #print(total, mean, stddev)

            #determine # of mutations to inject based on Gaussian  
            x = 0
            while x < 1: #make sure we're injecting at least 1 mutation 
                x = round(random.normal(loc=mean, scale=stddev))
                
            #VECTOR ADDITION BASED ON SIGNATURE
            injected_mutations = signature_weights * x
            injected_mutations = [round(n) for n in injected_mutations] #round to nearest whole number

            xs.append(x) #keep track of how many mutations we injected
            sds.append(stddev)

            #randomly remove x mutations before adding any mutations to maintain TMB of sample
            counter=0
            for i in range(sum(injected_mutations)):
                random_index = random.randint(0,len(counts)-1)
                print(counts[random_index])
                #print("Entering 2nd while loop")
                while counts[random_index] - 1 < 0: #pick a random index to remove a mutation, while ensuring that we don't have negative counts
                    random_index = random.randint(0,len(counts)-1) #pick another random index
                counts[random_index] = counts[random_index] - 1 

            #after we have subtracted mutations, we can now inject mutations as Gaussian noise      
            result = counts + injected_mutations #new count vector
            matrix.append(result)
    #print (matrix[1:10])

    result_df = pd.DataFrame(matrix).transpose()
    result_df.columns = df.columns

    #Ensure that TMB was maintained
    assert(result_df.iloc[:, 1:].to_numpy().sum() == df.iloc[:, 1:].to_numpy().sum())
    #assert(len(result_df.columns) == 248)
    #print(result_df.iloc[:, 1:].to_numpy().sum(), df.iloc[:, 1:].to_numpy().sum())
    return result_df


# passive_smoking = []
# samples = pd.read_csv("/data/khandekara2/dev/Sherlock-Lung/sherlock_1217_information.csv")
# for s, status in zip(samples["Tumor_Barcode"], samples["Passive_Smoking"]):
# 	if status == "Yes":
# 		passive_smoking.append(s)
# print(len(passive_smoking))


# percentages = []
# with open("/data/khandekara2/dev/Sherlock-Lung/scripts/SBS288N_simulations.swarm", "w") as f:
# 	input_matrix = "/data/khandekara2/dev/Sherlock-Lung/data/CRC_Manuscript_v1.SBS288.all"

# 	signature = pd.read_csv("/data/khandekara2/dev/Sherlock-Lung/data/SBS288N.txt", sep="\t")
# 	signature_weights = np.array(list(signature["SBS288N"])) #probability distribution of the signature over the channels

# 	percentages = [1, 5, 10, 15, 20]
# 	trials = [i for i in range(1,11)]
# 	sd=10

	for percentage in percentages:
		p = str(percentage)
		for trial in trials:
			t = str(trial)
			simMatrix = inject_signature(input_matrix, percentage, sd, signature_weights=signature_weights, samples=[])
			simMatrix.to_csv("/data/khandekara2/dev/Sherlock-Lung/data/CRC-simulations/SBS288N/SBS288N_injected_" + str(p) + "_trial" + str(t) + ".matrix.tsv", sep="\t", index=None)
			#_injected_1_trial9_
			project = "SBS288N_injected_" + str(p) + "_trial" + str(t)
			f.write("python3 /data/khandekara2/dev/Sherlock-Lung/data/simulations3/run.py /data/khandekara2/dev/Sherlock-Lung/data/simulations3/SBS288N_injected_" + str(p) + "_trial" + str(t) + ".matrix.tsv " + project+"\n") 

