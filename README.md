# RECALL Benchmarks

This benchmarks based on scaffold recall analyses. This repository consist of four different benchmarks to compare, to evaluate molecular generator.


### ▪ True positive recall all (TUPOR): 
the number of unique active scaffolds generated for the entire OS (UASo) in the ration to the total number of unique active scaffolds in RS (UASr). Quantifies the recall number of the unique active scaffolds for the entire structure generation process. TUPOR value ranges between 0 (no unique active scaffolds from RS were found in the OS) and 1 (all unique scaffolds from RS were found in the OS).  
   $$TUPOR = {{UASo} \over UASr}$$ 
   
* **Use case** : A high TUPOR value suggests that the model successfully recalls biologically active scaffolds, making it useful for generating molecules with pharmaceutical potential. Users can evaluate the model's accuracy in generating relevant scaffolds for drug discovery. 
### ▪ Set scaffold yield (SESY): 
value for the whole generated set, expressing how many different unique scaffolds (USo) are within their OS size (SSo). This metric reflects the diversity of the generated set by calculating how many different scaffolds are present in the OS. A value of 1 indicates that every scaffold in the OS is unique, while 0 means no valid scaffolds were generated. 

   $$SESY = {{USo} \over SSo}$$ 

* **Use case**: SESY is valuable for users aiming to strike a balance between exploration and exploitation. A high SESY value indicates that the generator explores diverse regions of chemical space by discovering novel scaffolds, which is essential for expanding the chemical space. On the other hand, a low SESY value suggests a focus on exploiting known scaffolds, refining existing chemical structures. By adjusting the SESY, users can fine-tune the generator's ability to either discover new molecular structures (exploration) or optimize known scaffolds (exploitation), depending on the needs of their project. 

### ▪ Absolute set scaffold recall (ASER): 
Absolute set scaffold recall (ASER) – value for the generated dataset expressing the total number of virtual compounds with active scaffold in the OS (CwASo) in relation to other scaffolds in OS (SSo - CwASo). It tells us how effectively the generator is exploring the desired region of chemical space, with higher values reflecting better alignment with biologically relevant scaffolds. 

   $$ASER = {{CwASo} \over SSo}$$ 

* **Use case**: A high ASER value indicates that the generator is targeting the most relevant chemical space for drug discovery, creating compounds that are more likely to exhibit biological activity. Users can rely on this metric to assess how well the generator is focusing on biologically active areas of molecular space. 

### ▪ Active scaffold ratio (ACR): 
Active scaffold ratio (ACR) – the number of unique active scaffolds generated for the entire OS (UASo) in the ration to other unique scaffolds in OS (USo - UASo). This metric calculates the ratio of active scaffolds to total scaffolds in the OS. A higher ASR value suggests that the model is generating a larger proportion of biologically active scaffolds compared to other types of scaffolds. 

   $$ACR = {{UASo} \over {{USo} - {UASo}}}$$ 

* **Use case**: ASR is useful for users prioritizing the biological relevance of generated compounds. A high ASR value means that the generator is more likely to produce active scaffolds, which is essential for applications that focus on optimizing lead compounds with desired biological properties. 


## Installation


# Example














$\frac{1}{2}$



# My calculation procedure:
1. Take the cIS and run the generator, 10 times because 5 times for Dissimilarity split and 5 times fo Similarity split
2. compound Output Set convert to one column set (used ***generated_set_convert_to_one_column.ipynb*** )
3. Used One column set for calculation metrics (used ***calculation_the_metrics.ipynb***)
4. Calculate the avarage value between this 5 calculations (used ***calculation_the_metrics.ipynb***)
5. Visualization (used ***visualization.ipynb***)

# Data:
Because all data(cOS) to large, I add to Github compress files.


# Content:
* All Folders are Bold
* All csv in Italic
* All .ipynb and .py in BOld and Italic
* dis = Dissimilarity split
* sim = Similarity split
* cIS = compound Input Set
* cRS = compound Recall Set
* cOS = compound Output Set

1. **data** -> consists of all data compound Input Sets, compound Output Sets, compound Recall Sets, Raw data, Calcelated Matrics:
   1. **chembl_data** -> the raw data from ChEMBL DB 31 for nuclear receptor and leukocyte elastase
   2. **input_recall_sets** -> compound Input Sets and compound Recall Sets for generators, divided depending on the target
      1. **Glucocorticoid_receptor** -> In ths folder you can see:
         1.  *Glucocorticoid_receptor.csv* -> the table with all active compounds and other information from DB
         2.  *Glucocorticoid_receptor_split_to_clusters_using_KMedoids.csv* -> after splitting to 5 clusters the table with active compounds with the column about cluster
         3.  *cIS_Glucocorticoid_receptor_dis_0.csv* -> the example for name cIS for Glucocorticoid receptro, dis -> the type of splitting (Dissimilarity split = dis, Similarity split = sim), 0-> the cluster (can change from 0 to 4), because we devided to 5 clusters, This csv is compound Input Set for all generators, consist of only one column with SMILES of active compounds
         4.  *cIS_Glucocorticoid_receptor_with_ID_dis_0.csv* -> cIS with ChEMBL_ID. the similar example like a previous but this csv consist of two column: SMILES for active compounds and ID for active compounds
         5.  *cIS_Glucocorticoid_receptor_with_p_chembl_dis_0.csv* -> cIS with p_chembl_value. the similar example like a previous but this csv consist of two column: SMILES for active compounds and p_chemble value for active compounds
         6.  *cIS_Molpher_Glucocorticoid_receptor_dis_0.csv* -> cIS for Molpher generator, because Molpher need another type of input data. dis = Dissimilarity split, 0 -> 0.cluster
         7.  *cRS_Glucocorticoid_receptor_dis_1.csv* -> compound Recall Set for Glucocorticoid receptor, dis = Dissimilarity split, 1 = 1. cluster
         8.  *cRS_Glucocorticoid_receptor_with_ID_dis_0.csv* -> cRS with ChEMBL_ID. the similar example like a previous but this csv consist of two column: SMILES for active compounds and ID for active compounds
         5.  *cRS_Glucocorticoid_receptor_with_p_chembl_dis_0.csv* -> cRS with p_chembl_value. the similar example like a previous but this csv consist of two column: SMILES for active compounds and p_chemble value for active compounds
      2. **Leukocyte_elastase** -> the same like in Glucocorticoid receptor but for Leukocyte elastase.The same type of csv files.
   3. **nuclear_receptor** -> raw data for Glucocorticoid receptor
      1. *25_pic50.csv*
   4. **output_sets** -> the folder consists of all compound Output Sets,divided depending on the target:
      1. **Glucocorticoid_receptro**: -> cOS for Glucocorticoid receptor: (the parameters what can change in the name of result, the name of generaor, type of split, number of cluster, all_columns/one_column)
         1. *cOS_Drug_Ex_dis_1_all_columns.csv* -> cOS for DrugEx, dis = Dissimilarity split, 1 = 1.cluster, all_columns -> the raw data after running the generator the it converted to one_columns data, because for calculation metrics we need only virtual compounds -> for this purpose existing the notebook ***generated_set_converte_to_one_column.ipunb***
         2. *cOS_Drug_Ex_dis_0_one_column.csv* -> cOS for DrugEx what we used for calculation metrics, the table consist of one column it's SMILES by virtual compounds
      2. **Leukocyte_elastase** -> the same.
   5. **protease** -> raw data for Leukocyte elastase
      1. *235_pic50.csv*
   6. **result_data_for_chosen_perfect_target** -> some csv flies which I save during selecting
   7. **results** -> the folder where the results are stored, after run the notebook (***calculate_the_metrics.ipynb***), divided depending on the target and type of scaffold and type pf slit and the name of generator. Resulting folder must create automaticly.
      1. **Glucocorticoid_receptor**: -> type of target
         1. **csk_scaffolds**: -> type of scaffolds
            1. **dis**: -> Dissimilarity split.
               1. **Molpher** -> the folder with results, consist of csv with calculated metrics and scaffolds in cOS and cRS and the table of occurance of Recall scaffolds:
                  1. *metrics_cluster_0_dis_Molpher.csv* -> only one row with results.
                  2. *count_of_occurance_cluster_0_dis_Molpher.csv* -> the helping table what I saved during metric calculation
                  3. *scaffolds_of_output_set_cluster_0_molpher.csv* -> the scaffolds what I saved during calculation metrics, maybe I can use them then
                  4. *scaffolds_of_recall_set_cluster_0_dis_Molpher.csv* -> the same but for compound Recall Set
               2. **DrugEx** -> same
               3. **Crem** -> same
               4. *Molpher_mean_csk_dis.csv* -> table with mean value of Dissimilarity split by Molpher
            2. **sim**: -> Similarity split.
         2. **murcko_scaffolds**:
            1. **dis**: same
            2. **sim**: same
2. **scripts** -> folder consist of helping scripts for calculation metrics.
   1. ***metrics.py*** -> script needed for calculation 3 importnant metrics (TUPOR, ASER, SESY). This script are used in ***calculation_the_metrics.ipynb***
   2. ***metrics_all.py*** -> old script with calculation all 7 metrics
   3. ***preprocesing.py*** -> helping script during creating sets (cIS, cRS). This script are used in ***generate_data_active.ipynb***
   4. ***visualization.py*** -> helping script durring visualization, but stil not the final version! This script are used in ***visualization.ipynb***
3. **img** -> folder for images
4. ***calculation_the_metrics.ipynb*** -> the notebook what I run for calculation metrics
5. ***checking_scaffolds.ipynb*** -> only check if some scaffolds not in the same time in cIS and cRS.
6. ***generate_data_active.ipynb*** -> the notebook for generating cIS and cRS for dis and sim
7. ***generate_data_non_active.ipynb*** -> the notebook ehat I need for generate dataset for non-active compounds. Non-active compounds we need for create QSPRpred.
8. ***select_the_biological_target.ipynb***-> the notebook when I do selecting the targets (Glucocortiocoid receptor and Leukocyte elastase) -> need some cleaning :)
9. ***visualization.ipynb*** -> the notebook when I do dome vizualization, but still not final version, because we need better vizualization than I have 
10. ***results_of_project.md*** -> here I try to write the results and all important informations
11. ***pca.py*** -> script what I need to visualization logistic PCA (https://github.com/brudfors/logistic-PCA-Tipping)
12.  ***generated_set_convert_to_one_column.ipynb*** -> converted COS to one column set what we need for inout file for metric calculation



