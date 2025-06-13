'''Import necessary libraries'''
import pandas as pd
from rdkit import Chem
import os
import numpy as np
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.Lipinski import NumAromaticRings
from rdkit.Chem.Lipinski import NumAliphaticRings
from sklearn.metrics import jaccard_score
from scipy.spatial.distance import euclidean, pdist, squareform
from sklearn_extra.cluster import KMedoids
import itertools
from itertools import combinations
from rdkit.Chem import inchi



class Preprocesing():
    def __init__(self, receptor_name, save_options = True):
        self.active_compounds = pd.DataFrame()
        self.active_compounds_with_clusters = pd.DataFrame()
        self.receptor_name = receptor_name
        self.save_options = save_options
    

    def MakeScaffoldGeneric_fixed(self, mol):
        '''Function wich need to fixed scaffold generic, because if existing SF6'''
        rxn = AllChem.ReactionFromSmarts\
                ('[*:0][S:1]([F:2])([F:3])([F:4])([F:5])[F:6]>>[*:0]-[C](-[F])(-[F])-[F].[*:1]([*:2])([*:3])([*:4])([*:5])[*:6]')
        for i in range(len(mol.GetSubstructMatches(Chem.MolFromSmiles("S(F)(F)(F)(F)F")))):
            products = rxn.RunReactants((mol,)) # tuple
            if len(products)>0:
                mol = products[0][0]
        return MakeScaffoldGeneric(mol)
    
    
    def duplicate_compounds(self, df):
        '''Function for find and work with duplicates'''

        duplicateRows = df[df['canonical_smiles'].duplicated()]
        for x in duplicateRows.index:
            print("In yours dataset we found some duplicates, we calculate mean value and then we delete another rows and leave only one")
            smiles = duplicateRows.canonical_smiles[x]  #smiles of duplicates
            duplicated_df = df[df.canonical_smiles==smiles] #find other duplicates row
            duplicated_df_index = duplicated_df.index #all index in duplicates dataframe
            new_value = duplicated_df.pchembl_value.mean() #new value, select all value in duplicates set and then calculate mean value
            inchikey = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(smiles)) #calculate the inchikey
            leave_id = df[df.stand_inchi==inchikey].index #index wich we leave in dataframe and others duplicqates we deleted
            df.loc[leave_id, 'pchembl_value'] = new_value #add new value
            duplicated_df_index = duplicated_df_index.drop(leave_id) # inbdex wich we need to delete
            for x in duplicated_df_index:
                df = df.drop([x])  #delete duplicates

        '''Only check if the numbers is the same'''
        array_inchi_key = [Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(x)) for x in df.canonical_smiles]
        if len(df.canonical_smiles.unique()) == len(set(array_inchi_key)) == len(df):
            print("Everything is okay")
        else:
            print("Something wrong with your duplicates")
        df.reset_index(level=0, inplace=True)
        df = df.drop(columns = ['index'])
        return df


    def load(self,type_receptor, name):
        '''Load all data needed, add header and then convert canonical smiles to scaffolds csk and calculate Morgan fingerprint with r=3, then save new dataset and save to contructor self.active_compounds'''
        dff = pd.read_csv(f"data/{type_receptor}/{name}",header=None)
        dff.columns =['molregno', 'stand_type', 
                                        'pchembl_value', 'stand_value',
                                        'canonical_smiles', 'stand_inchi',
                                       'chembl_id', 'tid','pref_name']


        dff['scaffolds_csk'] = [1 for x in range(len(dff))]
        delete_element = 0
        for x in range(len(dff)):
            
            remover = SaltRemover()
            res = remover.StripMol(Chem.MolFromSmiles(dff['canonical_smiles'][x]))
            #dff['canonical_smiles'][x] = Chem.MolToSmiles(res)
            dff.loc[x, 'canonical_smiles'] = Chem.MolToSmiles(res)
            if NumAromaticRings(Chem.MolFromSmiles(dff['canonical_smiles'][x])) == 0 and\
                NumAliphaticRings(Chem.MolFromSmiles(dff['canonical_smiles'][x])) == 0:
                print("Not contain ring")
                print(dff['canonical_smiles'][x])
                print(Chem.MolFromSmiles(dff['canonical_smiles'][x]))
                dff = dff.drop([x])
                delete_element = 1
            if delete_element != 1:
                try:
                    dff.loc[x, 'scaffolds_csk'] = MurckoScaffoldSmiles(\
                                             Chem.MolToSmiles(MakeScaffoldGeneric\
                                           (Chem.MolFromSmiles(dff['canonical_smiles'][x]))))
                except:
                    print("Faild to create scaffold_csk")
                    print("Index",x)
                    print(dff['canonical_smiles'][x])
                    print(Chem.MolFromSmiles(dff['canonical_smiles'][x]))
                    try:
                        mol = self.MakeScaffoldGeneric_fixed(Chem.MolFromSmiles(dff['canonical_smiles'][x]))
                        print(mol)
                        dff.loc[x, 'scaffolds_csk'] = MurckoScaffoldSmiles(Chem.MolToSmiles(mol))
                        print(Chem.MolFromSmiles(MurckoScaffoldSmiles(Chem.MolToSmiles(mol))))
                    except:
                        dff = dff.drop([x])
            delete_element = 0
        dff['mfp'] = [(AllChem.GetMorganFingerprintAsBitVect\
                    (Chem.MolFromSmiles(i),3, nBits=2048)) for i in dff['scaffolds_csk']]

        dff.reset_index(level=0, inplace=True)
        dff = dff.drop(columns = ['index'])
        #for glucocorticoid receptor we need to delete duplicates then, during create test/train sets
        if self.receptor_name != 'Glucocorticoid_receptor':
            dff = self.duplicate_compounds(dff)
        
        self.active_compounds = dff.copy()

        if self.save_options == True:
            if not os.path.exists(f"data/input_recall_sets/{self.receptor_name}/"):
                os.makedirs(f"data/input_recall_sets/{self.receptor_name}/")
            dff.to_csv(f'data/input_recall_sets/{self.receptor_name}/{self.receptor_name}.csv', index_label = False)
        return dff
    

    def split_data_for_clusters_KMedoids(self, dff, random_state_number=366):
        '''Function for spliting all data to 5 clusters. Using K-medoid algorithm with jaccard distance. We selected number of random state, because I need the similar number of compounds in each clusters'''

        '''Select only unique active scaffolds'''
        df = dff.drop_duplicates(subset='scaffolds_csk', keep="first")
        df.reset_index(level=0, inplace=True)
        df = df.drop(columns = ['index'])

        '''Selected morgan fingerprint, wich is was calculated in previous funbction'''
        fps = [x for x in df['mfp']]
        random_state = []
        #for glucocorticoid receptor random_state = 1160
        #for leukocyte elastase random_state = 366
        try_random_state = [random_state_number]
        for x in try_random_state: 
            #print("AAAAA")
            kmedoids = KMedoids(n_clusters=5,metric="jaccard", random_state = x, \
                                init='k-medoids++').fit(fps)
            labels = kmedoids.labels_

            print("Random_state: ", x)
            print("Number in 0 cluster: ", len(labels[labels==0]),"Number in 1 cluster: ",len(labels[labels==1]),\
                  "Number in 2 cluster: ",len(labels[labels==2]),"Number in 3 cluster: ",len(labels[labels==3]),"Number in 4 cluster: ",len(labels[labels==4]))
            random_state.append(x)
    
            dff = pd.DataFrame(index=np.arange(5), columns=np.arange(5))
            for x in range(5):
                for y in range(5):
                    dff[x][y] = jaccard_score(kmedoids.cluster_centers_[x],kmedoids.cluster_centers_[y])

            print("Matrix of similarity of each clusters:")
            print(dff)        

        df['clusters'] = labels

        self.active_compounds_with_clusters = df.copy()

        #print(df)
        non_active_scaffold = pd.read_csv(f"data/input_recall_sets/{self.receptor_name}/df_not_in_new_active_sets_new.csv")
        deleted_index = []
        for x in range(len(df)):

            if df.loc[x,'scaffolds_csk'] in list(non_active_scaffold.scaffolds_csk):
                deleted_index.append(x)

        for r in deleted_index:
            df = df.drop([r])
        df = df.reset_index(drop=True)
        #print(df)
        if self.save_options == True:
            df.to_csv(f'data/input_recall_sets/{self.receptor_name}/{self.receptor_name}_split_to_clusters_using_KMedoids.csv', index_label = False)

        return df
    

    #Next function about save data to sets for Molpher and ohter generators (DragEx)
    '''Data we splitting for 5 clusters and ih the '''
    def chembl_id_for_train_test_sets_dis(self,data_target,data_cluters, train, test):
        '''Selected chembl_id for following saving to sets for Molpher and other generators. molpher need only one copy and the most active, and other generated need all molecules, wich have the same scaffold'''
        '''This function for dissimilarity calculation, where 4 clusters it's out train set and 1 cluster it's the test set'''
        input_id_Molpher_train_dissimilar = []
        input_id_train_dissimilar = []
        input_id_test_dissimilar = []
        print('Traine index of cluster', train)
        print('Test index of cluster', test)
        count = 0
        for r in train:
            nazev = "clusters"
            a = data_cluters[data_cluters[nazev]==r]
            for x in a.index:
                smiles = a['scaffolds_csk'].loc[x]
                b = data_target[data_target['scaffolds_csk']==smiles]
                b = b.sort_values(by=['pchembl_value'], ascending=False)    
                input_id_Molpher_train_dissimilar.append(b['chembl_id'].iloc[:1].item())
                for y in b.index:
                    input_id_train_dissimilar.append(b['chembl_id'].loc[y])

        #test_data
        nazev = "clusters"
        a = data_cluters[data_cluters[nazev]==test]
        for x in a.index:
            smiles = a['scaffolds_csk'].loc[x]
            b = data_target[data_target['scaffolds_csk']==smiles]
            b = b.sort_values(by=['pchembl_value'], ascending=False)    
            for y in b.index:
                input_id_test_dissimilar.append(b['chembl_id'].loc[y])


        print("Set size for train set for Molpher:",len(input_id_Molpher_train_dissimilar))
        print("Set size fot train set for other generators:",len(input_id_train_dissimilar))
        print("Set size for test set:",len(input_id_test_dissimilar))
        return input_id_Molpher_train_dissimilar,input_id_train_dissimilar, input_id_test_dissimilar
    

    def chembl_id_for_train_test_sets_sim(self,data_target,data_cluters, train, perc_test):
        '''Selected chembl_id for following saving to sets for Molpher and other generators. molpher need only one copy and the most active, and other generated need all molecules, wich have the same scaffold'''
        '''This function for similarity calculation, where 80% in each clusters is train set and 20% is test set'''
        input_id_Molpher_train_similar = []
        input_id_train_similar = []
        input_id_test_similar = []
        count = 0
        for r in train:
            print("R",r)
            train_set = pd.DataFrame()
            nazev = "clusters"
            a = data_cluters[data_cluters[nazev]==r]
            df_split = np.array_split(a, 5)
            test_set = df_split[perc_test]
            print("Len_test",len(test_set))

            for x in range(5):
                if x!=perc_test:
                    train_set = train_set.append(df_split[x])
            print("Train_set",len(train_set))

            for x in train_set.index:
                smiles = a['scaffolds_csk'].loc[x]
                b = data_target[data_target['scaffolds_csk']==smiles]
                b = b.sort_values(by=['pchembl_value'], ascending=False)    
                input_id_Molpher_train_similar.append(b['chembl_id'].iloc[:1].item())
                for y in b.index:
                    input_id_train_similar.append(b['chembl_id'].loc[y])

            #test_data
            for x in test_set.index:
                smiles = a['scaffolds_csk'].loc[x]
                b = data_target[data_target['scaffolds_csk']==smiles]
                b = b.sort_values(by=['pchembl_value'], ascending=False)  

                for y in b.index:
                    input_id_test_similar.append(b['chembl_id'].loc[y])


        print("Set size for train set for Molpher:",len(input_id_Molpher_train_similar))
        print("Set size fot train set for other generators:",len(input_id_train_similar))
        print("Set size for test set:",len(input_id_test_similar))
        return input_id_Molpher_train_similar,input_id_train_similar, input_id_test_similar
    

    def sets_for_Molpher(self,df, data):
        '''Function for generated input set for Molpher, create pair between each others co if we have x compound then we get the x^2-x pair what we run to Molpher'''
        dff = pd.DataFrame(columns = ['start_id','stop_id','start_smiles','stop_smiles'])
        a = list(itertools.permutations(df,2))
        for x in a:
            start_id = x[0]
            stop_id = x[1]
            start_smiles = data[data['chembl_id']==start_id]\
                            ['canonical_smiles'].item()
            stop_smiles = data[data['chembl_id']==stop_id]\
                            ['canonical_smiles'].item()
            dff.loc[len(dff)] = [start_id,stop_id,start_smiles,stop_smiles]
        dff['id'] = [x for x in range(len(dff))]
        new_columns = ['id','start_id','stop_id','start_smiles','stop_smiles']
        dff = dff[new_columns]

        return dff
    

    def sets_for_other_gen_tran_and_test_with_pchembl(self,df, data):
        '''Function for generated sets for other generated, the output consist of only compound with pchemble value'''
        dff = pd.DataFrame(columns = ['smiles','p_chembl'])
        a = []
        b = []
        for x in df:
            chembl_id = x        
            smiles = data[data['chembl_id']==chembl_id]\
                            ['canonical_smiles'].item()
            p_chembl = data[data['chembl_id']==chembl_id]\
                            ['pchembl_value'].item()
            dff.loc[len(dff)] = [smiles,p_chembl]
            #this if check if datasets condist of duplicates
            if smiles in a:
                b.append(smiles)
            else:
                a.append(smiles)

        for x in b:
            indexs = dff[dff['smiles']==x].index
            index_1 = dff[dff['smiles']==x].index[0]
            r = dff[dff['smiles']==x]['p_chembl']
            p_chembl = r.mean()
            for y in indexs:
                if y != index_1:
                    dff = dff.drop(y, axis = 0)             

            dff.loc[index_1,'p_chembl'] = p_chembl
        dff = dff.reset_index()
        dff['id'] = [x for x in range(len(dff))]
        new_columns = ['id','smiles','p_chembl']
        dff = dff[new_columns]

        return dff
    

    def sets_for_other_gen_tran_and_test_with_ID(self,df, data):
        '''Function for generated sets for other generated, the output consist of only compound with ID'''
        dff = pd.DataFrame(columns = ['chembl_id','smiles'])
        a = []
        b = []
        for x in df:
            chembl_id = x        
            smiles = data[data['chembl_id']==chembl_id]\
                            ['canonical_smiles'].item()

            dff.loc[len(dff)] = [chembl_id,smiles]
            #this if check if datasets condist of duplicates
            if smiles in a:
                b.append(smiles)
            else:
                a.append(smiles)

        for x in b:
            indexs = dff[dff['smiles']==x].index
            index_1 = dff[dff['smiles']==x].index[0]
            for y in indexs:
                if y != index_1:
                    dff = dff.drop(y, axis = 0)             
        dff = dff.reset_index()
        dff['id'] = [x for x in range(len(dff))]
        new_columns = ['id','chembl_id','smiles']
        dff = dff[new_columns]

        return dff
    

    def sets_for_other_gen_tran_and_test(self,df, data):
        '''Function for generated sets for other generated, the output consist of only compounds'''
    
        a = []
        b = []
        for x in df:
            chembl_id = x        
            smiles = data[data['chembl_id']==chembl_id]\
                            ['canonical_smiles'].item()
            #this if check if datasets condist of duplicates
            if smiles in a:
                b.append(smiles)
            else:
                a.append(smiles)

        dff = pd.DataFrame(data = a)

        return dff
    

    def split_data_to_train_test_dis(self,data_target,data_clusters):
        '''The main function for create sets for dissimilarity type of project, using the previous function'''
        #data_target = self.duplicate_compounds(data_target)
        data_target = pd.read_csv(f"data/input_recall_sets/{self.receptor_name}/{self.receptor_name}_active_compounds.csv")

        for x in ([[1,2,3,4],0],[[0,2,3,4],1],[[0,1,3,4],2],[[0,1,2,4],3],[[0,1,2,3],4]):
            name = []
            '''Get the Chemble_ID for Molpher, other generators and tests set'''
            id_train_Molpher,id_train_all_gen,id_test_set = self.chembl_id_for_train_test_sets_dis\
                        (data_target,data_clusters,x[0],x[1])
            
            '''Only add SMILES based on Chemble_ID'''
            input_Molpher_train = self.sets_for_Molpher(id_train_Molpher, data_target)
            input_gener_train_with_ID = self.sets_for_other_gen_tran_and_test_with_ID\
                            (id_train_all_gen,data_target)
            input_test_with_ID = self.sets_for_other_gen_tran_and_test_with_ID\
                            (id_test_set,data_target)

            input_gener_train = self.sets_for_other_gen_tran_and_test\
                            (id_train_all_gen,data_target)
            input_test = self.sets_for_other_gen_tran_and_test\
                           (id_test_set,data_target)

            input_gener_train_with_p_chembl = self.sets_for_other_gen_tran_and_test_with_pchembl\
                            (id_train_all_gen,data_target)
            input_test_with_p_chembl = self.sets_for_other_gen_tran_and_test_with_pchembl\
                            (id_test_set,data_target)

            '''Option to save'''
            if self.save_options == True:
                print("Save options")
                #Input Set Molpher
                input_Molpher_train.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_Molpher_{self.receptor_name}_dis_{x[1]}.csv", index=False)

                #Input Set other generators 
                input_gener_train.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_{self.receptor_name}_dis_{x[1]}.csv", index=False,header=False)
                input_gener_train_with_ID.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_{self.receptor_name}_with_ID_dis_{x[1]}.csv", index=False)
                input_gener_train_with_p_chembl.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_{self.receptor_name}_with_p_chembl_dis_{x[1]}.csv", index=False,header=False)
                
                #Recall Set
                input_test_with_p_chembl.to_csv(f"data/input_recall_sets/{self.receptor_name}/cRS_{self.receptor_name}_with_p_chembl_dis_{x[1]}.csv",  index=False,header=False)
                input_test_with_ID.to_csv(f"data/input_recall_sets/{self.receptor_name}/cRS_{self.receptor_name}_with_ID_dis_{x[1]}.csv",  index=False)
                input_test.to_csv(f"data/input_recall_sets/{self.receptor_name}/cRS_{self.receptor_name}_dis_{x[1]}.csv",  index=False,header=False)

    
    def split_data_to_train_test_sim(self, data_target,data_clusters):
        '''The main function for create sets for similarity type of project, using the previous function'''
        #data_target = self.duplicate_compounds(data_target)
        data_target = pd.read_csv(f"data/input_recall_sets/{self.receptor_name}/{self.receptor_name}_active_compounds.csv")
        #similarity:
        for x in ([[0,1,2,3,4],0],[[0,1,2,3,4],1],[[0,1,2,3,4],2],[[0,1,2,3,4],3],[[0,1,2,3,4],4]):
            '''Get the Chemble_ID for Molpher, other generators and tests set'''
            input_id_Molpher_train_sim,\
            input_id_gen_train_sim,\
            input_id_test_sim = self.chembl_id_for_train_test_sets_sim\
                        (data_target,data_clusters,x[0],x[1])

            '''Only add SMILES based on Chemble_ID'''
            input_Molpher_train = self.sets_for_Molpher(input_id_Molpher_train_sim, data_target)
            input_gener_train_with_ID = self.sets_for_other_gen_tran_and_test_with_ID\
                            (input_id_gen_train_sim,data_target)
            input_test_with_ID = self.sets_for_other_gen_tran_and_test_with_ID\
                            (input_id_test_sim,data_target)

            input_gener_train = self.sets_for_other_gen_tran_and_test\
                            (input_id_gen_train_sim,data_target)
            input_test = self.sets_for_other_gen_tran_and_test\
                            (input_id_test_sim,data_target)

            input_gener_train_with_p_chembl = self.sets_for_other_gen_tran_and_test_with_pchembl\
                            (input_id_gen_train_sim,data_target)
            input_test_with_p_chembl = self.sets_for_other_gen_tran_and_test_with_pchembl\
                            (input_id_test_sim,data_target)
            
            '''Option to save'''
            if self.save_options == True:
                print("Save options")
                #Input Set Molpher
                input_Molpher_train.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_Molpher_{self.receptor_name}_sim_{x[1]}.csv", index=False)

                #Input Set other generators
                input_gener_train.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_{self.receptor_name}_sim_{x[1]}.csv", index=False,header=False)
                input_gener_train_with_ID.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_{self.receptor_name}_with_ID_sim_{x[1]}.csv", index=False)
                input_gener_train_with_p_chembl.to_csv(f"data/input_recall_sets/{self.receptor_name}/cIS_{self.receptor_name}_with_p_chembl_sim_{x[1]}.csv", index=False,header=False)
                
                #Recall Set
                input_test.to_csv(f"data/input_recall_sets/{self.receptor_name}/cRS_{self.receptor_name}_sim_{x[1]}.csv",  index=False,header=False)
                input_test_with_ID.to_csv(f"data/input_recall_sets/{self.receptor_name}/cRS_{self.receptor_name}_with_ID_sim_{x[1]}.csv",  index=False)
                input_test_with_p_chembl.to_csv(f"data/input_recall_sets/{self.receptor_name}/cRS_{self.receptor_name}_with_p_chembl_sim_{x[1]}.csv",  index=False,header=False)