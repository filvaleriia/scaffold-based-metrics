import pandas as pd
import argparse
from sklearn.preprocessing import MinMaxScaler



def connect_mean_value(type_cluster, type_scaffold, generators_name_list, receptor, subset):
    """
    Combines mean values from multiple generators into a single DataFrame.
    
    Parameters:
    - receptor: Target receptor
    - type_scaffold: Type of scaffold
    - type_cluster: Cluster type
    - generators_name_list: List of generator names
    - subset: Data subset to be analyzed

    Returns:
    - A combined DataFrame with mean values
    """
    
    # Define path to data
    link = f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"
    
    # List to store paths of mean value files
    link_mean = []
    for generator in generators_name_list:
        link_mean.append(f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv")

    # Load data and merge into a single DataFrame
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, ignore_index=True)

    # Save results to CSV files
    df.to_csv(f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}/mean_{type_scaffold}_{type_cluster}{subset}.csv", index=False)

    return df


def connect_mean_value_normalized(type_cluster, type_scaffold, generators_name_list, receptor, subset):
    """
    Loads and normalizes mean values using Min-Max scaling.
    
    Parameters:
    - receptor: Target receptor
    - type_scaffold: Type of scaffold
    - type_cluster: Cluster type
    - generators_name_list: List of generator names
    - subset: Data subset to be analyzed

    Returns:
    - Normalized DataFrame
    """

    # Define path to data
    link = f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}"

    # List to store paths of mean value files
    link_mean = [f"{link}/{generator}/{generator}_mean_{type_scaffold}_{type_cluster}.csv" for generator in generators_name_list]

    # Load data
    df_list = [pd.read_csv(f) for f in link_mean]
    df = pd.concat(df_list, ignore_index=True)

    # Normalize using Min-Max scaling
    scaler = MinMaxScaler()
    numeric_columns = df.select_dtypes(include=['number']).columns  # Select only numeric columns
    df[numeric_columns] = scaler.fit_transform(df[numeric_columns])  # Apply normalization

    # Save normalized results
    df.to_csv(f"data/results/{receptor}/{type_scaffold}_scaffolds/{type_cluster}/mean_{type_scaffold}_{type_cluster}{subset}_norm_min_max.csv", index=False)

    return df


def main():
    parser = argparse.ArgumentParser(description='Compute and visualize recall metric.')
     # Required arguments
    parser.add_argument('--type_cluster', type=str, required=True, help='Type of clustering (dis/sim)')
    parser.add_argument('--type_scaffold', type=str, required=True, help='Type of scaffold')
    parser.add_argument('--generator_list', type="+", required=True, help='Generator name')
    parser.add_argument('--receptor', type=str, required=True, help='Receptor name')
    parser.add_argument('--subset', type=str, required=True, help='Subset')
    
    args = parser.parse_args()
    
    connect_mean_value(args.type_cluster, args.type_scaffold, args.generator_list, args.receptor, args.subset)


if __name__ == "__main__":
    main()