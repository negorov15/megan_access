import pandas as pd
from io import StringIO
from get_lineage import get_lineage
from functools import reduce


def merge_data(list_of_datasets):
    """
    Merges the given datasets into a single dataset.
    The merge is done by the column with Taxa IDs.
    The output is saved in a csv file.
    """
    df_merged = reduce(
        lambda left, right: pd.merge(left, right, on=["Taxa"], how="outer"),
        list_of_datasets,
    )
    # Save the output in csv file
    pd.DataFrame.to_csv(df_merged, "merged.csv", sep="\t", index=False)

    print("Datasets were successfully merged!")


def dataset_modifier(input, name):
    """
    Modifies the given dataset by skipping the first 2 rows and deleting the column with weights assigned to the id.
    The output is a DataFrame with two columns: Taxa and the name of the sample.
    """
    # Read the response text into a DataFrame using pandas. StringIO allows pandas to treat input like file object, not file path.
    df = pd.read_csv(StringIO(input), skiprows=2, delimiter="\t", encoding="utf-8")

    # Index of the column with float numbers of reads
    column_index_to_delete = 2

    if column_index_to_delete < len(df.columns):
        # Drop the specified column
        df = df.drop(df.columns[column_index_to_delete], axis=1)
        # Set the names of the columns
        df.columns = ["Taxa", name]
        # Rewrite the actual file. Is not needed for now.
        # df.to_csv(input, sep='\t', index=False)
        return df


def otu_table(comparison_file):
    """
    Modifies the given dataset by deleting the column with Taxa IDs and setting the index of the dataframe.
    The output is a DataFrame with the index as OTU IDs and the columns as sample names.
    """
    # Read the comparison csv file
    merged_df = pd.read_csv(comparison_file, sep="\t")
    # Index of the column with Taxa IDs
    column_index_to_delete = 0
    # Drop the column with Taxa IDs
    if column_index_to_delete < len(merged_df.columns):
        merged_df = merged_df.drop(merged_df.columns[column_index_to_delete], axis=1)
    # Get the number of rows in the dataset
    num_of_rows = len(merged_df)
    # Set the index of the dataframe
    new_index = ["OTU" + str(i + 1) for i in range(num_of_rows)]
    merged_df.index = new_index
    return merged_df


def tax_table(dataframe):
    """
    Modifies the given dataset by deleting the column with float numbers of reads and setting the index of the dataframe.
    The output is a DataFrame with the index as OTU IDs and the columns as taxonomic information.
    """
    # Get the Dataframe with lineages
    df = get_lineage(dataframe)
    # Set the index for the dataframe
    new_index = ["OTU" + str(i + 1) for i in range(len(df))]
    df.index = new_index
    return df
