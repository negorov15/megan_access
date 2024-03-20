import pandas as pd
from ete3 import NCBITaxa

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

# Function gets a lineage for each taxa ID for tax matrix.
# Input: pandas DataFrame dataset
# Output: list of lineages
def get_lineage(dataset):
    # Read the input file
    df = pd.read_csv(dataset, delimiter="\t")
    list_of_lineages = []
    # Extract taxa ids and loop through the list
    tax_ids = df['Taxa']
    for id in tax_ids:
        # -2 = 'Not assigned'.
        # We just add 'not assigned' rank also to the list to avoid misunderstandings
        if id == -2:
            list_of_lineages.append('Not assigned')
            continue
        else:
            try:
                list_of_lineages.append(ncbi.get_lineage(id))
            except ValueError:
                continue
    # Call 'sort_ranks' function and store it to the variable 'sorted_lineages'
    sorted_lineages = sort_ranks(list_of_lineages)
    return translate_lineage(sorted_lineages)

# Helper function. Sorts ranks so that only major ranks left.
# Input: List of lineages
# Output: List of sorted lineages
def sort_ranks(data):
    major_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    sorted_ranks = []
    for lineage in data:
        # Temporary list for storing ranks.
        # Is required, so the final list would have separate lists of ranks associated with each element
        lineage_ranks = []
        # Check if the element of list has an actual lineage.
        # If it is the case, just add 'Not assigned' to the list
        if lineage == "Not assigned":
            sorted_ranks.append([lineage])
        else:
            # Get ranks for each lineage
            rank = ncbi.get_rank(lineage)
            for i in range(len(lineage)):
                # If the rank is not in the list of major ranks, we skip it
                if rank[lineage[i]] in major_ranks:
                    lineage_ranks.append(lineage[i])
                else:
                    continue
            # add the temporary list to the main list
            sorted_ranks.append((lineage_ranks))
    return sorted_ranks

# Function translates the lineages into names.
# Input: List of lineages
# Output: Translated list of lineages
def translate_lineage(list_of_lineages):
    # Empty list where all names will be stored
    names = []
    # Loop through each lineage list and get the names of each rank
    for lineage in list_of_lineages:
        # If the rank is 'not assigned' simply add the value to the list.
        if lineage == ["Not assigned"]:
            names.append(lineage)
        else:
            # get_taxid_translator(lineage) creates a dictionary,
            # so we need to get all values and store it to the list
            translate = ncbi.get_taxid_translator(lineage)
            names.append([translate[taxid] for taxid in lineage])
    # Create a DataFrame where:
    # Rows: taxa
    # Columns: taxonomy lineage associated with the taxa
    df = pd.DataFrame(names, columns=["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"])
    return df