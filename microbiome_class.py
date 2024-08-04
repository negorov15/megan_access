import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from skbio.stats.distance import permanova
from get_lineage import get_lineage
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
from modificator import otu_table, tax_table

# Class MicrobiomeDataAnalyzer.
# Contains: OTU table, Taxa table and corresponding metadata
class MicrobiomeDataAnalyzer:
    def __init__(self, OTU_table: str = "", Taxa_table: str = "",
                 Metadata: str = "") -> None:
        self.OTU_table=OTU_table
        self.Taxa_table=Taxa_table
        self.Metadata=Metadata

    # Data manipulation functions:
    # 1. prepare_dataset(self)
    # 2. abundance_table()
    # Helper-function. Transposes the given dataset.
    def prepare_dataset(self):
        # Transpose the OTU dataframe
        otumat_transposed = self.OTU_table.T

        # Reset indices in transposed OTU table
        otumat_transposed = otumat_transposed.reset_index()
        otumat_transposed = otumat_transposed.rename(columns={'index': 'SampleID'})

        # If the metadata includes multiple properties, you can delete them and use only 'Group' column
        filtered_metadata = self.Metadata.drop("Property", axis='columns')
        # Merge the OTU table with the metadata
        # input_data = pd.merge(otumat_transposed, self.Metadata, on='SampleID')
        input_data = pd.merge(otumat_transposed, filtered_metadata, on='SampleID')

        # In this step we move the column with group property right after the column with SampleID.
        column_to_move = input_data.pop("Group")

        # Insert column with insert(location, column_name, column_value)
        input_data.insert(1, "Group", column_to_move)

        return input_data

    # Create an abundance table.
    def abundance_table(self):
        # New DataFrame with relative abundances
        abundance_table = self.OTU_table.copy()
        for column in abundance_table:
            sum = abundance_table[column].sum()
            abundance_table[column] = abundance_table[column].apply(lambda x: x / sum)
        return abundance_table

    # Data visualization functions
    # 1. plot_rank()
    # 2. plot_top()
    # 3. plot_pcoa()
    # Plot the dataset by a certain rank
    def plot_rank(self, rank):
        taxa_mapping = self.Taxa_table[rank]
        aggregated = self.OTU_table.groupby(taxa_mapping, axis=0).sum()

        # Transpose the aggregated DataFrame to have taxa on the x-axis
        df_t = aggregated.T

        # Create a figure and axes with a specified size
        fig, ax = plt.subplots(figsize=(20, 15))

        # Plot the data with larger font sizes
        df_t.plot(kind='bar', stacked=True, ax=ax)

        ax.set_ylabel('Abundance', fontsize=40)  # Increase font size for y-axis label
        ax.tick_params(axis='x', labelsize=28, rotation=45)  # Increase font size for x-axis ticks
        ax.tick_params(axis='y', labelsize=28)  # Increase font size for y-axis ticks
        ax.legend(title=rank, fontsize='20', title_fontsize='25', loc='upper left', bbox_to_anchor=(1.05, 1),
                  borderaxespad=0.)  # Adjust legend font size

        plt.tight_layout()  # Adjust the layout
        plt.show()

    # Plot top OTUs from the samples
    def plot_top(self, top):
        # Compute total read counts for each OTU
        total_counts = self.OTU_table.sum(axis=1)

        # Identify the top OTUs
        top_otus = total_counts.nlargest(top).index

        # Extract the most specific taxonomic rank for the top OTUs from the taxonomy table
        top_otus_tax_rank = {}
        for otu in top_otus:
            tax_info = self.Taxa_table.loc[otu]  # Access the row corresponding to the OTU
            # Filter out empty or NaN values and get the last available taxonomic rank
            specific_rank = tax_info.dropna().iloc[-1] if not tax_info.dropna().empty else 'Unknown'
            top_otus_tax_rank[otu] = specific_rank

        # Filter the DataFrame for these top OTUs
        filtered_df = self.OTU_table.loc[top_otus]

        # Transpose the filtered DataFrame for plotting
        df_t = filtered_df.T

        # Create a figure and axes with a specified size
        fig, ax = plt.subplots(figsize=(20, 15))

        # Plot the data
        df_t.plot(kind='bar', stacked=True, ax=ax)

        # Update plot titles and labels with larger font sizes
        ax.set_title(f'Top {top} taxa by read count', fontsize=45)
        ax.set_ylabel('Read Count', fontsize=40)
        ax.tick_params(axis='x', labelsize=28, rotation=45)
        ax.tick_params(axis='y', labelsize=28)

        # Update the legend to use the most specific available taxonomic rank instead of OTU identifiers
        legend_labels = [top_otus_tax_rank[otu] for otu in top_otus]
        ax.legend(legend_labels, fontsize='20', title_fontsize='25', loc='upper left', bbox_to_anchor=(1.05, 1),
                  borderaxespad=0.)

        plt.tight_layout()
        plt.show()

    def plot_pcoa(self):
        # Compute beta diversity distance
        distance = self.beta_diversity()
        # Perform PCoA on the Bray-Curtis distance matrix
        pcoa_results = pcoa(distance)

        pcoa_samples = pcoa_results.samples

        groups = np.array(self.Metadata['Group'])
        treatment_day = np.array(self.Metadata['Property'])
        pcoa_samples.insert(0, "Group", groups)
        pcoa_samples.insert(1, "Property", treatment_day)

        # Set up the color palette
        unique_groups = pcoa_samples['Group'].unique()
        colors = plt.get_cmap('tab10')(range(len(unique_groups)))
        color_map = dict(zip(unique_groups, colors))

        # Initialize an empty set to keep track of groups already added to the legend
        seen_groups = set()

        # Plot each sample and annotate it with the sample ID
        for idx, row in pcoa_samples.iterrows():
            group = row['Group']
            pc1 = row['PC1']
            pc2 = row['PC2']
            property = row['Property']

            if group not in seen_groups:
                plt.scatter(pc1, pc2, color=color_map[group], label=group,
                            alpha=0.7)  # Include label only if group is unseen
                seen_groups.add(group)  # Mark the group as seen
            else:
                plt.scatter(pc1, pc2, color=color_map[group], alpha=0.7)  # No label for already seen groups

            plt.text(pc1, pc2, property, fontsize=9)

        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('PCoA of Beta Diversity')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.show()

    # Statistical analysis functions:
    # 1. t_test()
    # 2. anova_test()
    # 3. kruskal()
    # 4. permanova()
    # 5. beta_diversity()
    # Perform an independent t-Test
    def t_test(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()

        # Find all unique color groups in the dataset.
        unique_groups = input_data['Group'].unique()
        # Save the possible colors
        column_property = input_data['Group']

        # Perform t-Test for each OTU column
        t_test_results = {}
        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Property'

            # Filter the dataframe by color and drop missing values.
            groups = [input_data[column_property == property][otu].dropna() for property in unique_groups]

            t_stat, p_value = stats.ttest_ind(*groups)

            # Store the results
            t_test_results[otu] = (t_stat, p_value)

        # Display the results
        for otu, (t_stat, p_value) in t_test_results.items():
            print(f"{otu}: T statistic = {t_stat}, P-value = {p_value}")

    # Perform Anova-Test
    def anova_test(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()
        # Prepare groups for ANOVA
        # Find all unique groups in the dataset.
        unique_groups = input_data['Group'].unique()
        # Save the possible properties (for example: fur color of the subject)
        column_property = input_data['Group']

        # Perform ANOVA for each OTU column
        anova_results = {}
        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Property'

            # Filter the dataframe by color and drop missing values.
            groups = [input_data[column_property == property][otu].dropna() for property in unique_groups]

            f_value, p_value = stats.f_oneway(*groups)

            # Store the results
            anova_results[otu] = (f_value, p_value)

        # Display the results
        for otu, (f_value, p_value) in anova_results.items():
            print(f"{otu}: F-value = {f_value}, P-value = {p_value}")

    # Perform Kruskal test. Alternative to the ANOVA-test. Data does not have to be normally distributed.
    def kruskal(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()
        # Find all unique groups in the dataset.
        unique_groups = input_data['Group'].unique()
        # Save the possible properties (for example: fur color of the subject)
        column_property = input_data['Group']

        # Perform Kruskal for each OTU column
        kruskal_results = {}
        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Property'

            # Filter the dataframe by color and drop missing values.
            groups = [input_data[column_property == property][otu].dropna() for property in unique_groups]

            stat, p_value = stats.kruskal(*groups)

            # Store the results
            kruskal_results[otu] = (stat, p_value)

        # Display the results
        for otu, (stat, p_value) in kruskal_results.items():
            print(f"{otu}: H-statistic = {stat}, P-value = {p_value}")

    def wilcoxon_test(self):
        # Prepare the dataset
        input_data = self.prepare_dataset()

        # Assuming two unique groups represent paired data (e.g., "before" and "after")
        if len(input_data['Group'].unique()) != 2:
            raise ValueError("There must be exactly two groups for the Wilcoxon signed-rank test.")

        # Retrieve the unique groups for clarity
        group1, group2 = input_data['Group'].unique()

        # Perform Wilcoxon test for each OTU column
        wilcoxon_results = {}

        for otu in input_data.columns[2:]:  # Skipping 'SampleID' and 'Group'

            # Filter the dataframe by group and drop missing values
            group1_data = input_data[input_data['Group'] == group1][otu].dropna()
            group2_data = input_data[input_data['Group'] == group2][otu].dropna()

            # Ensure equal lengths before performing the test
            min_len = min(len(group1_data), len(group2_data))
            if min_len > 0:
                stat, p_value = wilcoxon(group1_data[:min_len], group2_data[:min_len])

                # Store the results
                wilcoxon_results[otu] = (stat, p_value)
            else:
                wilcoxon_results[otu] = (None, None)

        # Display the results
        for otu, (stat, p_value) in wilcoxon_results.items():
            print(f"{otu}: Wilcoxon statistic = {stat}, P-value = {p_value}")

    # Beta diversity
    def beta_diversity(self):
        # Drop na values
        dropped_df = self.OTU_table.dropna()
        # Transpsoe the OTU dataframe to make samples as rows
        otumat_transposed = dropped_df.T

        data = otumat_transposed.values
        ids = otumat_transposed.index

        # Compute Bray-Curtis dissimilarity
        b_div = beta_diversity(metric="braycurtis", counts=data, ids=ids)

        return b_div

    def permanova(self):
        distance = self.beta_diversity()
        grouping = self.Metadata['Group']
        permanova_results = permanova(distance, grouping)
        return permanova_results

# Example data
if __name__ == '__main__':

    # Object Alice
    otumat_alice = otu_table('alice.csv')
    taxmat_alice = tax_table('alice.csv')
    # Example metadata
    metadata_dict_alice = {
        'SampleID': ['Alice00-1mio.daa', 'Alice01-1mio.daa',
                     'Alice03-1mio.daa', 'Alice06-1mio.daa',
                     'Alice08-1mio.daa', 'Alice34-1mio.daa'],
        'Group': ['Without treatment', 'With treatment', 'With treatment',
                  'With treatment', 'Without treatment', 'Without treatment'],
        'Property': ['0-', '1+', '3+', '6+', '8-', '34-']
    }
    metadata_alice = pd.DataFrame(metadata_dict_alice)
    sample_alice = MicrobiomeDataAnalyzer(otumat_alice,taxmat_alice,metadata_alice)

    # Object Bob
    otumat_bob = otu_table('bob.csv')
    taxmat_bob = tax_table('bob.csv')
    metadata_dict_bob = {
        'SampleID': ['Bob00-1mio.daa', 'Bob01-1mio.daa',
                     'Bob03-1mio.daa', 'Bob06-1mio.daa',
                     'Bob08-1mio.daa', 'Bob34-1mio.daa'],
        'Group': ['Without treatment', 'With treatment', 'With treatment',
                  'With treatment', 'Without treatment', 'Without treatment'],
        'Property': ['0-', '1+', '3+', '6+', '8-', '34-']
    }
    metadata_bob = pd.DataFrame(metadata_dict_bob)
    sample_bob = MicrobiomeDataAnalyzer(otumat_bob, taxmat_bob, metadata_bob)

    print(sample_alice.permanova())
    print(sample_bob.permanova())
    print(sample_alice.t_test())
    print(sample_bob.t_test())
    sample_bob.plot_rank('Phylum')
    sample_bob.plot_top(5)