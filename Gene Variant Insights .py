import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the dataset
file_path_longevity = r'C:\Users\IoannisZografakis-Re\Downloads\3-longevity.csv'  # Adjust the file path accordingly
longevity_df = pd.read_csv(file_path_longevity)

# Step 1: Filter the Data for significant variants studied in multiple populations
significant_variants_df = longevity_df[longevity_df['Association'] == 'significant']
variant_counts = significant_variants_df.groupby('Variant')['Population'].nunique()
variants_with_multiple_populations = variant_counts[variant_counts >= 2].index

filtered_df = significant_variants_df[significant_variants_df['Variant'].isin(variants_with_multiple_populations)]

# Selecting specific columns
filtered_df = filtered_df[['id', 'Association', 'Population', 'Variant', 'Gene', 'PubMed']]

# Export the filtered data
output_file_path_longevity = r'C:\Users\IoannisZografakis-Re\Documents\filtered_longevity.csv'  # Adjust the path
filtered_df.to_csv(output_file_path_longevity, index=False)

# Step 2: Data Cleaning - Handle missing values in 'Variant' column
filtered_df['Variant'] = filtered_df['Variant'].fillna('Unknown')

# Step 3: Aggregation - Count significant variants per population
population_variant_counts = filtered_df.groupby('Population')['Variant'].nunique().reset_index(name='Variant_Count')

# Step 4: Visualization - Horizontal bar chart for variants per population
plt.figure(figsize=(10, 6))
colors = plt.get_cmap('tab20')(range(len(population_variant_counts)))

plt.barh(population_variant_counts['Population'], population_variant_counts['Variant_Count'], color=colors)
plt.title('Number of Significant Genetic Variants per Population')
plt.xlabel('Number of Variants')
plt.ylabel('Population')

plt.tight_layout()
plt.show()

# Step 5: Data Cleaning - Fill missing PubMed IDs with 0
filtered_df['PubMed'] = filtered_df['PubMed'].fillna(0).astype(int)

# Step 6: Create an artificial timeline from 2000 to 2023
num_records = len(filtered_df)
years = np.linspace(2000, 2023, num_records).astype(int)
filtered_df['Year'] = years

# Step 7: Aggregation - Calculate number of variants per year and 5-year rolling average
yearly_variant_counts = filtered_df.groupby('Year').size().reset_index(name='Variant_Count')
yearly_variant_counts['Rolling_Avg'] = yearly_variant_counts['Variant_Count'].rolling(window=5).mean()

# Step 8: Visualization - Line plot for annual counts and 5-year rolling average
plt.figure(figsize=(10, 6))

plt.plot(yearly_variant_counts['Year'], yearly_variant_counts['Variant_Count'], label='Annual Count', marker='o')
plt.plot(yearly_variant_counts['Year'], yearly_variant_counts['Rolling_Avg'], label='5-Year Rolling Avg', linestyle='--')

plt.title('Annual Count of Significant Genetic Variants and 5-Year Rolling Average')
plt.xlabel('Year')
plt.ylabel('Number of Significant Variants')
plt.legend()

plt.tight_layout()
plt.show()

# Step 9: Data Cleaning - Fill missing 'Variant' and 'Gene' with 'Unknown'
filtered_df['Variant'] = filtered_df['Variant'].fillna('Unknown')
filtered_df['Gene'] = filtered_df['Gene'].fillna('Unknown')

# Step 10: Aggregation - Calculate number of variants and unique genes per population
population_variant_gene_counts = filtered_df.groupby('Population').agg(
    num_variants=pd.NamedAgg(column='Variant', aggfunc='nunique'),
    num_unique_genes=pd.NamedAgg(column='Gene', aggfunc='nunique')
).reset_index()

# Step 11: Correlation Calculation between number of variants and unique genes
correlation = np.corrcoef(population_variant_gene_counts['num_variants'], population_variant_gene_counts['num_unique_genes'])[0, 1]

# Step 12: Visualization - Scatter plot to show relationship between variants and unique genes with 'o' symbol and different colors
plt.figure(figsize=(8, 6))

# Assign different colors to each point (bullet)
scatter_colors = sns.color_palette('husl', len(population_variant_gene_counts))

# Scatter plot with unique colors for each bullet (point)
for i, row in population_variant_gene_counts.iterrows():
    plt.scatter(row['num_variants'], row['num_unique_genes'], color=scatter_colors[i], s=100, marker='o')

plt.title(f'Relationship Between Significant Variants and Unique Genes (Correlation: {correlation:.2f})')
plt.xlabel('Number of Significant Variants')
plt.ylabel('Number of Unique Genes')

plt.tight_layout()
plt.show()
