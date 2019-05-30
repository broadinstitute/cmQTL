#!/usr/bin/env python
# coding: utf-8

# # Process Single Cell Profiles
# 
# **Gregory Way, 2019**
# 
# Load in specific columns:
# 
#   * `Cells_Neighbors_PercentTouching_Adjacent`
#   * `Cells_Neighbors_NumberOfNeighbors_Adjacent` 
#   
# to visualize distributions across batches and cell lines.
# 
# ## SQLite Tables
# 
# There are 4 tables in the database:
# 
# 1. Cells
# 2. Cytoplasm
# 3. Nuclei
# 4. Image
# 
# They are to be merged together.
# 
# The first three can be merged by the following columns: `["TableNumber", "ImageNumber", "ObjectNumber"]`
# 
# The `Image` table can be merged by `["TableNumber", "ImageNumber"]`

# In[1]:


import os
import glob
import pandas as pd
import sqlite3

from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects

import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


def get_sc_info(batch_id):
    conn = sqlite3.connect(sql_file_dict[batch_id])
    
    sc_cells_df = pd.read_sql_query(
        """
        select
        Cells_Neighbors_PercentTouching_Adjacent,
        Cells_Neighbors_NumberOfNeighbors_Adjacent
        from Cells
        """, conn)
    
    sc_cells_df.Cells_Neighbors_NumberOfNeighbors_Adjacent.hist(bins=12)
    plt.xlabel("NumberOfNeighbors_Adjacent - {}".format(batch_id))
    
    file = os.path.join("figures",
                        "single_cell_diagnostics",
                        "{}_adjacent.png".format(batch_id))
    plt.savefig(file, height = 4, width = 5)
    plt.show()
    plt.close()
    
    sc_cells_df.Cells_Neighbors_PercentTouching_Adjacent.hist(bins=12)
    plt.xlabel("NumberOfNeighbors_PercentTouching - {}".format(batch_id))
    
    file = os.path.join("figures",
                        "single_cell_diagnostics",
                        "{}_percenttouching.png".format(batch_id))
    plt.savefig(file, height = 4, width = 5)
    plt.show()
    plt.close()
    
    # Seaborn Joint Plot
    g = sns.jointplot("Cells_Neighbors_PercentTouching_Adjacent",
                      "Cells_Neighbors_NumberOfNeighbors_Adjacent",
                      data=sc_cells_df,
                      kind="hex")
    
    file = os.path.join("figures",
                        "single_cell_diagnostics",
                        "{}_neighbor_relationship.png".format(batch_id))
    plt.savefig(file, height = 4, width = 5)
    plt.show()
    plt.close()


# In[4]:


os.makedirs(os.path.join("figures", "single_cell_diagnostics"), exist_ok=True)
os.makedirs("data", exist_ok=True)
os.makedirs("results", exist_ok=True)


# ## Step 0 - Load Constants and Generate Objects

# In[5]:


batches = ["BR00103267", "BR00103268"]

cell_cols = ["TableNumber", "ImageNumber", "ObjectNumber"]
image_cols = ["TableNumber", "ImageNumber"]

bucket_dir = os.path.join("..", "..", "..", "..", "..", "..", "bucket")
backend_dir = os.path.join("projects", "2018_06_05_cmQTL", "workspace", "backend")
batch_dir = os.path.join(bucket_dir, backend_dir, "2019_05_13_Batch2")

sql_file_dict = {
    x: os.path.join(batch_dir, x, "{}.sqlite".format(x))
    for x in batches
}

sql_file_dict


# In[6]:


# Load example data for selected features
file = os.path.join(batch_dir, batches[0], 'BR00103267_normalized_variable_selected.csv')
example_df = pd.read_csv(file, sep=',')

cp_features = example_df.loc[:, ~example_df.columns.str.startswith("Metadata_")].columns.tolist()
print(len(cp_features))

example_df.head(2)


# In[7]:


use_cols = cell_cols + ["Metadata_Well", "Metadata_Plate"] + cp_features
len(use_cols)


# ## Step 1 - Visualize Single Cell Relationships

# In[8]:


# Generate and save plots
for batch in batches:
    get_sc_info(batch)


# ## Step 2 - Count how many cells per well and single cell type

# In[9]:


cell_count_list = []
for batch_id in batches:
    conn = sqlite3.connect(sql_file_dict[batch_id])
    
    # Read in single cell profiles
    sc_cells_df = pd.read_sql_query(
        """
        select
        TableNumber,
        ImageNumber,
        ObjectNumber,
        Cells_Neighbors_NumberOfNeighbors_Adjacent
        from Cells
        """, conn).merge(
        pd.read_sql_query(
            """
            select
            TableNumber,
            ImageNumber,
            Metadata_Well
            from Image
            """, conn),
        left_on=image_cols,
            right_on=image_cols,
            how='left'
        )

    # Get Cell Counts
    total_count_df = pd.DataFrame(
        sc_cells_df.groupby("Metadata_Well")
        ['ObjectNumber']
        .count()
        .reset_index()
        .rename({"ObjectNumber": "cell_count"}, axis='columns')
        .assign(batch_id=batch_id,
                sc_type="all")
    )
    cell_count_list.append(total_count_df)
    
    isolated_count_df = pd.DataFrame(
        sc_cells_df
        .query("Cells_Neighbors_NumberOfNeighbors_Adjacent == 0")
        .groupby("Metadata_Well")
        ['ObjectNumber']
        .count()
        .reset_index()
        .rename({"ObjectNumber": "cell_count"}, axis='columns')
        .assign(batch_id=batch_id,
                sc_type="isolated")
    )
    cell_count_list.append(isolated_count_df)
    
    colony_count_df = pd.DataFrame(
        sc_cells_df
        .query("Cells_Neighbors_NumberOfNeighbors_Adjacent >= 4")
        .groupby("Metadata_Well")
        ['ObjectNumber']
        .count()
        .reset_index()
        .rename({"ObjectNumber": "cell_count"}, axis='columns')
        .assign(batch_id=batch_id,
                sc_type="colony")
    )
    cell_count_list.append(colony_count_df)


# In[10]:


cell_count_df = pd.concat(cell_count_list)

file = os.path.join("results", "well_cell_counts.tsv")
cell_count_df.to_csv(file, sep='\t')


# ## Step 3 - Extract Single Cell Profiles

# In[11]:


for batch_id in batches:
    conn = sqlite3.connect(sql_file_dict[batch_id])
    
    # Read in single cell profiles
    sc_cells_df = pd.read_sql_query("""
        select *
        from Cells
        where Cells_Neighbors_NumberOfNeighbors_Adjacent == 0
        """, conn)

    # Load and Merge Single Cell Data
    sc_cells_df = (
        sc_cells_df
        .merge(
            pd.read_sql_query("""
                select *
                from Cytoplasm
                """, conn),
            left_on=cell_cols,
            right_on=cell_cols,
            how='left'
        )
        .merge(
            pd.read_sql_query("""
                select *
                from Nuclei
                """, conn),
            left_on=cell_cols,
            right_on=cell_cols,
            how='left'
        )
        .merge(
            pd.read_sql_query("""
                select *
                from Image
                """, conn),
            left_on=image_cols,
            right_on=image_cols,
            how='left'
        )
    ).loc[:, use_cols]

    # Write out the file to disk
    file = os.path.join("data", "{}_single_cell_isolated_profiles.tsv.gz".format(batch_id))
    sc_cells_df.to_csv(file, sep='\t', index=False)
    
    del sc_cells_df
    
    colony_cells_df = pd.read_sql_query("""
        select *
        from Cells
        where Cells_Neighbors_NumberOfNeighbors_Adjacent >= 4
        """, conn)

    # Load Single Cell Data
    colony_cells_df = (
        colony_cells_df.merge(
            pd.read_sql_query("""
                select *
                from Cytoplasm
                """, conn),
            left_on=cell_cols,
            right_on=cell_cols,
            how='left'
        ).merge(
            pd.read_sql_query("""
                select *
                from Nuclei
                """, conn),
            left_on=cell_cols,
            right_on=cell_cols,
            how='left'
        ).merge(
            pd.read_sql_query("""
                select *
                from Image
                """, conn),
            left_on=image_cols,
            right_on=image_cols,
            how='left'
        )
    ).loc[:, use_cols]
    
    file = os.path.join("data", "{}_single_cell_colony_profiles.tsv.gz".format(batch_id))
    colony_cells_df.to_csv(file, sep='\t', index=False)
    
    del colony_cells_df

