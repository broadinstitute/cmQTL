#!/usr/bin/env python
# coding: utf-8

# # Apply Cell Health models to cmQTL morphology profiles
# 
# Previously, we trained machine learning models to predict 70 different cell health phenotypes from Cell Painting morphology profiles ([Way et al. 2020](https://doi.org/10.1101/2020.07.08.193938)). Here, we apply three of the models to well-level aggregate profiles to separate data from the cmQTL project.
# 
# ## Cell Health models
# 
# The three models include:
# 
# 1. ALL - # cells (`cc_all_n_objects`)
# 2. ALL - Nucleus Area (`cc_all_nucleus_area_mean`)
# 3. % Dead Only (CASP-; DRAQ7+) (`vb_percent_dead_only`)
# 
# ## Input
# 
# Input to this notebook are normalized morphology profiles (prior to feature selection) for the cmQTL dataset.
# 
# ## Output
# 
# Predictions for the three cell health models per profile.
# 
# ## Discussion
# 
# See https://github.com/broadinstitute/cmQTL/issues/47 for more details.
# 
# ## WARNING
# 
# Currently, no predictions are output - we are working towards resolving the missing feature problem.
# There are 506 features measured in the Cell Health project not measured in cmQTL.

# In[1]:


import pathlib
import pandas as pd
from joblib import load

from pycytominer.cyto_utils import infer_cp_features


# ## Load and process data

# In[2]:


profile_file = pathlib.Path("data/cell-health/profiles/anycells.well.forCellHealth.tab")
profile_df = pd.read_csv(profile_file, sep="\t")

print(profile_df.shape)
profile_df.head(2)


# In[3]:


get_ipython().system(' md5 $profile_file')


# In[4]:


# Load dataset used for training to get order of coefficients
commit = "07e4b40c39dd27084be36fbef4d64c5654b2960f"
data_url = f"https://github.com/broadinstitute/cell-health/raw/{commit}/3.train/data/x_train_modz.tsv.gz"
cell_health_df = pd.read_csv(data_url, sep="\t")

cell_health_feature_order = infer_cp_features(cell_health_df)

print(cell_health_df.shape)
cell_health_df.head(2)


# In[5]:


cp_cols = infer_cp_features(profile_df)
metadata_cols = infer_cp_features(profile_df, metadata=True)

# Reindex the cmQTL profiles
profile_reindexed_df = profile_df.reindex(columns=metadata_cols + cell_health_feature_order)


# In[6]:


features_missing = set(cell_health_feature_order).difference(set(cp_cols))

print(f"There are a total of {len(cell_health_feature_order)} features in the Cell Health project.")
print(f"There are {len(features_missing)} features in the Cell Health project missing from the cmQTL project.")


# ## Load models

# In[7]:


model_dir = pathlib.Path("data/cell-health/models/")
model_dict = {
    "cc_all_n_objects": {"name": "ALL - # cells"},
    "cc_all_nucleus_area_mean": {"name": "ALL - Nucleus Area"},
    "vb_percent_dead_only": {"name": "% Dead Only (CASP-; DRAQ7+)"}
}


# In[8]:


for model in model_dict.keys():
    model_file = model_dir / pathlib.Path(f"cell_health_modz_target_{model}_shuffle_False_transform_raw.joblib")
    model_job = load(model_file)
    
    model_coef_df = (
        pd.DataFrame(
            model_job.coef_, columns=[model], index=cell_health_feature_order
        )
        .reset_index()
        .rename({"index": "cp_features"}, axis="columns")
        .assign(missing=False)
    )
    model_coef_df.loc[model_coef_df.cp_features.isin(features_missing), "missing"] = True
    model_coef_df = model_coef_df.assign(abs_coef=model_coef_df.loc[:, model].abs())
    
    model_dict[model]["model"] = model_job
    model_dict[model]["coef"] = model_coef_df
    
    # Check the total contribution lost by missing features
    print(model_coef_df.groupby(["missing"])["abs_coef"].sum())


# In[9]:


features_missing

