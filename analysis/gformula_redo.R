### Setup
library(data.table)

### Paths
data.path <- "data/bone_data.csv"

### Code

## Step 1 - generate person-day data from bone marrow transplant data
# Read in data
dt <- fread(data.path)

## Step 2 - estimate modeling coefficients used to generate probabilities

## Step 3 - sample with replacement from data

## Step 4 and 5 - run Monte Carlo sample for natural course, always and never exposed

## Step 6 - concatentate intervetion data sets and run Cox model
