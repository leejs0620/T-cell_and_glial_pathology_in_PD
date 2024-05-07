# T-cell_and_glial_pathology_in_PD
Code and data for The spatial landscape of glial pathology and adaptive immune response in Parkinson's Disease 

# Spatial cross-correlation
This is a Python module for spatial cross-correlation analyses

## Installation
```bash
$ git clone https://github.com/dalhoomist/T-cell_and_glial_pathology_in_PD.git
```

## Set up Environment
Requirement
```bash
Python: 3.6.13
Numpy: 1.18.5
Pandas: 1.1.5
```
[The GPU-accelerated computing module is currently not shareable. We aim to upload it as soon as possible.]

Optional requirements for GPU accelerated computing

-- Please ensure that your installation(NVIDIA driver, CUDA) is appropriate for your hardware and check if the GPUs are available.
```bash
CUDA: 11.1
TensorFlow: 1.15.4
```
You can build the environment with Anaconda(or miniconda):
```bash
$ conda create -n ssc python==3.6.13
$ conda activate scc
(scc) $ pip install -r requirements.txt
```

## Usage - sample data

```bash
(scc) $ python scc.py --n 5 --input sample_data/ --output sample_data/out/ --order 1
```
Parameter
- [--n]: The number of CPU cores.
- [--input]: Directory for both the enrichment matrices and adjacency matrices.
- [--output]: Directory for output
- [--order]: The order of adjacency matrices

## Data Format - input
- Adjacency matrices for input should be

|  | AGTGTGGTCTATTGTG-1 | GCTATCGCGGCGCAAC-1 |  ... |Barcode n|
| :----:| :----: | :----: |  :----: | :----: | 
|AGTGTGGTCTATTGTG-1|0|1|...|0|  
|AGTGTGGTCTATTGTG-1|0|0|...|1|
|...|...|...|...|...|
|Barcode n|0|0|...|0|

- Enrichment matrices for input should be

|  | AGTGTGGTCTATTGTG-1 | GCTATCGCGGCGCAAC-1 |  ... |Barcode n|
| :----:| :----: | :----: |  :----: | :----: | 
|Type_A|0.01306|0.00010|...|0.00022|  
|Type_B|0.33542|0.48310|...|0.07694|
|...|...|...|...|...|
|Type_n|0.05631|0.06172|...|0.04630|

## Data Directory - input
- Each sample should have matching names for both the enrichment matrices and adjacency matrices.
- You can set the path for the input folder, but ensure that all enrichment matrices are located within the 'enrich' folder.
- The folder containing adjacency matrices should follow the format 'adj_1', 'adj_3', etc., where the latter number indicates the order.

```bash
sample_data/
	├── enrich
    		├── Sample_A.csv
    		├── Sample_B.csv
    		└── ...
	├── adj_1
    		├── Sample_A.csv
    		├── Sample_B.csv
    		└── ...
	├── adj_2
    		├── Sample_A.csv
    		├── Sample_B.csv
    		└── ...
	└── out
```
## Data Directory - output

- Each sample's folder contains distributed computation outputs, which will eventually be concatenated and saved in the designated output path.

```bash


 sample_data/out/
	├── Sample_A_1
    		├── Sample_A_1_0.csv
    		├── Sample_A_1_1.csv
     		├── Sample_A_1_2.csv
    		└── ...
	├── Sample_B_1
    		├── Sample_B_1_0.csv
    		├── Sample_B_1_1.csv
                ├── Sample_B_1_2.csv
    		└── ...
	├── Sample_C_1
    		├── Sample_C_1_0.csv
    		├── Sample_C_1_1.csv
                ├── Sample_C_1_2.csv
    		└── ...
	├── Sample_A_1.csv --> (final output for Sample_A)
        ├── Sample_B_1.csv --> (final output for Sample_B)
        └── Sample_C_1.csv --> (final output for Sample_C)
```
