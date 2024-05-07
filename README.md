# T-cell_and_glial_pathology_in_PD
Code and data for The spatial landscape of glial pathology and adaptive immune response in Parkinson's Disease 

# Spatial cross-correlation
This is a Python module for spatial cross-correlation analyses

# Installation
```bash
$ git clone https://github.com/dalhoomist/T-cell_and_glial_pathology_in_PD.git
```

# Set up Environment
Requirement
```bash
Python: 3.6.13
Numpy: 1.18.5
Pandas: 1.1.5
```
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

# Usage(sample data)

```bash
(scc) $ python scc.py --n 5 --input sample_data/ --output sample_data/out/ --order 1
```
Parameter
- [--n]: The number of multi-processing
- [--input]: Directory for input matrices
- [--output]: Directory for output
- [--order]: The order of adjacency matrices

# Data Directory
sample_data/
	├── enrich
    		├── Sample_A.csv
    		├── Sample_B.csv
    		└── ...
	├── adj_1
    		├── Sample_A.csv
    		├── Sample_B.csv
    		└── ...
	├── adj2
    		├── Sample_A.csv
    		├── Sample_B.csv
    		└── ...
	└── out

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
	├── Sample_A_1.csv
  ├── Sample_B_1.csv
  └── Sample_C_1.csv
