## Table of Contents
1. [General Info](#general-info)
2. [Technologies](#technologies)
3. [Implementation](#implementation)
4. [Results](#results)

### General Info
***
Implementation of Border Apolarity Algorithm of Conner/Harper/Landsberg using
SageMath. Batch computing used to compute the many candidates. Used
computational cluster with SLURM as workload manager. Results will be cleaned
up and more streamlined in future.
## Technologies
***
A list of technologies used within the project:
* [SageMath](https://www.sagemath.org/): Version <9.0 
Recommend Installing LinBox as well.
## Implementation
***
It is sufficient to run Structure3_15.sh to reobtain results proving the border rank of the structure tensor of sl3 is greater than 15. The file borderapolarity3.sage contains appropriate functions and structure3_15.sage contains the driver code. It is recommended to run cases in parallel, in which case, a slurm file is used to distribute computations across a cluster.
## Results
***
Results for sl3 rank 15-- complete
Results for sl3 rank 16-- in progress
Results for sl3 rank 16 with flag condition-- in progress
Results for so4 rank 9-- complete
```
Side information: Previous versions in 'Previous Versions' folder
## Acknowledgements
***
Attribute: Initial version of borderapolarity.sage to Austin Conner
