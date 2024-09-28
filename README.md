# BOA Estimation Pipeline

## Introduction to BOA
**Breed of Origin of Allele (BOA)** identifies the breed from which a specific genetic marker originates. This is crucial for crossbred or composite livestock populations, as it provides insights into how different breeds contribute to different traits.

BOA helps researchers and breeders to:
- Pinpoint genetic contributions from different breeds to key traits.
- Analyze breed-specific marker effects for more accurate selection.
- Map breed composition across the genome.
- Identify signatures of selection, revealing regions under genetic pressure.

By understanding breed contributions at the allele level, BOA supports better-informed breeding decisions and helps optimize genetic improvement strategies.

## Pipeline Overview

### LAMP-LD: Estimating Breed of Origin

The BOA pipeline relies on **LAMP-LD** (Local Ancestry in Admixed Populations using Linkage Disequilibrium) to estimate breed of origin. LAMP-LD analyzes patterns of linkage disequilibrium in a reference set of animals too assign ancestry in a admixed set of animals. 

**Key parameters**:
- **Window Size**: Defines the size of the genomic regions (number of genetic markers) for ancestry analysis. A balance between resolution and computational efficiency.
- **Hidden States**: Correspond to the local ancestries on each chromosome within each window. More states increase complexity but can improve accuracy.

### Tuning BOA Pipeline

The **Tuning Pipeline** optimizes these parameters—window size and hidden states—to maximize the accuracy of breed origin assignments. This ensures the BOA estimation is fine-tuned to your dataset's specific characteristics.

**Comment from the creator**: Optimized does not always mean best. The best parameters should be optimal or near optimal, and they are dictated by your population, marker array, and project objectives. The tuning performed here is done using purebred and F1 individuals, which can introduce a bias toward larger window sizes. If your population has more historical recombination events, you may benefit from using a smaller window size that yields comparable results to the optimized window size.

Be sure to align your tuning with the specific structure of your population and the objectives of your project.


### Breed of Origin Pipeline

After deciding the best window size and number of hidden states, the **Breed of Origin Pipeline** performs the main BOA estimation. It assigns breed origin to alleles across the genome, producing detailed outputs that can be used for further research, such as identifying breed-specific marker effects or detecting regions under selection.

## BOA Software Used

The BOA pipelines are designed for **Linux-based HPC systems** with **SLURM** for job management. Core tools include:
- **LAMP-LD**: Found in the `softwares` directory, for estimating breed of origin.
- **PLINK v2 & v1.9**: Used for filtering and managing genomic data prior to BOA estimation.
- **R**: Utilized for additional processing and analysis of the results.

### System Requirements:
- Linux OS
- SLURM for HPC job scheduling
- Sufficient memory and computational resources to handle large genomic datasets

## Citations

If you use this repository in your work, please cite the following papers:

- **LAMP-LD**:  
  - Pasaniuc et al. (2009). _Inference of locus-specific ancestry in closely related populations_. [DOI: 10.1093/bioinformatics/btp197](https://doi.org/10.1093/bioinformatics/btp197)
  - Baran et al. (2012). _Fast and accurate inference of local ancestry in Latino populations_. [DOI: 10.1093/bioinformatics/bts144](https://doi.org/10.1093/bioinformatics/bts144)

- **BOA Estimation Pipeline**:  
  - Zayas, G., et al. (2024). _Breed of origin analysis in genome-wide association studies: enhancing SNP-based insights into production traits in a commercial Brangus population_. [DOI: 10.1186/s12864-024-10465-1](https://doi.org/10.1186/s12864-024-10465-1)

## How to Use the Repository

### Setting Up the Environment

To get started, clone the repository and run the setup script:
1. Clone the repository:
   ```bash
   git clone https://github.com/gzayasPR/BOA_Estimation.git
   cd BOA_Estimation

   ```
2. **Run the setup script**:
The `Setup.sh` script will handle the installation of the required software and tools, including LAMP-LD, PLINK, and R packages.
   ```bash
   bash Setup.sh   
   ```
### Running the Pipelines

Once the environment is set up, you can access the two pipelines located within the `code` directory:

- **Tuning BOA Pipeline**: This pipeline, located in `code/Tuning_BOA`, is run using the `Run_BOA_Tunning.sh` script to optimize the window size and number of hidden states for LAMP-LD.
- **Breed of Origin (BOA) Pipeline**: This pipeline, located in `code/Breed_of_Origin`, is run using the `Run_BOA.pipeline.sh` script to perform the actual breed of origin assignment across the genome using the tuned parameters.

For more detailed information about the pipelines, refer to the documentation found in the `doc` directory.

## Contributions & Issues

We welcome contributions from the community! If you'd like to contribute to this project:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Submit a pull request with a detailed description of the changes.

If you encounter any issues or have suggestions for improvements, feel free to submit an issue via the GitHub [Issues tab](https://github.com/gzayasPR/BOA_Estimation/issues).

For any inquiries, feel free to contact Gabriel Zayas at **gzayas97@ufl.edu** or **gazspr@gmail.com**.
