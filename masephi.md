# Marker Selection Strategies for Circulating Tumor DNA Guided by Phylogenetic Inference

**Authors:** Xuecong Fu¹, Zhicheng Luo¹, Yueqian Deng², William LaFramboise³, David Bartlett³, Russell Schwartz¹,²

¹Department of Biological Sciences, Carnegie Mellon University  
²Ray and Stephanie Lane Computational Biology Department, Carnegie Mellon University  
³Allegheny Health Network Cancer Institute

**Corresponding author:** russells@andrew.cmu.edu

## Abstract

**Motivation:** Blood-based profiling of tumor DNA ("liquid biopsy") offers prospects for noninvasive early cancer diagnosis, treatment monitoring, and clinical guidance, but requires computational advances to become a robust quantitative assay of tumor clonal evolution. We propose new methods to characterize tumor clonal dynamics from circulating tumor DNA (ctDNA) through two specific questions: 1) How to apply longitudinal ctDNA data to refine phylogeny models of clonal evolution, and 2) how to quantify changes in clonal frequencies that may indicate treatment response or tumor progression.

**Results:** We estimate a distribution over plausible clonal lineage models using bootstrap samples over pre-treatment tissue-based sequence data. We then refine these lineage models and implied clonal frequencies over successive longitudinal samples. The resulting framework poses optimization problems to select ctDNA markers that maximize utility for reducing phylogeny uncertainty or quantifying clonal frequencies. Testing on synthetic data showed effectiveness at refining distributions of tree models and clonal frequencies to minimize tree distance measures relative to ground truth. Application to real tumor data demonstrated effectiveness in refining clonal lineage models and assessing clonal frequencies.

**Availability:** https://github.com/CMUSchwartzLab/Mase-phi.git

## 1. Introduction

### 1.1 Background on Liquid Biopsy

Circulating free DNA (cfDNA) in human blood, particularly tumor-derived cfDNA at elevated levels compared to healthy cell DNA, established the potential for liquid biopsy - blood-based profiling of solid tumor genomics. This elevation occurs due to:
- Enhanced release of tumor cell DNA
- Abnormal clearance of DNA debris from cell death  
- Presence of circulating tumor cells in blood

Liquid biopsy offers possibilities for improving cancer diagnosis and treatment, including:
- Early prognosis
- Detecting residual disease and relapse
- Treatment monitoring across various cancer types

### 1.2 Technical Limitations and Current Approaches

The technology faces technical limitations primarily due to:
- Challenge of separating tumor signals from much larger numbers of healthy cells
- Need for highly sensitive genomic assays with low signal-to-noise ratio

Current molecular testing methods include:
- **Deep sequencing:** Costly and time-consuming for repeated use
- **Multiplex-PCR:** Limited to relatively few pre-selected markers
- **Droplet digital PCR (ddPCR):** Highly sensitive quantitation of somatic variants at low levels
- **Targeted sequencing:** Strategies for enriching tumor DNA

### 1.3 Current Clinical Applications vs. Research Goals

Current clinical applications focus primarily on:
- Prognosis or recurrence detection
- Rather than precise quantitative analysis of tumor genetics

**Research Gap:** While there is rich literature on tumor evolutionary trajectories from genomic assays, the need to work with low precision or limited blood-based marker sets makes it infeasible to incorporate longitudinal blood samples into current tumor phylogenetics methods straightforwardly.

### 1.4 Study Objectives

This work addresses two previously unaddressed challenges:

1. **Phylogeny Refinement:** How to leverage liquid biopsy data using small marker sets to refine phylogenetic models from primary tumor to:
   - Correct errors
   - Reduce inference uncertainty
   - Expand trees to accommodate variants/clones not seen in earlier samples

2. **Optimal Marker Selection:** For scenarios requiring selection of small marker subsets for high sensitivity profiling, which markers are most informative for:
   - Optimally refining tumor phylogeny models
   - Characterizing changes in clonal frequency or tumor heterogeneity over time

## 2. Methods

### 2.1 Overall Workflow

The methodology addresses marker selection for liquid biopsy assuming selection of small marker sets for high-precision assays (e.g., ddPCR). The goal is developing personalized assays allowing rapid longitudinal corrections on a patient-specific basis.

**Key Assumption:** The problem is solved for the general case of selecting from any observed marker for each specific subject, though conceptually identical when limited to choosing from predefined marker sets with available PCR probes.

### 2.2 Problem Formulation Framework

#### 2.2.1 General Probabilistic Framework

Given a candidate (bootstrapped) tree set **T** = {T^k}_{k=1,...,K}, we define:

- **E^k:** Clonal tree structure matrix where E^k_{i,j} = 1 if clone i is parent of node j, otherwise 0
- **M^k:** Mutation assignment matrix where M^k_{i,l} = 1 if mutation g_l belongs to node i, otherwise 0  
- **F^k:** Clonal frequency array F^k = (F^k_i)_{i=1,...,N} encoding variant allele frequency (VAF) of each mutation in each clone
- **G:** Set of mutations {g_l}_{l=1...N}
- **Q_n:** Selected n ≤ N gene markers {g_{f^n(j)}}_{j=1...n} via mapping f^n

**Objective:** Minimize E_{T^q∈T} E_{T^k∈T\{T^q}} log E_{(R^{f^n},S^{f^n})∼T^k} P(R^{f^n}, S^{f^n}|T^q)

Where:
- **R^{f^n}:** Random variables describing ddPCR read counts for probe set
- **S^{f^n}:** Binary matrix mapping biomarkers to clones in tree structure

#### 2.2.2 Likelihood Components

The likelihood splits into:
1. **Read count probability:** P(R^{f^n}|F^q, M^q)  
2. **Tree structure probability:** P(S^{f^n}|E^q)

**Read Count Modeling:**
Assuming high ddPCR sequencing depth, read counts follow normal distribution instead of binomial for computational convenience:

- r_j ∼ N(μ^k_j, (σ^k_j)²)
- μ^k_j = D·F^k_{M^k(f(j))}  
- (σ^k_j)² = D·F^k_{M^k(f(j))}(1 - F^k_{M^k(f(j))})

Where D is read depth.

**Structural Perturbation Modeling:**
Uses ancestor-descendant distance to model distance between subtrees:
- Define subtree S^k_{Q_n} of T^k given marker genes Q_n using ancestor-descendant matrix A
- For T^k: A^k_{i,j} = 1 if mutation i is ancestor of j
- Probability for each unit change of AD distance is λ

### 2.3 Selecting Markers to Minimize Phylogenetic Uncertainty

#### 2.3.1 Tree Distribution Estimation

Due to low signal-to-noise ratio of liquid biopsy samples, we expect poor results with standard phylogenetic inference tools. Therefore:

- Pose phylogenetic inference in terms of **distributions of trees** rather than single optimum
- Estimate distribution through **bootstrapping over sequence reads**
- Use existing tumor phylogeny methods designed to return single optimal tree

**Alternative Considered:** Bayesian phylogeny methods might provide more principled alternative, though designing efficient Bayesian sampler for non-trivial tumor phylogeny models is challenging.

#### 2.3.2 Objective Function

After combining read count and structural components:

```
min_{f^n} ∑^K_{q=1} ∑^K_{k=1,k≠q} [
  -1/2 ∑^n_{j=1} log(2π((σ^q_j)² + (σ^k_j)²))
  -1/2 ∑^n_{j=1} (μ^k_j - μ^q_j)²/((σ^q_j)² + (σ^k_j)²)
  + ∑^n_{i=0} ∑^n_{j=0} |Â^k_{i,j}(f) - Â^q_{i,j}(f)| log λ
]
```

This objective function evaluates utility of given marker set for reducing uncertainty measure in expectation over tree density and sequence read sampling.

### 2.4 Refining Tumor Phylogenetic Tree Distribution

#### 2.4.1 Statistical Hypothesis Testing Method

**Approach:** Accept or reject candidate trees based on whether measured marker frequencies are plausibly consistent with given tree topology.

**Assumptions:**
- Clone frequencies may change over time revealing some trees as implausible
- Set of clones and tree topology unchanged during follow-up

**Wald Test Implementation:**
For marker pair, test whether read counts at specific time point are consistent with tree topology using Wald test statistic:

```
W = δ̂/ŝê = (f̂_{g2} - f̂_r_{g1})/√(f̂_{g1}(1-f̂_{g1})/depth₁ + f̂_{g2}(1-f̂_{g2})/depth₂)
```

- Reject hypothesis when W > z_α with α = 0.05
- Apply Bonferroni correction for multiple marker pairs
- Remove tree structures rejected by any pairwise marker test

#### 2.4.2 Bayesian Update Method

**Motivation:** Statistical hypothesis testing may be too strict, allowing only accept/reject decisions. ctDNA data might provide evidence for/against trees without definitive acceptance/rejection.

**Bayesian Framework:**
```
P(S|R¹,R⁰) ∝ P(S,R¹|R⁰) = P(R¹|S,R⁰)P(S|R⁰) = P(R¹|S)P(S|R⁰) = ∫_f g(R¹,f|S)P(S|R⁰)
```

Where:
- **P(S|R¹,R⁰):** Updated weights
- **P(S|R⁰):** Original weights (bootstrap counts)
- **f:** Possible VAFs for mutation markers from liquid biopsy
- **g(R¹,f|S):** Probability density of clonal VAFs f and corresponding read counts given structure S

**Two-Marker Case:**
Four possible relationships between markers in tree structure:
1. **Case (a):** Marker 1 ancestor of marker 2 → f₁ > f₂
2. **Case (b):** Marker 2 ancestor of marker 1 → f₁ < f₂  
3. **Case (c):** Same clone → f₁ = f₂ = f
4. **Case (d):** Different branches → f₁ + f₂ < 1

### 2.5 Selecting Markers to Track Subclonal Population Changes

#### 2.5.1 Deterministic Trees Approach

**Problem Statement:** Given most likely tree index k̂, select n gene markers to maximize sum of weighted tracked clones in T.

**Weight Definition:** Set weights as estimated clonal fractions from previous time point, maximizing estimated fraction of tumor tracked.

**Mathematical Formulation:**
- Create clonal frequency array F̄ = (F̄^k_i)_{i=1,...,Clone_num} from most recent VAF estimates
- F̄^k_i = F^k_i - ∑_j F^k_j where j ∈ Children(i)
- Binary output vector z identifies chosen markers
- x = M^k̂z where x_j = 1 indicates node j tracked

**Information Assumptions:**
1. **Partial Information:** Can only measure markers without accounting for inheritance
2. **Complete Information:** Must measure marker of clone and its children's clones to identify actual clonal fraction

**Complete Information Optimization:**
```
max_z t^T F̄^k̂
subject to: y^k̂ = σ(Ê^k̂ ⊙ x^k̂ - Ê^k̂)
where t = abs(1 - ȳ^k̂)
```

#### 2.5.2 Candidate Tree Set Extension

**Objective for Tree Distribution:**
```
max_z ∑_{k̄∈S} (x^k̄)^T F̄^k̄  (partial information)
max_z ∑_{k̄} (t^k̄)^T F̄^k̄  (complete information)
```

Where S ⊂ {1,...,K} represents subset of candidate trees for optimization.

### 2.6 Tracking Tumor Subclonal Population

**Method:** Infer clone frequency using mean VAF of markers in given clone minus sum of mean VAFs of markers inferred for child clones.

**Key Assumption:** Only choosing markers in copy number neutral regions, treating VAF as proxy for cancer cell fraction (CCF).

### 2.7 Implementation Details

**Optimization:** All probe optimization problems formulated as integer linear programs (ILPs)
- Binary variable z ∈ {0,1}^{N×1} where z_i = 1 if i-th mutation selected
- Constraint: ∑^N_{i=1} z_i = n for marker number constraint
- Solver: Python with Gurobi

**Phylogeny Methods Tested:**
1. **PhyloWGS:** MCMC-based method
2. **deconv:** Simplified ILP method (TUSV-ext restricted to SNVs)

## 3. Simulations

### 3.1 Simulation Framework

**Tree Structure Generation:**
- Parameterized by number of subclones and maximum degree of each subclone
- Randomized subclone distributions for total tree structure
- Beta distribution to generate true allele frequencies for each subclone

**Frequency Modeling:**
- Dirichlet probability distribution to randomize clonal frequencies with tumor tissues
- Different frequencies at different tumor sites or tissue vs. blood

**Masking Strategy:**
- Only subclones nearer root observed in tissue samples
- Remaining clones only observable in liquid biopsy samples
- Normalize fractions of observed clones so fractions sum to one

**Liquid Biopsy Simulation:**
- Add 0.9 frequency to normal cell (dilution effect)
- Subsequent Dirichlet random variable for remaining subclone frequency
- Poisson random variable for mutation number assignment

**Read Generation:**
- Assign read depth for tissue and blood samples
- Poisson random variable for total read count per variant per sample
- Binomial distribution with VAF probability for variant/reference reads

**ddPCR Simulation:**
- Assign droplet collection number
- Determine mutant-positive droplets based on normalized real data
- Use mean droplets and known VAF to simulate mutant positive detection

### 3.2 Phylogeny Method Comparison

**Methods Compared:**
- **PhyloWGS:** MCMC-based tumor phylogenetics
- **deconv:** SNV-only version of TUSV-ext

**Distance Metrics:**
- **CASet:** Tumor phylogeny-specialized distance measure
- **DISC:** Alternative tumor phylogeny distance measure

**Key Finding:** PhyloWGS yields good results even with single tissue sample, with accuracy improving with more tissue samples. deconv performs poorly by CASet measure but comparably to PhyloWGS by DISC measure.

## 4. Results

### 4.1 Simulated Data Results

#### 4.1.1 Phylogeny Inference Accuracy

**Setup:** Comparison of PhyloWGS and deconv using paired liquid biopsy and multi-regional tissue samples.

**Results:** 
- PhyloWGS consistently outperforms deconv, especially with CASet distance metric
- Performance improves with increasing number of tissue samples
- Substantial room for improvement in tree inference, supporting value of refinement using subsequent ctDNA samples

#### 4.1.2 Full Pipeline Testing

**Bootstrap Approach:**
- Generated 100 bootstrapped samples per simulation case
- Inferred tree for each bootstrap replicate using PhyloWGS
- Applied marker selection methods on bootstrap tree sets

**Tree Refinement Results:**
- Using limited marker numbers can shift tree distribution toward ground truth
- Almost all cases yield refined trees with reduced weighted distance to ground truth
- Weights of best trees relative to alternatives improved

**Marker Selection Comparison:**
1. **Tree uncertainty reduction markers:** Substantially better at refining tree structure
2. **Clonal frequency estimation markers:** Poor performance at tree refinement (as expected)
3. **Random markers:** Poor performance at tree topology updating

#### 4.1.3 Clonal Fraction Monitoring

**Results:** 
- Increasing ability to characterize clonal frequencies with larger numbers of optimally chosen markers
- Random marker sets typically unstable, requiring more markers for minimal accuracy
- Demonstrated effectiveness across 10 simulation cases

### 4.2 Real Data Application: TRACERx Lung Cancer Case

#### 4.2.1 Dataset Description

**Case:** CRUK0044 from TRACERx study
- Three primary tumor multi-regional samples sequenced
- Six consecutive temporal liquid biopsy samples
- Focus on tree structure-based marker selection (based on simulation results)

#### 4.2.2 Tree Refinement Results

**Methodology:**
- Iterative marker selection for each time point
- Selected subset of larger marker set profiled by TRACERx study
- Used selected markers to update empirical tree distributions

**Key Findings:**
1. **Reinforcement:** Blood samples reinforced tree structure originally most frequently observed in bootstrap replicates
2. **Confidence Increase:** ctDNA allowed increased confidence in inference by providing evidence against other possible topologies
3. **Post-surgery dynamics:** First sample after surgery had least ctDNA remaining, making tree distribution adjustment less stable
4. **Stabilization:** Tree distribution stabilized after first few blood draws with enhanced confidence

#### 4.2.3 Clonal Population Tracking

**Marker Selection Strategy:**
- Used most probable tree structure as presumed correct tree
- Selected top five markers occurring in at least two of three multi-regional sample sets

**Temporal Dynamics Observed:**
1. **Overall trend:** Gradual then accelerating expansion of several subclones after surgery
2. **Early detection:** Trend apparent by day 91, though some dominant clones present at negligible levels until day 174
3. **Clonal expansion:** Sharp rise of initially rare subclone 5, becoming dominant by day 259 post-relapse

**Technical Issues:**
- Inference of negative clonal frequency for once-common clone 1
- Negative frequency indicates markers for descendants observed at higher frequency than clone markers
- Impossible result if tree and marker frequencies exact, but occurs due to imprecision in marker frequencies or tree topology error

## 5. Discussion

### 5.1 Contributions and Significance

**Primary Contributions:**
1. **Computational framework** for interpreting ctDNA data for refining phylogenetic tree models and tracing clonal frequencies
2. **Optimization methods** for selecting limited high-sensitivity markers for liquid biopsy in cancer progression tracking
3. **Proof of concept** demonstrating potential to bring liquid biopsy more effectively to clonal lineage studies and clinical applications

**Performance Validation:**
- **Simulated data:** Methods yield good accuracy at both tasks with limited markers, substantially improved over random selection or markers selected for different tasks
- **Real data:** Demonstrated potential for precise, quantitative tracking of clonal dynamics changes over time

### 5.2 Limitations and Future Directions

#### 5.2.1 Current Limitations

**Modeling Framework:**
- Largely proof of concept requiring extension in multiple ways
- Could adapt to more sophisticated Bayesian models of tree space
- Limited to two-marker Bayesian updates (more rigorous k-marker generalization left for future work)

**Biological Scope:**
- Present application focused on SNVs only
- Need to consider structural variations (SVs) - likely high impact and effective probes
- Copy number alterations (CNAs) often mechanism of tumor driver genes, can confound SNV VAF interpretation

#### 5.2.2 Future Extensions

**Technical Improvements:**
1. **Bayesian modeling:** More principled but tractable strategies for larger marker sets
2. **Objective functions:** Alternative functions more specifically tuned to clinical questions
3. **Multi-variant types:** Integration of SNVs, SVs, and CNAs

**Clinical Translation:**
- Considerable work needed to leverage new informatics capabilities for public health interventions
- Clinical treatment improvements for cancer care
- Validation in larger patient cohorts and different cancer types

### 5.3 Broader Impact

**Computational Methods Value:**
- Demonstrates power of computational methods to improve:
  - Marker selection
  - Clonal lineage reconstruction  
  - Clonal dynamics profiling
- Enables more precise and quantitative assays of tumor progression

**Clinical Potential:**
- Framework provides foundation for personalized liquid biopsy assays
- Enables quantitative monitoring of treatment response and tumor progression
- Supports clinical decision-making through precise clonal dynamics characterization

## 6. Implementation and Availability

**Software:** Available at https://github.com/CMUSchwartzLab/Mase-phi.git

**Technical Requirements:**
- Python with Gurobi solver for ILP optimization
- Compatible with PhyloWGS and TUSV-ext phylogeny methods
- Supports ddPCR and other high-sensitivity marker assay platforms

## 7. Acknowledgments and Funding

**Funding Sources:**
- Gift from Highmark Healthcare
- National Human Genome Research Institute of the National Institutes of Health (R01HG010589)

**Acknowledgments:** Oana Carja and Thomas Rachman for helpful discussions

---

**Contact Information:**
- **Corresponding Author:** Russell Schwartz (russells@andrew.cmu.edu)
- **Repository:** https://github.com/CMUSchwartzLab/Mase-phi.git