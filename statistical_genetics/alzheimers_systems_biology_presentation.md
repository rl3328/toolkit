# Integrated Systems Approach to Late-Onset Alzheimer's Disease
### Zhang et al. · *Cell* 153(3): 707–720 · 2013

---

## Overview

| | |
|---|---|
| **Disease** | Late-onset Alzheimer's disease (LOAD) |
| **Approach** | Integrative gene network analysis |
| **Data** | 1,647 post-mortem brain samples, 549 subjects |
| **Transcripts profiled** | 39,579 per sample |
| **SNPs genotyped** | 838,958 per subject |
| **Brain regions** | Prefrontal cortex (PFC), Visual cortex (VC), Cerebellum (CB) |
| **Published** | *Cell*, April 25, 2013 |

---

## Background & Motivation

Late-onset Alzheimer's disease (LOAD) is the most common neurodegenerative disease, with incidence expected to **double by 2050**. Despite decades of research, no effective disease-modifying treatment exists.

**The problem with prior approaches:**
- Most LOAD research relied on mouse models of *early-onset* familial AD, based on three rare mutations in *APP*, *PSEN1*, and *PSEN2*
- Drug discovery focused almost exclusively on the **amyloid-β pathway** — with very high attrition rates
- GWAS identified risk loci in immunity (*CLU*, *CR1*, *CD33*), lipid processing (*APOE*, *ABCA7*), and endocytosis (*PICALM*, *BIN1*), but these explain only a small fraction of heritability

**The key question:** Can we clarify LOAD pathology by zooming out to the *pathway and network level*, rather than focusing on individual genes?

---

## Study Design — Three-Step Systems Approach

```
Step 1: Build co-expression networks
        ↓
   Gene expression profiling (microarray) + SNP genotyping
   → Identify top 1/3 most variable genes (~13,193 per region)
   → Compute pairwise co-expression (topological overlap)
   → Hierarchical clustering → gene modules

Step 2: Detect network remodeling (LOAD vs. healthy)
        ↓
   Compare LOAD and non-demented networks
   → Compute MDC (Modular Differential Connectivity) per module
   → Classify: Gain of Connectivity (GOC) vs. Loss of Connectivity (LOC)
   → Functional enrichment + correlation to neuropathology

Step 3: Identify causal regulators
        ↓
   Construct directed Bayesian networks using cis-eSNPs as anchors
   → Rank modules for LOAD relevance
   → Identify key causal regulators via combined scoring
   → Validate in mouse microglia overexpressing TYROBP
```

---

## Mathematical Methods

### 1. Co-expression Network Construction (WGCNA-based)

Gene co-expression networks were built using the top one-third (~13,193) most variably expressed genes per brain region. Pairwise **topological overlap** was computed as the co-expression similarity measure, then hierarchical clustering defined discrete gene modules.

This approach captures not just which genes are expressed, but which genes are **co-regulated together** — revealing the modular architecture of the transcriptome.

---

### 2. Modular Differential Connectivity (MDC)

> **Core novel metric of the paper.**

To detect how LOAD *reorganizes* the network, the authors defined MDC as the ratio of average gene-gene connectivity in LOAD vs. healthy brains, for any given module Ω:

$$\delta_\Omega(x, y) = \frac{\displaystyle\sum_{i=1}^{N-1} \sum_{j=i+1}^{N} k_{ij}^{x}}{\displaystyle\sum_{i=1}^{N-1} \sum_{j=i+1}^{N} k_{ij}^{y}}$$

Where:
- $k_{ij}$ = connectivity (co-expression strength) between genes $i$ and $j$
- $x$ = LOAD state
- $y$ = non-demented (healthy) state
- $N$ = number of genes in the module

**Interpretation:**

| MDC value | Meaning |
|---|---|
| MDC > 1 | Gain of Connectivity (GOC) — genes become more co-regulated in LOAD |
| MDC < 1 | Loss of Connectivity (LOC) — genes become less co-regulated in LOAD |
| MDC ≈ 1 | Conserved — network structure unchanged |

**Results:** At 10% FDR, **54% of modules showed GOC** and **4.5% showed LOC**. Crucially, only 6% of these changes overlapped with standard differential gene expression — confirming that MDC captures a *previously untapped* disease signal.

---

### 3. Bayesian Network Construction with cis-eSNP Priors

To convert undirected co-expression networks into **directed causal networks**, the authors used brain *cis* expression SNPs (eSNPs) as causal anchors.

**Why eSNPs?**
Cis-eSNPs are genetic variants in linkage disequilibrium with (or identical to) causal variants that affect the expression of a nearby gene. Because DNA sequence is upstream of gene expression, eSNPs serve as **natural perturbations** to infer directionality:

$$\text{DNA variant (eSNP)} \rightarrow \text{gene A expression} \rightarrow \text{gene B expression}$$

The study detected **11,318 unique cis eSNP–transcript pairs** across the three brain regions at 10% FDR — the largest brain eSNP dataset reported at the time.

These eSNPs were used as Bayesian priors to orient edges between nodes, producing a **directed probabilistic Bayesian network** for each co-expression module.

---

### 4. Causal Regulator Ranking Score

To identify the *most influential* genes within a module, a combined ranking score was defined that integrates both regulatory strength and differential expression:

$$G_j = \prod_i g_{ji}$$

Where the discriminant value $g_{ji}$ for gene $j$ in sample $i$ is:

$$g_{ji} = \frac{\max_i(r_{ji}) + 1 - r_{ji}}{\displaystyle\sum_j r_{ji}}$$

- $r_{ji}$ = expression rank of gene $j$ in sample $i$
- Genes with **high $G_j$** control many downstream nodes AND are strongly differentially expressed in LOAD

This score was applied to rank causal regulators of the top immune/microglia module. **TYROBP scored highest** across all brain regions.

---

### 5. Module–Trait Correlation via PCA

To link network modules to clinical LOAD outcomes, the authors used **principal component analysis (PCA)** to reduce each module to a single representative expression axis (the first PC), then tested its correlation to 25 LOAD neuropathology traits using multiple regression.

Significance was assessed with the **false discovery rate (FDR)** method to correct for the large number of module–trait tests.

The immune/microglia module correlated to the greatest number of LOAD traits, including:
- Frontal cortex atrophy (r = 0.27, FDR = 0.018)
- Neostriatum atrophy (r = 0.28, FDR = 3.3×10⁻⁹)
- Ventricular enlargement (r = 0.17, FDR = 0.031)

---

### 6. Functional Enrichment — Fisher's Exact Test (FET)

Over-representation of canonical pathways in modules was measured using **Fisher's Exact Test**, corrected for the number of modules and functional categories tested:

$$p = \frac{\binom{K}{k}\binom{N-K}{n-k}}{\binom{N}{n}}$$

Where:
- $N$ = total genes in background
- $K$ = genes in the pathway of interest
- $n$ = genes in the module
- $k$ = overlap between module and pathway

This revealed that the immune/microglia module had the **most statistically significant functional enrichment** of all modules.

---

### 7. Cross-species Validation — Distance-dependent Prediction

To validate the causal network, mouse microglia overexpressing full-length or truncated (dominant-negative) TYROBP were RNA-sequenced. Differentially expressed genes were mapped onto the ~8,000-node human Bayesian network.

The key statistical finding: a **strong negative correlation** between the fraction of confirmed differentially expressed targets and their network path distance from TYROBP:

$$r = -0.82, \quad p = 4 \times 10^{-7}$$

This means genes *closer* to TYROBP in the network were *more likely* to change expression in the mouse experiment — validating that causal network structure predicts real biological responses.

---

## Key Findings

### Finding 1 — Widespread network remodeling in LOAD

LOAD brains showed massive reconfiguration of gene co-expression patterns. Of 111 modules identified, **54% gained connectivity** and **4.5% lost connectivity**. Only 6% of these changes overlapped with standard differential expression analysis — network disruption is largely an independent signal.

**GOC modules** were enriched for: immune response, unfolded protein, vascular system, neurogenesis, glucose homeostasis, synaptic transmission

**LOC modules** were enriched for: nerve myelination, GABA metabolism, neurotrophin signaling, cell cycle

---

### Finding 2 — Immune/microglia module is the top LOAD network

The immune/microglia module ranked #1 among all modules based on:
- Strongest functional enrichment (microglia, pathogen phagocytosis)
- Highest correlation to LOAD neuropathology traits
- Strongest eSNP enrichment in PFC
- Significant gain of differential connectivity (MDC 1.49–100, FDR < 0.001)

The module was dissected into five immune pathway families:

| Family | Function |
|---|---|
| Complement | Classical immune activation |
| Fc receptors | Antibody-mediated phagocytosis |
| MHC | Antigen presentation |
| Cytokines/Chemokines | Inflammatory signaling |
| Toll-like receptors | Innate immune recognition |

---

### Finding 3 — TYROBP is the master causal regulator

Bayesian inference ranked **TYROBP** (DAP12) as the top causal regulator of the immune/microglia network:
- Upregulated in LOAD (1.18-fold, p = 0.028; replicated 1.17-fold, p = 5.1×10⁻⁵)
- More differentially expressed than classical microglia markers AIF1 and CD68
- Physically associates with TREM2 — independently identified as a major LOAD risk gene
- Overexpression in mouse microglia recapitulated human TYROBP network changes

---

## Summary & Significance

| Contribution | Impact |
|---|---|
| MDC metric | New way to detect disease-associated network rewiring, independent of differential expression |
| Immune/microglia module | Implicates microglial dysfunction — not just amyloid — as central to LOAD |
| TYROBP as hub regulator | Data-driven therapeutic target with downstream network leverage |
| TREM2–TYROBP convergence | Network approach independently converged on the same pathway as subsequent human genetics studies |
| Cross-species validation | Causal network structure predicts gene perturbation responses (r = −0.82) |
| Disease specificity | Immune module disruption absent in Huntington's disease — specific to LOAD |

---

## Conclusions

This study demonstrates that **systems-level network analysis** of large human brain datasets can uncover disease mechanisms invisible to single-gene or pathway approaches. The causal network framework:

1. Identified the **immune/microglia system** as the most disease-relevant molecular network in LOAD
2. Pinpointed **TYROBP** as its key causal regulator — a bottleneck gene controlling hundreds of downstream targets
3. Validated the **predictive power** of causal network structure across species
4. Offered a new **drug discovery framework**: target network hubs that restore the system to a healthy state, rather than single amyloid-pathway proteins

> "The accumulated data suggest that TYROBP may be a therapeutic target in prevention of neuronal damage in LOAD."
> — Zhang et al., 2013

---

*Zhang B, Gaiteri C, Bodea LG, et al. Integrated systems approach identifies genetic nodes and networks in late-onset Alzheimer's disease. Cell. 2013;153(3):707–720. doi:10.1016/j.cell.2013.03.030*
