# The circadian clock time tunes axonal regeneration 

> Published in [Cell Metabolism](https://doi.org/10.1016/j.cmet.2023.10.012)
> November 2023

Di Virgiliis, F., **Mueller, F.** et al. (2023) The circadian clock time tunes axonal regeneration. *Cell Metabolism*  35(12) pp.2153-2164.

## Abstract
Nerve injuries cause permanent neurological disability due to limited axonal regeneration. Injury-dependent and -independent mechanisms have provided important insight into neuronal regeneration, however, common denominators underpinning regeneration remain elusive. A comparative analysis of transcriptomic datasets associated with neuronal regenerative ability revealed circadian rhythms as the most significantly enriched pathway. Subsequently, we demonstrated that sensory neurons possess an endogenous clock and that their regenerative ability displays diurnal oscillations in a murine model of sciatic nerve injury. Consistently, transcriptomic analysis showed a time-of-day-dependent enrichment for processes associated with axonal regeneration and the circadian clock. Conditional deletion experiments demonstrated that Bmal1 is required for neuronal intrinsic circadian regeneration and target re-innervation. Lastly, lithium enhanced nerve regeneration in wild-type but not in clock-deficient mice. Together, these findings demonstrate that the molecular clock fine-tunes the regenerative ability of sensory neurons and propose compounds affecting clock pathways as a novel approach to nerve repair.

## RNA-sequencing analysis 

### **Circadian dataset**

    Methodology: Peripheral (sciatic nerve) injury at *Zeitgeber* (ZT) 8 and ZT20 in C57bl/6 mice (ages 6-8 weeks); 72 hr post-injury, L4-L6 dorsal root ganglia were extracted. A neuronal enriched cell population was selected and processed for RNA-sequencing using Illumina NovaSeq 6000 platform.

    RNA-sequencing: 
    - Convert using bcl2fastq conversion software
    - Quality control: FastQC-v.0.11.9
    - Trimmed for low-quality reads and adapters: TrimGalore!-v0.6.6
    - Mapping and alignment: Salmon-v1.6.0 in mapping-based mode using M27 assembly
    - Reformat for differential expression analysis: Tximeta-v1.12.4
    - Differential expression analysis: DESeq2-v1.34.0 using R-v4.1.2

    Deposited in NCBI GEO under the accession code: [GSE235687](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235687).

### **Previously published dataset**

    Datasets include: 
    - [GSE97090](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97090): Sham, sciatic nerve injury (SNA), laminectomy (LAM), dorsal column injury (DCA).
    - [GSE59547](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi): Sham, sciatic nerve crush (SNC).
    - [GSE125793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi): environmental enrichment (EE), EE+SNA
    -  Available upon request: intermittent fasting (IF) vs *ad libitum* (AL).


## Figures and associated scripts
- Figure 1: Circadian rhythms are common regulators of axonal regeneration.

- Figure 2: Transcriptional anlaysis in DRG after injury identifies time-of-day-dependent processes associated with axonal regeneration.

- Figure 5: Bmal1-dependent H3k27ac after injury and long-lasting time-of-day-dependent regeneration of DRG neurons.

- Suppl. Figure 2: *related to Figure 2 and 3*. PCA plots and odd ratio correlational analysis of the time-of-day transcriptomic analysis of DRG neurons after injury. 

- Table S2: Extended table reporting gene expression data (gene name, fold change, p-value) of
148 DE Clock associated genes (CAGs) and Regeneration associated genes (RAGs) found in the
149 RNA-seq of DRG after an injury performed at ZT20 vs ZT8 (see also Figure 2).

- Table S3: List of RAGs predicted to contain the Bmal1 motif obtained from in silico
153 transcription factor binding site (TFBS) analysis.


