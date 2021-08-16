Construction of mosaic reference standards             
==========================================
<br/>

![Fig1자산 3@3x](https://user-images.githubusercontent.com/77031715/129295890-578a2439-049a-4543-91a2-8b78a8cc7cfe.png)
<br/>
<br/>
Introduction
------------
A robust pipeline for constructing mosaic reference standards. The reference standards are consisting of XXXX mosaic SNVs and INDELs in a wide range of variant allele frequencies, from 0.1% to 56%. Negative controls of wild-type positions and germline vairants are accompanied with abundant positive controls. The reference stands were constructed by mixing genetic materials of pre-genotped, six normal cell lines 
<br/>
   
   
<br/>

## 1. Genotyping six cell lines

  #### A. Germline variant calling

    - Alignment and preprocessings
    - Strelka2 
    - DeepVariant 
    - CNVkit   
  
  #### B. Generating candiates 

    - Positive controls 
      + Mutually exclusive germline variants    
    
    - Negative controls   
      + Common wildtype positions (Set A)   
      + MRC5 germine variants (Set B)   
    
## 2. Generating Set B 

    - Down-sampling of MRC5 (39 times with random seeds) 
    - Extraction of reads embedding positive controls from Set A
    - Replacement of extracted reads to MRC5 down-sampled data

## 3. Finalizing the reference standards

  #### A. Pileup for post-filters

  #### B. Apply post-filters

    - Sequencing coverage
    - Variant coverage
    - High-quality alternative alleles

    
    
    
