# Bioinformatics analysis of E3 ubiquitin ligase family
### Multi-template Homology Auto modeling 

Materials for semester project in Bioinformatics Institute (spring 2021)

__Students__:   Pyankov I., Shemyakina A.

__Supervisor__: Popov P. (Skoltech)

**Table of Contents**

[TOC]

## Abstract

This work is inspired by the development process of Proteolysis-targeting chimeras (PROTACs) and related molecules that induce targeted protein degradation by the ubiquitin-proteasome system. They represent a new therapeutic modality and are the focus of great interest, however the progress is hindered by the low efficiency of protein crystallography that provides E3 ubiquitin ligases' 3D structures required in the initial steps of PROTAC development. The automated _in silico_ modeling tool could assist in expanding the number of enzymes available for development of targeted protein degradation systems.

## Objective

Develop script for performing Multi-template Homology Auto Modeling of a target protein.

## Plan

1) Develop the script.
2) Test script performance with E3 ubiquitin ligase target.

## Results

You can find the example result of script performance (with target protein human E3 ubiquitin-protein ligase TRIM69, EC 2.3.2.27) in the __results__ folder.

## Methods

Homologues search - mafft-homologs L-INS-i [[1]](#1);
modeling - rosetta_scripts Application [[2]](#2), FastRelax Mover [[3]](#3), energy landscape [[4]](#4);
quality assessment - Ornate [[5]](#5).

## Requirements

Compatibility is guaranteed for followed python packages versions:

	numpy==1.20.3
	pandas==1.2.3
	lxml==4.5.0
	requests==2.25.1
	urllib.request==3.8
	Bio==1.78
	matplotlib==3.3.4

Mafft version: v7.453
Rosetta version: rosetta.source.release-275 r275 2021.07+release.c48be26
Ornate has no version, but requiers tensorflow==1.14.0 or sooner (thus, python 3.5-3.7)

## User Guide

#### Description

The script consists of three parts.

First part downloads SWISS-MODEL Repository and/or World Wide Protein Data Bank, performs homologues search in the databases and downloads pdb structures of potential templates.

Second part processes the templates: chooses the correct chain from pdb file, selects templates above the identity percent threshold, calculates target coverage by remaining templates. Then the script follows the [RosettaCM tutorial](https://www.rosettacommons.org/demos/latest/tutorials/rosetta_cm/rosetta_cm_tutorial) and generates model(s).

Third part evaluates the resulting model(s) and generates the csv file with scores mean and sd and picture with score distribution across target length for each model.

#### Launch

In the beginning script requests in the command line the required information. During execution script informs user about target coverage by selected templates and waits for reply whether to proceed or not.

First part generates its output files in provided working directory. The results of the second part are contained in the __Modeling__ folder. The results of the third part are contained in the __score__ folder (Ornate output) and  __score_result__ folder (summary of Ornate output) .

## References
<a id="1">[1]</a> 
Katoh K. et al. MAFFT version 5: improvement in accuracy of multiple sequence alignment //Nucleic acids research. – 2005. – Т. 33. – №. 2. – С. 511-518. doi: 10.1093/nar/gki198

<a id="2">[2]</a> 
Fleishman SJ, Leaver-Fay A, Corn JE, Strauch E-M, Khare SD, Koga N, Ashworth J, Murphy P, Richter F, Lemmon G, Meiler J, and Baker D.  (2011).  RosettaScripts: A Scripting Language Interface to the Rosetta Macromolecular Modeling Suite.  PLoS ONE 6(6):e20161.  doi: 10.1371/journal.pone.0020161.

<a id="3">[3]</a>
Khatib F, Cooper S, Tyka MD, Xu K, Makedon I, Popovic Z, Baker D, and Players F.  (2011).  Algorithm discovery by protein folding game players.  Proc Natl Acad Sci USA 108(47):18949-53.  doi: 10.1073/pnas.1115898108.

<a id="4">[4]</a> 
Maguire JB, Haddox HK, Strickland D, Halabiya SF, Coventry B, Griffin JR, Pulavarti SVSRK, Cummins M, Thieker DF, Klavins E, Szyperski T, DiMaio F, Baker D, and Kuhlman B.  (2020).  Perturbing the energy landscape for improved packing during computational protein design..  Proteins "in press".  doi: 10.1002/prot.26030.

<a id="5">[5]</a> 
Pages G., Charmettant B., Grudinin S. Protein model quality assessment using 3D oriented convolutional neural networks. bioRxiv. – 2018. doi: 10.1093/bioinformatics/btz122
