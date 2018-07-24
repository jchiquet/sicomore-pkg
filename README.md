# SIComORe
Selection of Interaction effects in COmpressed Multiple Omics REpresentation

## Description 
From a set of input matrices and phenotype related to the
    same set of individual, sicomore is a two-step method which 
    1. find and select groups of correlated variables in each input matrix
    which are good predictors for the common phenotype; 2. find the most
    predictive interaction effects between the set of data by testing for
    interaction between the selected groups of each input matrix

## Vignettes

Two vignettes to help use the package are available in the directory /vignettes. The first vignette "singleData.pdf" describes how to use getHierLevel() function for finding and selecting relevant groups of variable while the second "twoDataInteraction.pdf" how to use sicomore() function to detect interactions between 2 datasets.

## Authors
Christophe Ambroise : christophe.ambroise@genopole.cnrs.fr
Julien Chiquet : julien.chiquet@inra.fr 
Florent Guinot : guinotflorent@gmail.com
Marie Szafranski : marie.szafranski@@ath.cnrs.fr

## Maintainer
Florent Guinot : guinotflorent@gmail.com
