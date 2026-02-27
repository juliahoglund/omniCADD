#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
:Author: Christian Gross
:Contact: c.gross@tudelft.nl
:Date: 06-08.2018

This files goes through the ancestral sequence 
and compares it with the reference genome of interest
and the known population variants to identfy derived variants.

The identification is based on five criteria:
Cases		:1		2		3		4		5
Ancestor	:C		C		C		C		C
Reference	:C		A>0.9	A		A>0.9	A
Alternative	:A>0.9	G		G>0.9	C		A

Case 1: The reference genome and ancestor are the same at position n, 
		but in the population data, there is a nearly fixed derived variant with >0.9
Case 2: The reference genome and ancestor are different at position n, 
		and the reference allele has a frequency >0.9
Case 3: The reference genome and ancestor are different at position n, 
		but there is a third allele in the population with a frequency >0.9
Case 4: The reference genome and ancestor are different, but the alternative allele
		in the population data is not the one in the reference genome, i.e. is most 
		likely ancestral. Here the population reference allele at >0.9 is what is desired
Case 5: The ancestral and reference allele are different at position n

Case X: The reference genome and ancestor are different, but there is no data on that position
		in the population file. If a true variant and true missingness, it can be assumed to be 100% fixed
		in the reference species / population, and hence be monomorphic in the reference.
		not yet implemented

:Edited by: Seyan Hu
:Date: 19-10-2022
:Extension and modification: Julia Höglund
:Date: 2023-08-08
:Edited by: Job van Schipstal
:Date: 27-10-2023
# ...existing code...
