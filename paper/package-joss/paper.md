---
title: 'strait: Bounds tests for reported summary statistics'
tags:
  - R
  - metascience
  - research integrity
  - error detection
  - trustworthiness assessment
  - summary statistics
authors:
  - name: Ian Hussey
    orcid: 0000-0001-8906-7559
    affiliation: 1
affiliations:
  - name: Institute of Psychology, University of Bern, Switzerland
    index: 1
date: 24 July 2026
bibliography: paper.bib
---

# Summary

When a measure is bounded, such as a 1-to-7 Likert item, a 0-to-63 depression inventory, or a percentage, the reported mean constrains the standard deviation (SD) that the underlying data could have produced. The two summary statistics are not independent: pushing observations toward the extremes to inflate the SD also moves the mean, so for a given mean, sample size and scale range only a limited interval of SDs is arithmetically attainable. A reported SD outside that interval cannot have arisen from real data, whatever the study.

`strait` computes those bounds in closed form and turns them into a consistency verdict for a reported (mean, SD, *n*) triple. Constraints are nested, and each one narrows the interval: the scale limits alone; the limits and *n*; the limits, *n* and the mean; additionally the data being integers; the observed minimum and maximum, when reported; and, for multi-item scales, a reported Cronbach's $\alpha$. The package handles the fact that published values are rounded, deferring the granularity tests GRIM and GRIMMER [@brown2017grim; @anaya2016grimmer] to `scrutiny` [@jung2024scrutiny] rather than reimplementing them, and reports which bound binds and which test failed. Companion functions apply the check across a data frame, transform results to a scale-free percent-of-maximum-possible metric so that studies using different instruments can be pooled, and visualise the feasible region of means and SDs.

# Statement of need

Problematic studies — those with serious concerns about the trustworthiness of their data or findings, whether from misconduct or from critical error — distort evidence synthesis and can lead to harmful recommendations. Carlisle estimated that 26% of randomised controlled trials submitted to one anaesthesia journal were problematic [@carlisle2021analysis], and risk-of-bias tools are not designed to detect problems of this kind. Because such studies are often not otherwise identifiable, checks that test whether reported numbers are internally possible have become a practical necessity.

This need is now formalised. INSPECT-SR (INveStigating ProblEmatic Clinical Trials in Systematic Reviews), developed through a Delphi process with 71 experts and consensus meetings and endorsed by Cochrane, includes as item 4.8 of its Results domain: *"Are the means and variances of integer data impossible?"* [@wilkinson2025inspect]. That item attracted among the highest ratings in the Delphi, with a median importance of 9 (IQR 8 to 9) and 94.9% consensus, and was retained in 13 of 14 consensus votes. `strait` implements exactly that check.

Existing software covers neighbouring ground but leaves a gap. `scrutiny` [@jung2024scrutiny] implements the granularity tests GRIM and GRIMMER, which ask whether a reported mean or SD is consistent with integer data at a given sample size, but not whether the mean and SD are jointly attainable given the scale limits. `rsprite2` [@wallrich2023rsprite2] and `unsum` [@jung2025unsum] reconstruct candidate datasets by search — SPRITE [@heathers2018sprite] and CLOSURE respectively — which is powerful and more complete, but the cost grows with sample size and scale width, and a null result establishes impossibility only after the search space has been exhausted. `strait` fills the gap with closed-form bounds: each check is an arithmetic evaluation rather than a search, so screening thousands of reported triples is immediate, and an out-of-bounds verdict is a proof of impossibility rather than the outcome of an enumeration. The two approaches compose naturally, with the bounds used as a fast screen and reconstruction reserved for the cases that survive it.

The bounds themselves have been derived repeatedly and independently for ninety years — in mathematics [@popoviciu1935; @bhatia2000better], survey sampling [@muilwijk1966note], statistics education [@petocz2005upper], constraint programming [@pesant2005spread], psychology [@mestdagh2018sidelining] and medicine [@nowbar2014discrepancies] — usually in partial form, with a sharp ceiling but no floor or the reverse. `strait` implements them as a single nested framework with both bounds sharp at every level, including results that are novel to the accompanying methodological article: a count-parity correction to the mean-conditional ceiling, a quasi-integer minimum that is defined at every reported mean and so frees the lower-bound check from requiring the mean to pass GRIM, and bounds conditional on a reported Cronbach's $\alpha$.

# Functionality

`sd_bounds()` returns the smallest and largest sample SD consistent with whichever constraints are supplied, together with the rule that binds each bound. `brimmer()` turns those bounds into a consistent/inconsistent verdict for a reported triple, unrounding the reported values, naming the failing tests, and adding percent-of-maximum-possible transforms of the mean and SD; `brimmer_map()` applies it row-wise across a data frame. `sd_bounds_sample()` constructs a dataset that attains a bound, and `sd_bounds_curve()`, `umbrella_data()`, `plot_sd_region()` and `plot_umbrella()` build and draw the feasible region — as a continuous band, or, for strictly integer data, as the lattice of attainable reported values. The single-purpose bound primitives are also exported. The package is validated against exhaustive enumeration and, for the tuples that pass every test, against the CLOSURE algorithm implemented in `unsum` [@jung2025unsum].

# Acknowledgements

I thank the INSPECT-SR project team, and Lukas Jung for `scrutiny` and `unsum` and for discussion of the underlying algorithms.

# References
