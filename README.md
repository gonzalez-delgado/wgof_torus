# wgof_torus
Goodness-of-fit testing approaches on the two-dimensional flat torus based on Wasserstein distance.

The goal of wgof_torus is to provide some practical approaches to perform two-sample goodness-of-fit tests for measures on the two-dimensional flat torus based on Wasserstein distance. These techniques are introduced and fully detailed in [arXiv] and consist on:

1. Testing the equality of the one-dimensional marginals on the circle.
2. Upper bounding p-values.

The file test_functions.R includes the required functions to implement both approaches, whose documentation can be found in documentation.pdf. The file examples.R includes some reproducible examples, where the two techniques are implemented under some simple null and alternative hypotheses. 
