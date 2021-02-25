# weightedClustSuite

Please Read the README_Instructions_Showcase.pdf for an introduction, documentation and examples of the package use!!

Thanks for checking out my opensource package weightedClustSuite, developed by myself.


This package was created to easily implement weighted point density clustering and validation. In my package,
one can easily give their datapoints weights and then see the results shown and validated with multiple
unsupervised clustering algorithms. This is all done through a single clustering object making the data, ex.
different clustering assignments, very easy to access.


The package operates through my implementation of manipulating the gaussian density estimates and distance
matrices with the weights, prior to the bulk of each clustering algorithm. This performs more EFFICIENTLY
and can allow for much LARGER weighted datasets than the brute force alternative of manually adding each
points weight to the dataset.
