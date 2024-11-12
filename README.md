# Master Thesis
The code and reports associated with my work as part of my M2 research thesis. 

Work based on the paper: "From Elephant to Goldfish (and Back): Memory in Stochastic Volterra Processes," written by researchers: Ofelia Bonesini, Giorgia Callegaro, Martino Grasselli, and Gilles Pagès.
You can download it from the following arxiv link: https://arxiv.org/abs/2306.02708

Short description of the files and directories of interest:
- from_elephant_to_goldfish.pdf : The research paper on which the majority of my work is based.
- quadratic_rough_Heston_model.pdf : A useful paper for understanding the quadratic rough Heston model in depth.
- CPP: The C++ implementation of the numerical scheme of the stochastic process proposed as an example in the paper.
- Python : The Jupyter notebooks presenting the mathematical details as well as the implementations of the studied algorithms.

The intended reading order for the notebooks is as follows:
- fbm.ipynb: For a first example of a Volterra process.
- kernel.ipynb: To understand convolution kernels and their effects on a typical stochastic process.
- modele_financier.ipynb: Implementation of the example presented in the research paper.
- visualisation_convergence.ipynb: Extending the calculations to develop an algorithm for visualizing the convergence of the numerical scheme for the Volterra process.
