# polymer-stats
Numerical modelling of the electroelasticity of polymer chains consisting of dielectric monomers and polar monomers.
Chains consisting of dielectric monomers can be simulated in either the fixed end-to-end vector (`mean_field.jl`) or fixed force ensembles (`mcmc_eap_chain.jl`).
The fixed end-to-end vector ensemble is simulated using mean-field theory and numerical integration; whereas the fixed force ensemble is simulated using the Markov chain Monte Carlo method.
As of now, chains consisting of polar monomers and dipole-dipole interactions (for either chain type) can only be simulated in the fixed force ensemble (`mcmc_eap_chain.jl`). The mean-field theory simulations tend to be faster, but more physics can be captured using the Markov chain Monte Carlo method.

## Getting started
This research code is written in the [julia language](https://julialang.org/).
In the ``scripts`` directory, there is a script for installing package dependences:
    
    julia scripts/install_dependencies.jl
    
The main file for the mean-field theory simulations is ``mean_field.jl``.
You can see options by running:

    julia mean_field.jl --help
    
The theoretical details of this code can be found in the first reference (in particular, sections 3 and 4.1) in the Citations section below.

The main file for the Markov chain Monte carlo method is ``mcmc_eap_chain.jl``.
You can see options by running:

    julia mcmc_eap_chain.jl --help

An experimental algorithm for clustering of the chain is implmeneted in ``mcmc_clustering_eap_chain.jl``.
You can see options by running:

    julia mcmc_clustering_eap_chain.jl --help
    
It has been verified against analytical results, so I am relatively confident in its correctness.
The expectation is that the clustering algorithm will vastly improve the convergence rate for dielectric chains with monomer-monomer interactions, but this has not been studied yet.

## Questions
Documentation of this work is still in progress. However, I am very receptive to questions and comments via email. You can email questions to grasingerm@gmail.com.

## Citations
There are a few papers associated with polymer-stats. 
If you found this code useful for your research, please be kind enough to cite it, using the DOIs 10.1039/D0SM00845A, 10.1073/pnas.2102477 and 10.1016/j.jmps.2021.104658, or the following BiBTeX entries:


    @article{grasinger2020statistical,
        title={Statistical mechanical analysis of the electromechanical coupling in an electrically-responsive polymer chain},
        author={Grasinger, Matthew and Dayal, Kaushik},
        journal={Soft Matter},
        volume={16},
        number={27},
        pages={6265--6284},
        year={2020},
        publisher={Royal Society of Chemistry},
        doi={10.1039/D0SM00845A}
    }
    
    @article {grasinger2021flexoelectricity,
        author = {Grasinger, Matthew and Mozaffari, Kosar and Sharma, Pradeep},
        title = {Flexoelectricity in soft elastomers and the molecular mechanisms underpinning the design and emergence of giant flexoelectricity},
        volume = {118},
        number = {21},
        year = {2021},
        doi = {10.1073/pnas.2102477118},
        publisher = {National Academy of Sciences},
        issn = {0027-8424},
        journal = {Proceedings of the National Academy of Sciences}
    }

    @article{grasinger2021statistical,
        title={Statistical mechanics of a dielectric polymer chain in the force ensemble},
        author={Grasinger, Matthew and Dayal, Kaushik and deBotton, Gal and Purohit, Prashant K},
        journal={Journal of the Mechanics and Physics of Solids},
        pages={104658},
        year={2021},
        publisher={Elsevier}
    }
