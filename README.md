# polymer-stats
Numerical modelling of the electroelasticity of polymer chains consisting of dielectric monomers and polar monomers.
Chains consisting of dielectric monomers can be simulated in either the fixed end-to-end vector (`mean_field.jl`) or fixed force ensembles (`mcmc_eap_chain.jl`).
The fixed end-to-end vector ensemble is simulated using mean-field theory and numerical integration; whereas the fixed force ensemble is simulated using the Markov chain Monte Carlo method.
As of now, chains consisting of polar monomers and dipole-dipole interactions (for either chain type) can only be simulated in the fixed force ensemble (`mcmc_eap_chain.jl`). The mean-field theory simulations tend to be faster, but more physics can be captured using the Markov chain Monte Carlo method.

## Getting started
This research code is written in the [julia language](https://julialang.org/).
The main file for the mean-field theory simulations is ``mean_field.jl``.
You can see options by running:

    julia mean_field.jl --help
    
The main file for the Markov chain Monte carlo method is ``mcmc_eap_chain.jl``.
You can see options by running:

    julia mcmc_eap_chain.jl --help

## Questions
Documentation of this work is still in progress. However, I am very receptive to questions and comments via email. You can email questions to grasingerm at gmail dot com.

## Citations
There are a few papers associated with polymer-stats. 
If you found this code useful for your research, please be kind enough to cite it, using the DOIs 10.1039/D0SM00845A and 10.1073/pnas.2102477, or the following BiBTeX entries:


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
    
    @article{grasinger2021flexoelectricity,
      title={Flexoelectricity in soft elastomers: molecular mechanisms, materials design and the emergence of giant flexoelectricity},
      author={Grasinger, Matthew and Mozaffari, Kosar and Sharma, Pradeep},
      journal={Proceedings of the National Academy of Sciences},
      year={in-press},
      publisher={National Acad Sciences},
      doi={10.1073/pnas.2102477}
    }

    @phdthesis{grasinger2019multiscale,
      title={Multiscale Modeling and Theoretical Design of Dielectric Elastomers},
      author={Grasinger, Matthew},
      year={2019},
      school={Carnegie Mellon University}
    }
