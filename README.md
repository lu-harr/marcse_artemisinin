### A repo for models of prevalence ofr molecular markers of artemisinin and partner drug resistance

#### What the code does

There should (eventually) be (at least) three models in here. For now, there is a model of kelch13 mutation 
prevalence only. This model takes as input published data from the WWARN artemisinin resistance surveyor 
dataset [include link] of prevalence (number of people tested **and** number of people with mutant *P. falciparum* 
parasites) of WHO or WWARN-validated or -associated markers for artemisinin resistance in the propeller region
of the kelch13 gene of *P. falciparum* malaria.

We model prevalence as a Gaussian Process with a binomial observation model and Circular Matern kernel in 
space [inc. refs]. [Say something about uncertainty]

#### What's in here

#### See also

- Flegg et al., 2022
- Flegg et al., 2024
- Foo et al., 2024

#### TODOs

- grab data from K13 surveyor
- write down equations for existing binomial model and for extensions
- do some summaries for myself of the number/spread of the different mutations
- I'm all about getting this on the HPC
- put some links into README
- think about partner drugs
- think about how partner drug outputs combine with K13 outputs to make "overall ACT molecular resistance map"


