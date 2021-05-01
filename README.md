# A biobjective approach for finding majority-minority districts

Section 2 of the Voting Rights Act enacts creating majority-minority districts to provide minority groups with opportunities to elect representatives of their choice. However, a majority-minority district is created if at least the following conditions (“Gingles prongs'') are satisfied: (i) compactness and numerosity (say more than t% of the citizen voting-age population), (ii) political cohesion, and (iii) sufficiency of white majority vote to defeat the preferred candidate of minorities.  Furthermore, if the first condition is not satisfied, then there is no requirement on creating a majority-minority district. Although the second and third Gingles prongs are qualitative measures, the first prong can be captured by mathematical models and algorithms.

Gingle prongs are not the only criteria that should be respected. There are three more basic rules that must be considered in evaluation of districting plans:

#### 1. Population balance.
Districts must have roughly equal populations (ideal population) to reflect "one person-one vote" concept. By Section 2 of Article I of the U.S. Constitution and Karcher v. Dagget case, any variation from the ideal population must be justified by a legitimate goal; e.g., creating compact districts.
#### 2. Contiguity. 
A district is contiguous if one can move between any two points of the district without leaving it.
#### 3. Compactness. 
One of the important districting principles that the Supreme Court considers to beat racial gerrymandering practices. Although there is no accepted operational definition of compactness, there are multiple quantitative compactness scores; e.g., number of cut edges.

## How to create compact majority-minority districts?

One might be interested in finding the maximum number of majority-minority districts such that the base rules of districting (i.e., population balance, contiguity, and compactness) are respected. Although I am a big fan of Mixed Integer Programming approaches to solve this biobjective optimization problem, I prefer to use [GerryChain](https://gerrychain.readthedocs.io/en/latest/) for simplicity purposes. [GerryChain](https://gerrychain.readthedocs.io/en/latest/) employs [Markov chain Monte Carlo (MCMC)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo) for creating feasible maps. GerryChain uses two main strategies to move from a feasible solution to another one:

#### 1. Flip.
A strategy in which a node on the boundary of a district is assigned to an adjacent district.

#### 2. Recommbination.
A strategy in which (i) d adjacent districts (d >= 2) are unified; (ii) a spanning tree is built on the unified nodes; and (iii) d-1 edges of the spanning trre are removed to creat d new districts. 

In this experiment, I ran GerryChain for 10,000 iterations (as suggested by members of [MGGG](https://mggg.org)) on some county-level and tract-level instances. Further, I considered a population deviation of 1% and a minority population threshold of t = 30%. Dual graphs and shapefiles are borrowed from [Daryl DeFord's website](https://people.csail.mit.edu/ddeford/dual_graphs.html) and [Eugene Lykhovyd's website](https://lykhovyd.com/files/public/districting), respectively.

Finally, I found the [Pareto Frontier](https://en.wikipedia.org/wiki/Pareto_efficiency) that represents the "best" plans with respect to number of cut edges and number of majority-minority districts. 

![Figure 1](heur_MS_county_pareto.png?raw=true "heur_MS_county_pareto")![Figure 2](heur_MS_tract_pareto.png?raw=true "heur_MS_tract_pareto")





