# metaSPPstu: solving the Set Packing Problem (student version)
Implementation in Julia (compliant Julia v1.x) of tools related to the Set Packing Problem (SPP) for pedagogical purposes.

## Le projet
Exercices d'implémentations réalisés dans le cadre du cours de métaheuristiques.

## Membres du binôme
- Le Normand Léna
- Sachot Clara

------

## Le répertoire

### Data
Ce dossier contient une instance "test", didactic.dat, et une sélection de 10 innstances numériques de SPP au format OR-libraryn choisies parmis une [lste d'instances disponibles en ligne](https://www.emse.fr/~delorme/SetPackingFr.html)

 - didactic.dat
 - pb_100rnd0100.dat
 - pb_100rnd0300.dat
 - pb_200rnd0100.dat
 - pb_200rnd0500.dat
 - pb_200rnd1600.dat
 - pb_500rnd0100.dat
 - pb_500rnd1700.dat
 - pb_1000rnd0100.dat
 - pb_1000rnd0200.dat
 - pb_2000rnd0100.dat


Ce sont sur ces instances que nous avons testé nos programmes.

### Documentation
Ce dossier contient la version PDF du rapport de ce projet.

### Codes sources
Lorsqu'une fonction demande un fname, il faut la lancer avec le chemin d'une instance de fichier, exemple : `greedy("../Data/pb_100rnd0300.dat")`


 **EI1.jl**
  - greedy(fname) : construction gloutonne
  - resoudreSPP(fname) : construction avec descente plus profonde en kp échange
  - experimentationSPP() : résoudre le dadactic.dat + les 10 instances de SPP

**EI2.jl**
 - grasp(fname, alpha, repeat) : construction GRASP avec descente. Les paramètres alpha et repeat ont une valeur par défaut à 0.7 et 12 respectivement
  - ReactiveGrasp(fname, m, maxIteration, Nalpha) : résoudre une instance avec Reactive GRASP. Les paramètres m, maxIteration et Nalpha ont une valeur par défaut de respectivement 10, 85 et 15.
  - experimentationSPP() : résoudre le dadactic.dat + les 10 instances de SPP en utilisant Reactive GRASP

**EI3.jl**
 -  genetic_algorithm(fname, population_size, generations, mutation_rate, crossover_rate) : résoudre une instance avec l'algorithme génétique. Les paramètres population_size, generations, mutation_rate et crossover_rate ont une valeur par défaut de respectivement 200, 50, 0.9 et 0.01.
 - experimentationSPP() : résoudre le dadactic.dat + les 10 instances de SPP en utilisant l'algorithme génétique 

 **main.jl**
  - solveGLPK(fname) : résoudre l'instance avec GLPK

#### Ce qu'il faut pour les lancer
**main.jl**
 - `using GLPK`
 - `include("loadSPP.jl")`
 - `include("setSPP.jl")`

**Chaque fichier demande**
 - `include("loadSPP.jl")`
 - `using LinearAlgebra`

 **Spécifique à EI2.jl**
 - `include("EI1.jl")`
 - `using LinearAlgebra`
 - `using InvertedIndices`
 - `using Distributions`
 - `using StatsBase`

**Spécifique à EI3.jl**
- `include("EI2.jl")`






- `loadSPP.jl` : lecture d'une instance de SPP au format OR-library
- `setSPP.jl` : construction d'un modèle JuMP de SPP
- `experiment.jl`: protocole pour mener une expérimentation numérique avec sorties graphiques


------
