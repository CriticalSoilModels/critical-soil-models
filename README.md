# critical-soil-models
Modern fortran library of common geotechnical constitutive models. This library is being made:
- [x] Interface with  - For testing, calibrating, and debugging their constitutive model.

- [ ] Expose an api that let's users develop and test there consitutive models extremely quickly. Users should be able to use predifined types and through that solvers to quickly test test and modifiy there constitutive relations.

- [ ] Have an interface that is usable from:
  - [ ] abaqus 
  - [ ] Other common numerical solvers (Open an issue and suggest an interface)

- [ ] Be able to integrate with other constitutive model libraries

*Warning*: This library is just starting and we'll have to see how the progress comes along. The API is subject to change and this code has not been tested.

## Design Philosophy
Geotechnical engineers spend way too much time looking for, reimplementing, testing, debugging, constitutive relations that already exist. This sucks time away from running the larger simulation, actually getting results, sleeping/graduating/etc.

Hence, there are two overarching goals for this repo. The first is to collect, abstract and accelerate the development of constitutive models. The second, is to provide a library disconnected from any one PDE solver. 

As pointed out in (*Kayenta Paper Add Here*) it is extremly difficult to figure out if a bug is caused by your constitutive model or the numerical technique, when the constitutive strongly embeded in the PDE solver. A split framework allows the constitutive relations to crows source testing and validation to a much higher degree. 

## Soil Models that we aim to include:
(Assme the model is in 3D unless otherwise specified)
- [x] **Linear Elastic**
  - [ ] Paper
  - [ ] Unit Tested
  - [ ] Abstracted

- [ ] **Bingham**
  - [ ] Paper
  - [ ] Unit Tested
  - [ ] Abstracted

- [ ] **Modified Cam Clay - Standard version**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted

- [x] **Mohr-coulomb strain softening**
  * This model can act as an associative (yield and plastic surface are the same) or non-associate model (yield and plastic surface are independent) depending on the input parameters.
  - [ ] Paper: The repo was originally written by A. Yerro-Colom (use [Yerro (2015)](https://upcommons.upc.edu/handle/2117/102412) )
  - [ ] Unit Tested
  - [ ] Abstracted

- [x] **Mohr-Coulomb with strain rate effects**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted
  
- [ ] **NorSand**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted
  
- [ ] **PM4Sand (2D - Plain Strain)**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted

- [ ] **Drucker-Prager**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted

Soil models that are lower on the list:
- [ ] **Original Cam-Clay**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted
  
- [ ] **Mohr-Coulomb with Cap**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted
  
- [ ] **Hypoplastic models**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted

## The general plan:
- [x] Integrate them with [incremental-driver](https://github.com/CriticalSoilModels/Incremental_Driver)
  * Allows to constitutive models to be tested and calibrated in a single element framework
  * Turns out the first steps of this were really easy thanks to ```fpm```
- [ ] Collect repos that already exist
- [ ] Put them into a unified framework
- [ ] Build unit tests for them
- [ ] Implement the ones that don't exist
- [ ] Steamline the process for future users to implement there own models
  - [ ] Expose an api for extending the abstract model (?) class so that users can use the same solver functions but be able to model there own things.
- [ ] Users should be able to generate, run, and analyze data all in fortran. Why? Fortran is goated :sunglasses:
- [ ] Use the tensor types from ttb to ease calculations
- [ ] Add interface for C++/C/python/Julia to be able to be able to call the constitutive relations

## Citations

If you use this repo in your research or a publish manuscript please cite this library and the paper the implementation of the soil model was taken from.

Please refer to the Citation.cff file for the citation inforamtion for this repo.

## License

The critical-soil-models source code and related files and documentation are distributed under the [GNU GENERAL PUBLIC LICENSE](https://github.com/CriticalSoilModels/critical-soil-models/blob/main/LICENSE).

## Other Notes
There are other libraries of constitutive models out there. The ones that I've come across are:
* (Write list here)

The plan is to include, licenses and time permitting, as many of the implementations in this library.

Soil consitutive models are generally based on a specific set of stress and strain invariants. Currently, the implementation of the invariants are included in this library. In the future, they'll be moved to there own repo (maybe?) and given a python interface. The values should be compared to the Julia library that implements them with automatic differentiation (Fortran not having built in AD is not elite. Hey you, yeah you, drop out of geotech and go implement Automatic differention in the Lfortran compiler. Yeah that's a good idea.) This way users won't have to implement the invariants again inevitablly making a mistake. Instead we all have the same errors together. Yay!.