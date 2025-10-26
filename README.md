# critical-soil-models
Modern fortran library of common geotechnical constitutive models. This library is being made:
- [x] Interface with [incremental-driver](https://github.com/CriticalSoilModels/Incremental_Driver) - For testing, calibrating, and debugging their constitutive model.

- [ ] Expose an api that let's users develop and test there consitutive models extremely quickly. Users should be able to use predifined types and through that solvers to quickly test test and modifiy there constitutive relations.

- [ ] Have an interface that is usable from:
  - [ ] abaqus 
  - [ ] Other common numerical solvers (Open an issue and suggest an interface)

- [ ] Be able to integrate with other constitutive model libraries

*Warning*: This library is just starting and we'll have to see how the progress comes along. The API is subject to change and this code has not been tested.

## Soil Models that we aim to include:
(Assme the model is in 3D unless otherwise specified)
- [x] **Linear Elastic**
  - [ ] Paper
  - [ ] Unit Tested
  - [ ] Abstracted
    
- [ ] **Modified Cam Clay - Standard version**
  - [ ] Paper: 
  - [ ] Unit Tested
  - [ ] Abstracted

- [x] **Mohr-coulomb strain softening**
  * This model can act as an associative (yield and plastic surface are the same) or non-associate model (yield and plastic surface are independent) depending on the input parameters.
  - [ ] Paper: 
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
- [ ] Collect repos that already exist
- [ ] Put them into a unified framework
- [ ] Build unit tests for them
- [ ] Integrate them with incremental driver
- [ ] Implement the ones that don't exist
- [ ] Steamline the process for future users to implement there own models
  - [ ] Expose an api for extending the abstract model (?) class so that users can use the same solver functions but be able to model there own things.


## Design Philosophy
The aim is to modularize the code that the same solvers can be used for the different soil models. As much as possible the solvers should be soil model agnostic.

Also where applicable the code should be parallizable (if the user wants to)

Also intially the tensor types being used here and the formulas to calculate the invariants will be included in this library but once that function is stable they'll be moved into there own library. This will allow them to be tested individually and be used in additional codes more easily.

Another goal is to write all of this code in fortran. Including the final plotting and data analysis. Let's see how this goes.


## License

The critical-soil-models source code and related files and documentation are distributed under the [GNU GENERAL PUBLIC LICENSE](https://github.com/CriticalSoilModels/critical-soil-models/LICENSE).

