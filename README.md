# critical-soil-models
Modern fortran library of common geotechnical constitutive models. This library is being made to interface with:
- [ ] incremental-driver - For testing and calibrating the soil models
- [ ] pumat/fumat - For plotting and analyzing the results
- [ ] Have an interface that is usable from abaqus and other common numerical solvers
- [ ] Be able to integrate with other constitutive model libraries

*Warning*: This library is just starting and we'll have to see how the progress comes along. The API is subject to change and this code has not been tested.

Soil Models that we aim to include:
- [x] Linear Elastic
- [ ] Associative Mohr-Coulomb - Paper: [N/A]()
- [ ] Non-associate Mohr-Coulomb - Paper: [N/A]()
- [ ] Modified Cam Clay - Standard version - Paper: [N/A]()
- [x] Mohr-coulomb strain softening 
  - Paper: [N/A]()
- [ ] Mohr-Coulomb with strain rate effects - Paper: [N/A]()
- [ ] NorSand - Paper: [N/A]()
- [ ] PM4Sand (possibly?) - Paper: [N/A]()
- [ ] Drucker-Prager - Paper: [N/A]()

Soil models that are lower on the list:
- [ ] Original Cam-Clay - Paper: [N/A]()
- [ ] Mohr-Coulomb with Cap - Paper: [N/A]()
- [ ] Hypoplastic models - Paper: [N/A]()

The general plan:
- [ ] Collect repos that already exist
- [ ] Put them into a unified framework
- [ ] Build unit tests for them
- [ ] Integrate them with incremental driver
- [ ] Implement the ones that don't exist
- [ ] Steamline the process for future users to implement there own models
  - [ ] Expose an api for extending the abstract model (?) class so that users can use the same solver functions but be able to model there own things.
    * Honestly not sure what finite strain and small strain models in the same repo is going to look like.

The aim is to modularize the code that the same solvers can be used for the different soil models. As much as possible the solvers should be soil model agnostic.

Also where applicable the code should be parallizable. There might not be too much a need in the actual constitutive models since looks are generally small and time dependent.

Also intially the tensor types being used here and the formulas to calculate the invariants will be included in this library but once that function is stable they'll be moved into there own library. This will allow them to be tested individually and be used in additional codes more easily.

Another goal is to write all of this code in fortran. Including the final plotting and data analysis. Let's see how this goes.

## Compiling
A `fpm.toml` file is provided for compiling critical-soil-models with the [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm). To install fpm, I recommend installing it with conda. You can install conda [here](https://www.anaconda.com/docs/getting-started/miniconda/install) 

To build the program with fpm use
```
fpm build
```

This will build the program in debug mode. Which is likely what you want so that the debugger can step into your umat.

To run the unit tests:

```
fpm test
```

To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```ford ford.md```


The latest API documentation can be found [here (Not Made yet)](). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## License

The critical-soil-models source code and related files and documentation are distributed under a permissive free software [license](https://github.com/CriticalSoilModels/Incremental_Driver/LICENSE) (BSD-style).

