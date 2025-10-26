# Notes

* State updaters should take procedure arguments. Should be able to define your own routines and pass then into the general solver framework.
* Need to think about this. 
  * Options: Finish the implementation of the associative mohr-coulomb model and then use those types
  * Option 2: Move mc_ss to using types, move mc_sr to using types, compare the implementation of ortiz_simo
* Add funciton for calculating elastic moduli
  * Might be able to convert every to G and nu and then translate from there. will have to see
  * ```G = calc_moduli("G", "moduli_1", "moduli_2")