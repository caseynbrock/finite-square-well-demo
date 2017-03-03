# finite-square-well-demo

quick and dirty demo of pseudopotential concept using 1D finite well/barrier scattering problem. Follows discussion in Richard Martin Electronic Structure textbook problem ???. See also Scherrer problem ???.

The main feature of this demo is an animation of a plane wave scattering past a finite square well. The well depth (a) and width (V0) can be adjusted. Normalized to ???. The transmission coefficient is given by ???. For the animation, the incident wave is normalized to ???. Full solution given by ???.

Pseudopotential concept:
Given a well of width a and depth V0, other values of V0 and a can be found that do not alter the transmission coefficient (transmitted/incident). So we have a \*pseudo\* potential that essentially creates the same scattering scattering properties outside the well. The same concept applies to pseudopotentials in the atomic case, where scattering properties outside some core radius need to match the all-electron case.



## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Requires Python 2, matplotlib, scipy, numpy

### Running example

To run with default paramters:

```
Python test.py
```

This should open several windows running the scattering animation. Notice the well depths are different but the transmission is the same in each case.

### Running with your own paramters

Test.py can be modified to run the demo with different parameters.

Edit a, V0 to change well shape. Edit E to change the incident plane wave energy.

The arguments to find\_some\_pseudopotentials are Vmin, Vmax, and steps, which control how the code searches for pseudo values of V0. Decreasing Vmin will tend to find lower values of V0. Increasing Vmax will tend to find higher values of V0 (although there may be an upper limit). Increasing steps will tends to increase the number of pseudo V0 values found. What these actually control are guesses fed in to the broyden1 zero finder, which are calculated as np.linspace(Vmin, Vmax, steps).

### Known Issues

Running a large number of animations simultaneously will bog down computer because each is run in a separate thread. Each animation is called using subprocess, which isn't ideal but was the quickest way to run an arbitrary number of simultaneous animations in separate windows.


## Authors

* **Casey Brock** - *Initial work* - [caseynbrock](https://github.com/caseynbrock)



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



## Acknowledgments

* Richard Martin's Electronic Structure text book
* UC Boulder PHET's Quantum Tunneling and Wave Packets Demo
