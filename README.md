# finite-square-well-demo

This is a quick and dirty demo of the pseudopotential concept using 1D finite well quantum scattering problem. It follows discussion in Richard Martin's Electronic Structure textbook, problem 11.2. See also Scherrer's QM book, section 4.1.

The core feature of this demo is an animation of a plane wave scattering past a finite square well or barrier, with the real part shown in blue and the imaginary part shown in green. The well depth (*a*) and width (*V0*) can be adjusted, as well as the energy (*E*) of the plane wave. The plane waves are normalized so the magnitude of the incident plane wave is 1.

Pseudopotential concept:
Given a well of width *a* and height *V0* (barrier if postive, well if negative), other values of *V0* and *a* can be found that do not alter the transmission coefficient (transmitted/incident). So we have a \*pseudo\* potential that reproduces the true potential's scattering properties away from the well or barrier. The same concept applies to pseudopotentials in the atomic case, where scattering properties outside some core radius need to match the scattering from the true atomic potential.

## Getting Started

### Prerequisites

Requires Python 2, matplotlib, scipy, numpy

### Running demo with default parameters:

```
python test.py
```

This should open several windows running the scattering animation. Notice the well depths are different but the transmission is the same in each case. Often the wells are too deep to fit  in the animation window, but the relative depth can be inferred from the wavelength of the plane wave inside the well, with deeper wells causing shorter wavelengths.

### Running with your own parameters

test.py can be modified to run the demo with different parameters.

Edit a, V0 to change well shape. Edit E to change the incident plane wave energy.

The arguments to find\_some\_pseudopotentials are Vmin, Vmax, and steps, which control how the code searches for pseudo values of V0. Decreasing Vmin will tend to find lower values of V0. Increasing Vmax will tend to find higher values of V0 (although there may be an upper limit). Increasing steps will tends to increase the number of pseudo V0 values found. What these actually control are guesses fed in to the broyden1 zero finder, which are calculated as np.linspace(Vmin, Vmax, steps).

### Known Issues

Running a large number of animations simultaneously will bog down computer because each animation is run in a separate thread. Each animation is called using subprocess, which isn't ideal but was the quickest way to run an arbitrary number of simultaneous animations in separate windows.


## Authors

* **Casey Brock** - *Initial work* - [caseynbrock](https://github.com/caseynbrock)



## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



## Acknowledgements

* UC Boulder PHET's Quantum Tunneling and Wave Packets Demo
