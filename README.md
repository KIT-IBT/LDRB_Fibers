# Adapted LDRB method for ventricular fibers

This is a MATLAB implementation of the Laplace–Dirichlet Rule-Based (LDRB) algorithm for generating ventricular fibers. The original algorithm [[1](#1)] was adapted to eliminate a discontinuity of fibers in the free walls and to yield a fiber rotation that is directly proportional to the transmural Laplace solution (approximately linear across the wall).
&nbsp;  
&nbsp;  
![Fibers](https://user-images.githubusercontent.com/31965820/107220831-1b4ebd00-6a13-11eb-97a3-c883126692b7.jpg)

## Dependencies

Run [dependencies/install.sh](dependencies/install.sh) to install the dependencies [[2](#2),[3](#3)].

## Example

Run [examples/example.m](examples/example.m) to generate fibers for an exemplary biventricular geometry.

## License

All source code is subject to the terms of the GPL-3.0 License.  
Copyright 2021 Steffen Schuler, Karlsruhe Institute of Technology.

## Contact

Steffen Schuler  
Institute of Biomedical Engineering  
Karlsruhe Institute of Technology  
www.ibt.kit.edu

## References

<a id="1">[1]</a> [Bayer, J. et al., 2012. A novel rule-based algorithm for assigning myocardial fiber orientation to computational heart models. Ann Biomed Eng.](https://doi.org/10.1007/s10439-012-0593-5)  
<a id="2">[2]</a> [Jacobson, A., 2018. gptoolbox: Geometry processing toolbox.](https://github.com/alecjacobson/gptoolbox)  
<a id="3">[3]</a> [Schuler, S., 2020. vtkToolbox: A MEX interface to the VTK library.](https://github.com/KIT-IBT/vtkToolbox)  
