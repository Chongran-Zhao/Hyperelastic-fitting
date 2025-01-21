# Hyperelastic Material Deformation Calibration Code

## Overview
This repository contains code used to **calibrate** the deformation of hyperelastic materials. The calibration process employs two basic algorithms:

- For models with parameters **≤ 2**, only uniaxial tensile tests are used to calibrate and obtain optimized parameters. Experimental data from other tests are used for prediction and validation.  
  This approach is implemented in the `micromechanical_model/` directory.

- For models with parameters **≥ 3**, a single loading condition (e.g., uniaxial tensile) is insufficient to obtain a unique parameter set. Therefore, we calibrate the model using multiple sets of experimental data, known as **simultaneous fitting**.  
  This approach is implemented in the `AB/`, `Ogden/`, `Hill/`, and `Yeoh/` directories.

For more detailed discussions on the calibration of different hyperelastic models, please refer to:

- Dal, H., Açıkgöz, K., and Badienia, Y. (2021). "On the Performance of Isotropic Hyperelastic Constitutive Models for Rubber-Like Materials: A State of the Art Review." *ASME Applied Mechanics Reviews*, 73(2): 020802.
- Ogden, R., Saccomandi, G., & Sgura, I. (2004). "Fitting hyperelastic models to experimental data." *Computational Mechanics*, 34, 484–502.
- Xiang, Y., Zhong, D., Rudykh, S., Zhou, H., Qu, S., and Yang, W. (2020). "A Review of Physically Based and Thermodynamically Based Constitutive Models for Soft Materials." *ASME Journal of Applied Mechanics*, 87(11): 110801.

---

## Author

**Chongran Zhao**  
*Southern University of Science and Technology, China*  
*Email: [chongranzhao@outlook.com](mailto:chongranzhao@outlook.com)*  
*Date: Jan. 21, 2024*

---

## Table of Contents
- [Experimental Data](#experimental-data)
- [Continuum Basis](#continuum-basis)
- [Calibration Details and Material Models](#calibration-details-and-material-models)
- [Hill's Hyperelastic Model with Generalized Strains](#hills-hyperelastic-model-with-generalized-strains)
- [Evaluation Functions (NMAD, MSD, R²)](#evaluation-functions-nmad-msd-r²)

---

## References

### Experimental Data
We have collected several sets of experimental data from various research papers, including:

- Treloar, L. R. (1944). *Stress-strain data for vulcanized rubber under various types of deformation.* *Rubber Chemistry and Technology*, 17(4), 813-825.
- Kawabata, S., Matsuda, M., Tei, K., & Kawai, H. (1981). "Experimental survey of the strain energy density function of isoprene rubber vulcanizate." *Macromolecules*, 14, 154-162.
- Meunier, L., Chagnon, G., Favier, D., Orgéas, L., & Vacher, P. (2008). "Mechanical experimental characterisation and numerical modelling of an unfilled silicone rubber." *Polymer Testing*, 27, 765-777.
- Jones, D. F., & Treloar, L. R. G. (1975). "The properties of rubber in pure homogeneous strain." *Journal of Physics D: Applied Physics*, 8, 1285–1304.
- Kawamura, T., Urayama, K., & Kohjiya, S. (2001). "Multiaxial deformations of end-linked poly(dimethylsiloxane) networks. 1. Phenomenological approach to strain energy density function." *Macromolecules*, 34(23), 8252-8260.
- James, A. G., Green, A., & Simpson, G. M. (1975). "Strain energy functions of rubber. I. Characterization of gum vulcanizates." *Journal of Applied Polymer Science*, 19(7), 2033-2058.
- Katashima, T., Urayama, K., Chung, U. I., & Sakai, T. (2012). "Strain energy density function of a near-ideal polymer network estimated by biaxial deformation of Tetra-PEG gel." *Soft Matter*, 8(31), 8217-8222.

---

### Continuum Basis
For a comprehensive understanding of the **continuum basis**, please refer to:

- Holzapfel, G. A. (2002). *Nonlinear Solid Mechanics: A Continuum Approach for Engineering Science.*

The strain energy is decomposed into isochoric and volumetric parts as follows:
$$
\Psi(\mathbf{C}) = \Psi_\mathrm{ich}(\tilde{\mathbf{C}}) + \Psi_\mathrm{vol}(J),
$$
For incompressible materials, a Legendre transformation on \( J \) converts the Helmholtz free energy into Gibbs free energy:
$$
G(\mathbf{C}) = G_\mathrm{ich}(\tilde{\mathbf{C}}) + G_\mathrm{vol}(P).
$$
For more details, see:

- Liu, J., & Marsden, A. L. (2018). "A unified continuum and variational multiscale formulation for fluids, solids, and fluid–structure interaction." *Computer Methods in Applied Mechanics and Engineering*, 337, 549-597.

---

### Calibration Details and Material Models
For insights into **calibration details** and a review of **material models**, see:

- Dal, H., Açıkgöz, K., & Badienia, Y. (2021). *On the Performance of Isotropic Hyperelastic Constitutive Models for Rubber-like Materials: A State of the Art Review.* *Applied Mechanics Reviews*, 73(2), 020802.

---

### Hill's Hyperelastic Model with Generalized Strains
Detailed information on the **Hill's hyperelastic model** with **generalized strains** can be found in:

- Liu, J., Guan, J., Zhao, C., & Luo, J. (2024). *A Continuum and Computational Framework for Viscoelastodynamics: III. A Nonlinear Theory.* *Computer Methods in Applied Mechanics and Engineering*, 430, 117248. [DOI: 10.1016/j.cma.2024.117248](https://doi.org/10.1016/j.cma.2024.117248)

---

### A New Micro-Macro Transition for the Micro-Mechanical Model
I have reproduced part of the results from:

- Zhan, L., Wang, S., Qu, S., Steinmann, P., & Xiao, R. (2023). "A new micro–macro transition for hyperelastic materials." *Journal of the Mechanics and Physics of Solids*, 171, 105156.

These results are included in the `micromechanical_model/` directory, along with some experimental data. Other models are calibrated simultaneously, such as those in the `AB/`, `Hill/`, `Ogden/`, and `Yeoh/` directories.

---

### Evaluation Functions (NMAD, MSD, R²)
For an explanation of the **evaluation functions** used (NMAD, MSD, R²), please visit:

- [Fitness Functions Used by MCalibration](https://polymerfem.com/fitness-functions-used-by-mcalibration/)

---

*Please note that the links provided are for reference only and may require a stable internet connection for access. If you encounter any issues with the links, please verify their legitimacy and try again later.*

*I'm glad to communicate and cooperate with you if you have any questions. Feel free to contact me!*