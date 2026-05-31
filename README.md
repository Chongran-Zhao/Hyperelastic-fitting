# Hyperelastic Fitting

MATLAB code for calibrating incompressible hyperelastic material models from
experimental deformation data. The project now uses one shared driver and one
shared source tree instead of separate per-model script folders.

## Project Layout

- `Driver.m` - main example script for selecting data sets, material models,
  bounds, fitting, and prediction plots.
- `src/` - shared fitting, data loading, objective, plotting, and model code.
- `src/material_models/` - built-in constitutive models.
- `src/tools/` - tensor, kinematics, quadrature, evaluation, and plotting
  helpers.
- `data/` - experimental data and the Lebedev quadrature table used by the
  non-Gaussian chain-network model.

Legacy model-specific folders have been removed. New work should be routed
through `Driver.m` and the shared `src/` implementation.

## Requirements

- MATLAB
- Optimization Toolbox, for `lsqnonlin`

No external MATLAB packages are required.

## Quick Start

Open MATLAB in the repository root and run:

```matlab
Driver
```

Edit `Driver.m` to choose the fitting data, material model, parameter bounds,
and optional prediction data. A typical workflow is:

```matlab
fitting_cases = {};
fitting_cases = Add_exp_data_sets(fitting_cases, 'James', 'UT');

models = [];
models = Add_material_model(models, ...
    Zhan_NonGaussian([1.0, 100.0]), ...
    [0.0, 0.0], [Inf, Inf]);

[models, fit] = start_fit(models, fitting_cases);
plot_simultaneous_fit(models, fitting_cases);
```

Call `Add_exp_data_sets()` with no inputs to print the registered data-set
calls available in the current code.

## Available Models

The built-in model constructors live in `src/material_models/`:

- `Neo_Hookean`
- `Mooney_Rivlin`
- `Yeoh`
- `Ogden`
- `Arruda_Boyce`
- `Hill_SH`
- `Hill_Hencky`
- `Hill_CR`
- `Hill_CZ`
- `Hill_DN`
- `Zhan_Gaussian`
- `Zhan_NonGaussian`

Models can be combined by calling `Add_material_model` multiple times. Bounds
are attached to each added model and then assembled into one fitting vector.

## Data Sets

Experimental data are stored under `data/`. Registered data include uniaxial
tension (UT), equibiaxial tension (ET), pure shear (PS), uniaxial extension
(UE), and biaxial tension (BT) cases from:

- Treloar (1944)
- Kawabata et al. (1981)
- Meunier et al. (2008)
- Kawamura et al. (2001)
- Jones and Treloar (1975)
- James, Green, and Simpson (1975)
- Katashima et al. (2012)

To add a new data set, place the files under `data/` and register the case in
`src/Add_exp_data_sets.m`.

## Extending

To add a material model, create a constructor in `src/material_models/` that
returns a struct with these fields:

- `name`
- `parameters`
- `parameter_names`
- `energy`
- `S`
- `P`
- `set_parameters`

Then add it in `Driver.m` with `Add_material_model`.

## References

- Dal, H., Açıkgöz, K., & Badienia, Y. (2021). *On the Performance of
  Isotropic Hyperelastic Constitutive Models for Rubber-Like Materials: A State
  of the Art Review.* **ASME Applied Mechanics Reviews**, 73(2), 020802.
- Ogden, R., Saccomandi, G., & Sgura, I. (2004). *Fitting hyperelastic models
  to experimental data.* **Computational Mechanics**, 34, 484-502.
- Xiang, Y., Zhong, D., Rudykh, S., Zhou, H., Qu, S., & Yang, W. (2020). *A
  Review of Physically Based and Thermodynamically Based Constitutive Models
  for Soft Materials.* **ASME Journal of Applied Mechanics**, 87(11), 110801.
- Holzapfel, G. A. (2002). *Nonlinear Solid Mechanics: A Continuum Approach for
  Engineering Science.*
- Liu, J., Guan, J., Zhao, C., & Luo, J. (2024). *A Continuum and Computational
  Framework for Viscoelastodynamics: III. A Nonlinear Theory.* **Computer
  Methods in Applied Mechanics and Engineering**, 430, 117248.
- Zhan, L., Wang, S., Qu, S., Steinmann, P., & Xiao, R. (2023). *A new
  micro-macro transition for hyperelastic materials.* **Journal of the
  Mechanics and Physics of Solids**, 171, 105156.

## Author

- Chongran Zhao
- Southern University of Science and Technology, China
- [chongranzhao@outlook.com](mailto:chongranzhao@outlook.com)
