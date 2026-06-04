# Hyperelastic Fitting

MATLAB code for calibrating incompressible hyperelastic material models from
experimental deformation data. The project now uses one shared driver and one
shared implementation layout instead of separate per-model script folders.

## Project Layout

- `driver.m` - main example script for selecting data sets, material models,
  bounds, fitting, and prediction plots.
- `src/` - shared fitting, data loading, objective, plotting, and model code.
- `material_models/` - built-in constitutive models.
- `tools/` - tensor, kinematics, quadrature, evaluation, and plotting helpers.
- `data/` - experimental data and the Lebedev quadrature table used by the
  non-Gaussian chain-network model.

Legacy model-specific folders have been removed. New work should be routed
through `driver.m` and the shared implementation folders.

## Requirements

- MATLAB
- Optimization Toolbox, for `lsqnonlin`

No external MATLAB packages are required.

## Quick Start

Open MATLAB in the repository root and run:

```matlab
driver
```

Edit `driver.m` to choose the fitting data, material model, parameter bounds,
and optional prediction data. A typical workflow is:

```matlab
fitting_cases = {};
fitting_cases = add_exp_data_sets(fitting_cases, 'James', 'UT');

models = [];
models = add_material_model(models, ...
    Zhan_NonGaussian([1.0, 100.0]), ...
    [0.0, 0.0], [Inf, Inf]);

[models, fit] = start_fit(models, fitting_cases);
plot_simultaneous_fit(models, fitting_cases);
```

Call `add_exp_data_sets()` with no inputs to print the registered data-set
calls available in the current code.

## Available Models

The built-in model constructors live in `material_models/`. All model
constructors return a struct with energy, second Piola-Kirchhoff stress, first
Piola-Kirchhoff stress, parameter names, and a `set_parameters` callback used
by the fitting code.

Classical incompressible hyperelastic models:

| Constructor | Parameters | Notes |
| --- | --- | --- |
| `Neo_Hookean` | `[mu]` | Isochoric first-invariant model. |
| `Mooney_Rivlin` | `[C1, C2]` | Isochoric first- and second-invariant model. |
| `Yeoh` | `[C1, C2, C3]` | Cubic polynomial in `I1_bar - 3`. |
| `Ogden` | `[mu1, alpha1, mu2, alpha2, ...]` | Any number of `[mu, alpha]` pairs. |
| `Arruda_Boyce` | `[mu, N]` | Eight-chain model using the inverse Langevin approximation. |

Zhan micro-macro transition models:

| Constructor | Parameters | Notes |
| --- | --- | --- |
| `Zhan_Gaussian` | `[mu]` | Closed-form Gaussian chain-network model based on the principal stretches of `U_bar`. |
| `Zhan_NonGaussian` | `[mu, N]` | Non-Gaussian chain-network model evaluated by Lebedev sphere quadrature. |

Guan Gaussian-chain micro-macro models:

| Constructor | Chain strain `E_hat` | Parameters |
| --- | --- | --- |
| `Guan_Gaussian_GL` | Green-Lagrange, `1/2*(lambda_n^2 - 1)` | `[mu, strain parameters...]` |
| `Guan_Gaussian_Biot` | Biot, `lambda_n - 1` | `[mu, strain parameters...]` |
| `Guan_Gaussian_SH2_3` | Seth-Hill `m = 2/3`, `3/2*(lambda_n^(2/3) - 1)` | `[mu, strain parameters...]` |
| `Guan_Gaussian_Hencky` | Hencky, `log(lambda_n)` | `[mu, strain parameters...]` |
| `Guan_Gaussian_SH` | Seth-Hill, `(lambda_n^m_hat - 1)/m_hat` | `[mu, m_hat, strain parameters...]` |

For the Guan models, the constructor suffix specifies the chain-scale strain
`E_hat`. The second constructor argument selects the macroscopic generalized
strain `E_bar`, for example:

```matlab
models = add_material_model(models, ...
    Guan_Gaussian_Hencky([0.1, 1.0], 'SH'), ...
    [0.0, -Inf], [Inf, Inf]);
```

For `Guan_Gaussian_SH`, `m_hat` is the chain-scale Seth-Hill parameter in
`E_hat`. Any generalized-strain parameters for `E_bar` follow after `m_hat`;
when `E_bar` also uses the Seth-Hill family, the parameter vector is
`[mu, m_hat, strain_m]`.

Generalized-strain model:

| Constructor | Parameters | Notes |
| --- | --- | --- |
| `Hill_GenStrain` | `[mu, strain parameters...]` | Hill-type model with `W = mu * E:E`. |

The available macroscopic generalized strain families are shared by
`Guan_Gaussian_GL`, `Guan_Gaussian_Biot`, `Guan_Gaussian_SH2_3`,
`Guan_Gaussian_Hencky`, `Guan_Gaussian_SH`, and `Hill_GenStrain`:

| Family | Constructor argument | Extra parameters |
| --- | --- | --- |
| Seth-Hill | `'SH'` | `[m]` |
| Hencky | `'Hencky'` | `[]` |
| Biot | `'Biot'` | `[]` |
| Curnier-Rakotomanana | `'CR'` | `[m, n]` |
| Curnier-Zysset | `'CZ'` | `[m]` |
| Darijani-Naghdabadi | `'DN'` | `[m, n]` |

Models can be combined by calling `add_material_model` multiple times. Bounds
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
`src/add_exp_data_sets.m`.

## Extending

To add a material model, create a constructor in `material_models/` that
returns a struct with these fields:

- `name`
- `parameters`
- `parameter_names`
- `energy`
- `S`
- `P`
- `set_parameters`

Then add it in `driver.m` with `add_material_model`.

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
- Guan, J., Li, X., Yuan, H., & Liu, J. (2025). *Hyperelastic modeling based on
  generalized Landau invariants and multi-stage calibration.* **Journal of the
  Mechanics and Physics of Solids**, 106338.
- Zhan, L., Wang, S., Qu, S., Steinmann, P., & Xiao, R. (2023). *A new
  micro-macro transition for hyperelastic materials.* **Journal of the
  Mechanics and Physics of Solids**, 171, 105156.

## Author

- Chongran Zhao
- Southern University of Science and Technology, China
- [chongranzhao@outlook.com](mailto:chongranzhao@outlook.com)
