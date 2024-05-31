# ShellOpt (version 0.1.0)

ShellOpt is a MATLAB package for the numerical simulation and shape optimization of thin elastic shells using methods from Riemannian geometry and isogeometric analysis.

## User guide

Please consult the README file in the folder `packages` to install the required third-party packages for the software.

The file `./packages/geopdes/inst/space/@sp_vector/op_KL_shells_tp.m` has to be modified as follows:
- Replace `t_coeffs = t_coeff (x{:});` in line 52 with `t_coeffs = t_coeff;`.

Initialize the search path by running the script `init_path.m`.

Runnable MATLAB scripts are indicated by the prefix `script_`.
