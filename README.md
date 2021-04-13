# occbin-n
Generalization of Occbin package to handle an arbitrary number of regime-switching constraints

## Installation
`solve_n_constraints` is the main function, and `mkdatap_anticipated_n` is an auxiliary function. They require [Dynare](https://www.dynare.org/) to be installed and configured, and Dynare's subdirectory "\matlab\occbin\" must be on the MATLAB path. The functions were written using MATLAB 2018b and Dynare 4.6.2.

## Usage
The user must familiarize yourself with Guerrieri and Iacoviello (2015) "[OccBin: A toolkit for solving dynamic models with occasionally binding constraints easily](https://www.sciencedirect.com/science/article/pii/S0304393214001238)", Journal of Monetary Economics.

The function `solve_n_constraints` leaves a significant degree of freedom to the user, and the necessary details are in the function description, but I may suggest the following way of usage that requires minimal effort.

`solve_n_constraints` handles arbitrary binary regime-switching constraints, but for simplicity, consider the following example. Suppose an optimizing agent is subject to inequality constraints `f_i>=0` for `i=1,...,n` and for all contingencies. Here we focus on a specific contingency. Let `lambda_i` denote the normalized Lagrange multiplier for constraint `i`. The Karush---Kuhn---Tucker conditions will involve the complementary slackness conditions `f_i>=0`, `lambda_i>=0`, and `f_i*lambda_i=0` for all `i`. Without loss of generality, suppose that in the steady state, `lambda_i>0`, and thus `f_i=0` for all `i`---all constraints are binding.

To use `solve_n_constraints`, you need at least two Dynare mod files: a reference mod file, where all constraints are in the reference---steady-state---regime, and an additional mod file, where constraint `1` is in the alternative regime and constraints `2,...,n` are in the reference regimes.

Without further loss of generality, suppose `n=2`.

### Reference mod file
Suppose the reference mod file has a name `name.mod`. You need to format it as follows.
1. Before the model declaration, insert a line
```
@#define regimes=[0,0]
```
If `n=3`, the line must be `@#define regimes=[0,0,0]`, and so on.

2. In the model declaration, declare the complementary slackness conditions as follows:
```
// Regime switching
@#if regimes[1]==0
    0=f_1;
@#else
    0=lambda_1;
@#endif
@#if regimes[2]==0
    0=f_2;
@#else
    0=lambda_2;
@#endif
```
If `f_i` is a model expression, it of course must be declared elsewhere in the model declaration.

### Additional mod file
Copy `name.mod` to `name_10.mod`. (If `n=3`, use `name_100.mod`, and so on.) The naming convention is important. Modify `name_10.mod` as follows.
1. Remove everything except sections `var`, `varexo`, `parameters` (can keep only declaration, no need to assign values), `model`, `shocks`, and the line `@#define regimes=[0,0]`.
2. Change `@#define regimes=[0,0]` to `@#define regimes=[1,0]`.

If you follow this approach, `solve_n_constraints` will create mod files for all other regimes automatically. Moreover, you do not need to supply the `rspace` argument, that is, you can set `rspace=[]`.
