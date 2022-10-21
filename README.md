# Causal Regularization

Files for reproducing the simulation and experimental sections of **Causal Regularization: A trade-off between in-sample risk
and out-of-sample risk guarantees** (Kania & Wit).

### Install required packages

Execute the following command to install the necessary packages

```
Rscript requirements.R
```

## Simulations and experiments

To reproduce the figures presented in the paper, please execute the
following commands.

### Convergence of out-of-sample risk at square-root rate) (Section 7.1, Figure 2)

```
Rscript ./experiments/finite_sample_bound/finite_sample_bound.R
```

### Bootstrap confidence interval (Section 7.2, Figure 3)

```
Rscript ./experiments/bootstrap_cis/bootstrap_cis.R
```

### Causal regularization vs causal Dantzig (Section 7.3, Figure 5 and Figure 6)

For figure 5, execute:

```
Rscript ./experiments/out_of_sample_risk/out_of_sample_risk.R
```

For figure 6, execute:

```
Rscript ./experiments/cr_vs_cd/cr_vs_cd.R
```

### Fulton fisher market (Section 8.1, Figure 7  and 8)

```
Rscript ./experiments/fish/fish.R
```

### Gene knockout experiments (Section 8.2, Figure 1)

```
Rscript ./experiments/genes/genes.R
```





