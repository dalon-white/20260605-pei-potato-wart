# PEI Potato Wart Border Sampling Analysis

This repository provides a complete R analysis workflow to design and simulate a zero-acceptance border sampling strategy for detecting **pest A (fungus)** in **commodity B** shipments.

## Project layout

- `R/`: modular analysis code (`config`, `distributions`, `sampling_rules`, `upstream_screening`, `border_sampling`, `simulation`, `metrics`, `plotting`, `report`)
- `config/default.yaml`: default assumptions and tunable parameters
- `scripts/run_cli.R`: command line entry point
- `tests/`: focused tests for sample-size math and simulation behavior
- `outputs/`: generated CSV outputs, figures, and markdown report

## Methods

### Shipment-level zero-acceptance sample size

- Binomial approximation:  
  `n >= log(alpha) / log(1 - p0)`
- Hypergeometric finite-population exact zero-detection probability:  
  `P0 = choose(N-D, n) / choose(N, n)` with `D = ceiling(p0*N)`
- Sensitivity-adjusted approximation with effective sensitivity `Se_eff`:  
  `P0 = (1 - p0 * Se_eff)^n`

### Border policy

- Any positive border test rejects the shipment.
- Border sampling is asymptomatic-only.
- Tiered rule minimizes tests while keeping assurance:
  - Small `N_asym`: test all
  - Medium/large `N_asym`: hypergeometric/binomial floor with sensitivity inflation

### Annual simulation

`simulate_year()` performs Monte Carlo at shipment level, including:
1. Shipment count generation (overdispersed)
2. Shipment mass and unit-size draws → unit count
3. Beta-binomial prevalence/clustering
4. Upstream visual screening conditioning
5. Border sampling and detection outcome
6. Annual aggregation of risk and test burden

## Usage

Install required packages:

```r
install.packages(c("yaml","dplyr","tidyr","purrr","ggplot2","stringr","tibble","readr","scales"))
```

Run end-to-end:

```bash
Rscript scripts/run_cli.R --config config/default.yaml --seed 1234 --outputs_dir outputs
```

Optional scenario override:

```bash
Rscript scripts/run_cli.R --scenarios 0.001,0.01,0.1
```

## Outputs

- `outputs/csv/shipment_level.csv`
- `outputs/csv/annual_summary.csv`
- `outputs/figures/required_n_vs_N.png`
- `outputs/figures/oc_curve.png`
- `outputs/figures/annual_risk_sensitivity.png`
- `outputs/report.md`

## Assumptions and limitations

- Analytic sample-size functions are explainable approximations; finite-population and workflow complexities are captured in simulation.
- Upstream screening performance and asymptomatic fraction are configurable and materially affect annual risk.
- Results should be re-run as new data on prevalence, clustering, test performance, and shipment composition become available.
