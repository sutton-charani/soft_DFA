README
================
Nicolas Sutton-Charani & Francis Faux
2024-02-16

#### Data definition

``` r
x <- Nile
```

``` r
devtools::source_url(paste0("https://raw.githubusercontent.com/sutton-charani/",
                            "soft_DFA/main/code/my_lib_soft_dfa.R"),
                     sha1="b8030baadb314f3747d2bc0de3e150c8e74c9495")
```

#### Precise DFA

``` r
dfa_complexity(x)
```

    ## [1] 1.04962

#### Soft DFA

``` r
soft_dfa(x)
```

    ##   slope_min slope_max       mass
    ## 1 0.8128825  1.427578 0.22465438
    ## 2 0.6712844  1.705121 0.33698157
    ## 3 0.2772593  1.961915 0.42684332
    ## 4      -Inf       Inf 0.02564103
