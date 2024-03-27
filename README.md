# Discrete-time Kalman filter implementation methods
 
This repository contains MATLAB functions for various Kalman filter (KF) implementation methods. They are given in a posteriori form, i.e., no data are assumed to be known at the initial step and, hence, the time update comes first. 

# References
Each code (implementation method) includes the exact reference where the particular algorithm was published. 
If you use these codes in your research, please, cite the corresponding articles mentioned in the codes.  

# Remark
The codes have been presented here for their instructional value only. They have been tested with care but are not guaranteed to be free of error and, hence, they should not be relied on as the sole basis to solve problems. 

# Steps to reproduce
- `Test_KFs` is the script that performs Monte Carlo runs for solving filtering problem by various KF implementations.
- `Test_LLF` is the script for calculating the negative log LF by various filtering algorithms. 
- `Illustrate_XP` is the script that illustrates the obtained estimates and the diagonal entries of the error covariance matrix (over time). You can find its call at the end of the script above, which is commented. Just delete this comment sign.
- `Illustrate_LLF` is the script that illustrates the negative log LF calculated by various filtering algorithms. 
- `Simulate_Measurements` stands for simulating the state-space model and generating the measurements for the filtering methods.

When the state is estimated, the resulted errors should be the same for all implementation methods because they are mathematically equivalent to each other. Their numerical properties differ, but the ill-conditioned test examples are not given here. 

# List of the KF implementation methods
**Conventional algorithms:**
| Function | Description |
| ---: | :--- |
| `Riccati_KF_standard` | Conventional implementation, original method[^1] |
| `Riccati_KF_Joseph` | Conventional Joseph stabilized implementation[^2] |
| `Riccati_KF_Swerling` | Conventional implementation based on Swerling's formula[^3] |
| `Riccati_KF_seq`      | Sequential Kalman Filter (component-wise measurement update)[^4] |


### Square-root algorithms 
Cholesky factorization-based methods:
 -  `Riccati_KF_SRCF_QL`   is the Square-Root Covariance Filter (SRCF) with lower triangular factors by Park & Kailath (1995), <a href="http://doi.org/10.1109/9.384225">DOI</a> 
 -  `Riccati_KF_SRCF_QR`   is the SRCF with upper triangular factorsby Park & Kailath (1995), <a href="http://doi.org/10.1109/9.384225">DOI</a> 
 -  `Riccati_KF_SRCF_QR_seq` is the Sequential SRCF with upper triangular factors by Kulikova (2009), <a href="http://dx.doi.org/10.1134/S0005117909050129">DOI</a>  
 -  `Riccati_KF_eSRCF_QL`  is the Extended SRCF with lower triangular factors by Park & Kailath (1995), <a href="http://doi.org/10.1109/9.384225">DOI</a>  
 -  `Riccati_KF_eSRCF_QR`  is the Extended SRCF with upper triangular factors by Park & Kailath (1995), <a href="http://doi.org/10.1109/9.384225">DOI</a> 

Singular value decomposition (SVD) factorization-based methods:
 -  `Riccati_KF_SVDSR`     is the SVD-vased Filter by L. Wang et.al. (1992), <a href="http://doi.org/10.1109/CDC.1992.371522">DOI</a>
 -  `Riccati_KF_SVD`       is the SVD-based Covariance Filter by Kulikova & Tsyganova (2017), <a href="http://doi.org/10.1049/iet-cta.2016.1282">DOI</a>
 -  `Riccati_KF_SVDe`      is the "economy size" SVD-based Covariance Filter by Kulikova et.al. (2021), <a href="10.1016/j.cam.2019.112487">DOI</a>

[^1]: Kalman, R.E. (1960) A new approach to linear filtering and prediction problems. Journal of basic Engineering. 1960 Mar, 82(1):35-45. <a href="https://doi.org/10.1115/1.3662552">DOI</a>
[^2]: Bucy, R.S. and Joseph, P.D. Filtering for Stochastic Processes, with Applications to Guidance. New York, John Wiley & Sons, 1968.
[^3]: Swerling, P. (1959) First order error propagation in a stagewise differential smoothing procedure for satellite observations, Journal of Astronautical Sciences, V.6, 46--52. 
[^4]: Grewal, M.S. and Andrews, A.P. Kalman filtering: theory and practice using MATLAB. Prentice-Hall, New Jersey, 4th edn., 2015. 


