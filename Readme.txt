Works： 
The article has been published in Electronics and is open access. Those who are interested can access it at https://doi.org/10.3390/electronics15132822.
Ye xu, Peng Li*, , Wenhui Wang, Youpeng Sun, et.l “Extension Difference Mapping Based PMBM Filter for Non-Ellipsoidal Extended Target Tracking”


Author: Ye xu, Peng Li*, , Wenhui Wang, Youpeng Sun, et.l
Year: 2026

The "scenario" in Line 6 controls the operating scenarios and parameters of the filters.
The "numMC" in Line 20 controls the number of Monte Carlo runs for all filters.


THE Switchin  Lines    36,            122,            198          301        400      are related to the 
                    EDM-PMBM,     GGIW-PMBM  ,  Smooth-pmbm   , KS-PMBM   ,  BGGIW-PMBM and GBePMB   respectively. 
While the value of "if" is 0,   the filter is turned off; when the value is 1, the filter is turned on.


THE Switchin  Lines  574 is control both the BGGIW-PMB and GBePMB.
    The "i_filter" in Line 365, When "i_filter = 1", run the    BGGIW-PMB  filter. 
                                When "i_filter = 2", run the     GBePMB   filter.

*EDM-PMBM || This Work
https://doi.org/10.3390/electronics15132822

*GGIW-PMBM From DOI: 10.23919/FUSION43075.2019.9011181
Xia, Yuxuan, et al. "Extended target Poisson multi-Bernoulli mixture trackers based on sets of trajectories." 2019 22th International Conference on Information Fusion (FUSION). IEEE, 2019.

*KS-PMBM From DOI: 10.1109/JSEN.2024.3521436
Li, Peng, et al. "Poisson multi-bernoulli mixture filter for multiple extended object tracking using kolmogorov–smirnov test." IEEE Sensors Journal 25.4 (2025): 6541-6555.

*BGGIW-PMBM and *GbePMB From https://doi.org/10.1016/j.dsp.2025.105204
Chen, Cheng, Jinlong Yang, and Jianjun Liu. "Adaptive Poisson multi-Bernoulli filter for multiple extended targets with Gamma and Beta estimator." Digital Signal Processing 163 (2025): 105204.
