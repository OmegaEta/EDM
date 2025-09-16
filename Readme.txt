Works：
PMBM Filter for Non-ellipsoidal Extended Object Tracking using Extension Difference Mapping

Author: Peng Li, Ye xu, Wenhui Wang*, Congzhe You, Wenqi Geng
Year: 2025



The "scenario" in Line 7 controls the operating scenarios and parameters of the filters.


The "numMC" in Line 17 controls the number of Monte Carlo runs for all filters.


THE Switchin  Lines    32,            113,             184    are related to the 
                    EDM-PMBM,     GGIW-PMBM  ,      KS-PMBM  respectively. 
While the value of "if" is 0,   the filter is turned off; when the value is 1, the filter is turned on.


THE Switchin  Lines  274 is control both the BGGIW-PMB and GBePMB.
    The "i_filter" in Line 365, When "i_filter = 1", run the    BGGIW-PMB  filter. 
                                When "i_filter = 2", run the     GBePMB   filter.


