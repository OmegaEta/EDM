Works：
Extension Difference Mapping Based PMBM Filter for Non-Ellipsoidal Extended Target Tracking 

Author: Peng Li, Ye xu, Wenhui Wang*, Congzhe You, Wenqi Geng
Year: 2026

The "scenario" in Line 6 controls the operating scenarios and parameters of the filters.


The "numMC" in Line 20 controls the number of Monte Carlo runs for all filters.


THE Switchin  Lines    36,            122,            198          301        400      are related to the 
                    EDM-PMBM,     GGIW-PMBM  ,  Smooth-pmbm   , KS-PMBM   ,  BGGIW-PMBM and GBePMB   respectively. 
While the value of "if" is 0,   the filter is turned off; when the value is 1, the filter is turned on.


THE Switchin  Lines  574 is control both the BGGIW-PMB and GBePMB.
    The "i_filter" in Line 365, When "i_filter = 1", run the    BGGIW-PMB  filter. 
                                When "i_filter = 2", run the     GBePMB   filter.


