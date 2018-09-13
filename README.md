# MALT: Meta Adaptive Lasso Tree

![alt text](https://upload.wikimedia.org/wikipedia/commons/thumb/7/7d/Malt_en_grain.JPG/1920px-Malt_en_grain.JPG)

Meta-analysis with individual participant data (IPD-MA) has become increasingly popular because of its advantages in facilitating consistent evaluations and adjustment and increasing statistical power. However, properly handling the study heterogeneity is challenging and usually requires novel statistical methods. How to borrow information across different studies in spite of heterogeneity is the core problem in IPD-MA. In our motivating study, the Alliance for Clinical Trials in Oncology Network is formed by the merging of three major cooperative groups created by NCI, Cancer and Leukemia Group B (CALGB), North Central Cancer Trials Group (NCCTG), and the American College of Surgeons Oncology Group (ACOSOG). The 17 Alliance legacy trials ( 8 CALGB trials, 4 NCCTG trials and 5 ACOSOG trials) have various sample sizes, between 33 and 5300. Each trial was conducted for different population, therefore, the survival outcomes (time to recurrence or death) and the distributions of patients over the cancer characteristics can be considerably different across trials. Some trials even focus on one particular level of T or N. In such situation, although it is still reasonable to assume an invariant cancer stating structure across trials,  it is notably difficult to have a consistent and robust evaluation of it. 

We propose to improve cancer staging by Lasso Tree in meta-analysis by using an adaptive group fused Lasso penalization with a Cox proportional hazard model.

The functionality of the penalty can be explained from two aspects: (1) The "adaptive fuse" part: an intrinsic L1 penalty is applied to the difference of neighboring coefficients. Due to the sparsity-enforcing feature of L1 penalty, the neighboring coefficients that are close will coalesce and therefore lead to cancer staging. (2) The "group" part: The same neighboring pairs of coefficients from different trials are grouped together. With such a setup, we utilize the group-wise-selecting property of group Lasso to ensure the same staging structures across all trials. In our modeling, the auxiliary heterogeneity of trials is naturally absorbed into the semi-parametric baseline function, while the only constant information across trials, the staging structure, is extracted by the grouping of the coefficients.  

## Installation

library(devtools)

install_github("WangTJ/malt")
