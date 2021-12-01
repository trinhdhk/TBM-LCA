# Introduction

**Tuberculous meningitis (TBM)** is the most severe form of
tuberculosis. Diagnosing TBM is notoriously challenging, with
microbiological confirmation requiring identification of *Mycobacterium
tuberculosis* in paucibacillary cerebrospinal fluid (CSF). In addition
to widely used confirmatory methods of CSF testing for *M. tuberculosis*
such as Ziehl-Neelsen (ZN) staining, GeneXpert MTB/RIF (Xpert), and
mycobacterial culture, additional parameters may increase the likelihood
of a diagnosis of TBM. Such parameters are illustrated in the uniform
case definition for TBM (Marais et al. 2010) where, in the absence of
positive microbiological tests, an increased certainty of TBM is
assigned in the presence of particular clinical, CSF, and imaging
findings, or with evidence of non-neurological *M. tuberculosis*.

The uniform case definition, albeit widely adopted in clinical practice,
does have its own disadvantages, which mostly born from its expertise
base. The categorised representation of TBM risk into *Definitely no*t,
*Possible, Probable,* and *Definite* rather than actual probability can
be hard to interpret, especially when used in estimating diagnostic
tests’ sensitivities and specificities and treatment’s response (1).
Secondly, scores of risk factors are consensually chosen not
statistically built could lead to overestimation and underestimation of
the risks (2). Apart from the extreme Definite and Definitely not, which
are micro-biologically confirmed, clinical diagnoses outputted from this
scoring system is highly uncertain, leave a lot of room to the
clinicians’ decision (3). Lastly, the fact that it depends on “slow”
laboratorial assays (culturing) likely delays the diagnosis of patients,
which ocassionally might not be fruitful.

## Analysis objectives

The main aims of this analysis are to:

1.  Estimate the latent performances, namely sensitivities and
    specificities, or current TBM confirmation tests, taking into
    account the uncertainty in current diagnosis,
2.  Re-adjust the current scoring system basing on statistical model
    which output an estimation of probability of having TBM, on the
    individual level

As secondary objectives, this analysis also aim to:

1.  Build a simplified scoring system which only needs minimal (ideally
    clinical-only) information but has the capacity to *approximate* the
    full system’s output
2.  Estimate a latent representation of patients’ bacillary burden given
    that they get TBM which may impact the tests results and patient’s
    prognosis.

## Prior knowledges in the field

Estimated sensitivity and specificity of confirmation tests based on the
TBM case definition (Marais et al. 2010) are listed in table
@ref(tab:confirmation-test-tbl) (Nhu et al. 2013).

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

# Methods

## Population to analyse

Data used for this analysis was extracted from an ongoing observational
study conducted at the Hospital for Tropical Diseases (HTD) - Vietnam
from the end of August 2017 to the end of January 2021. This is a large
centre in southern Vietnam, provides secondary and tertiary treatment
for a wide range of tropical infections (Thwaites et al. 2002). In the
main study, enrolled population includes patients at least 16 years old
with suspected neurological infection, admitted to Viet-Anh Ward, and
underwent lumbar puncture at baseline as a routine diagnostic procedure.
They were monitored and clinical data were prospectively collected
during treatment time, until either discharge or death. Exclusion
criteria include contra-indication of lumbar puncture and the invalidity
of informed consent. A subset of the same data was randomised and
analysed in a past study in which performance of Xpert MTB/RIF Ultra and
Xpert MTB/RIF was compared (Donovan et al. 2020).

## Test procedures and data collection

At baseline, after informed consents were provided, demographic and
history information were curated. Patients also underwent clinical
examination and laboratorial investigations according to main study
protocol, including: blood test, sputum, and lumbar puncture, unless
contra-indicated. HIV-related information was asked and optionally
tested unless the patients were to enrol in a subsequent Tuberculous
meningitis randomised controlled trial also conducted at HTD. Patients
with either an known HIV infection or a positive subsequent HIV test
would be considered **HIV positive**.

During lumbar puncture, at least 3mL (ideally 6mL) of cerebrospinal
fluid (CSF) were taken; if lower amount was collected, the tests would
still be done, with collected volume noted. Tests might include standard
haematological and biochemical bio-markers, confirmation markers for
differential diagnoses, gene expression and profiling, if applicable.

TBM confirmation tests used in this analysis are **Ziehl–Neelsen stained
smear (ZN Smear), Culturing in Mycobacteria Growth Indicator Tube**
(**MGIT**)**,** and **GeneXpert MTB/RIF** or **Xpert Ultra MTB/RIF
(Xpert)**; as the last two tests’ performances are comparable (Donovan
et al. 2020), we considered them to be the same. The confirmation tests
were done if TBM were suspected clinically, either by the TBM definition
scoring system (Marais et al. 2010) or by clinicians’ judgment and no
other diagnosis was confirmed. In order to perform these tests, CSF
samples were centrifuged at 3000g for 15 minutes (Donovan et al. 2020).
<!--# Joe, can you write some introduction about the tests and procedure, as in XpertUltra paper -->

All patients underwent appropriate treatment regimen according to
national and local guidelines, depending on the diagnosis, without any
interference ignited by the study. During the treatment, patients might
continue to have additional standard-of-care lumbar puncture for
diagnosis or follow-up. At the time of discharge or death, all patients
received a final diagnosis. If at least one ZN Smear, MGIT, or Xpert,
was positive at anything time during the follow-up, the patient would be
considered **confirmed TBM**; otherwise, **suspected TBM**, unless the
patient recovered without anti-tuberculosis chemotherapy - where they
could be reassigned to **not TBM**.

Although tests could be done several times as mentioned, in this
analysis, in order to robustly evaluate the sensitivities of
confirmation tests, only the first samples at baseline were curated.
Results of haematological, biochemical, ZN Smear, MGIT, and Xpert, all
must have come from one first CSF sample, together with which a quick
blood glucose should also be taken, as listed in table
@ref(tab:test-bl-tbl). In general, the choice of predictors are based on
the TBM definition scoring system (Marais et al. 2010). CSF Oeosinophil
count was additionally included as it is a strong bio-marker for
oeosinophilic meningitis, a condition usually caused by the parasites;
while CSF erythrocyte count was added as a marker for a traumatic
procedure. We also made an adjustment for Past TB contact, in which we
hypothesise that by changing the question into “Past *noticeable*
contact with TB patients within the past recent year”, we would be able
to imply all “Unknown” answer as a “No”. In the main study, no brain
imaging was taken, hence not included.

## Statistical Model

All data preparation, cleaning, and processing were performed on
statistical package , version 4.1.0 (R Core Team 2021). The model was
developed on the probabilistic language via the interface , version 2.27
(Stan Development Team 2021a). Plotting was done using package (Gabry et
al. 2019), (Defazio and Campbell 2020), and (Yan 2021). Some other
packages used include: (Chang 2020), (Harrell Jr 2021).

In the model, three aforementioned TBM confirmation tests (ZN Smear,
MGIT, and Xpert) were used as manifest variables. Linear predictors for
the latent class prevalence were **history and demographic information**
(HIV status, Age, Past noticeable TB contact), **clinical signs and
symptoms** (Days from onset to admission, Systemic symptoms suggestive
of tuberculosis, Focal neurological deficit, Cranial nerve palsy,
Glasgow Coma Score - GCS), **Imaging chest data** (Pulmonary TB, Miliary
TB), **CSF criteria (**lymphocyte count, neutrophil count, oeosinophil
count, erythrocyte count, glucose - with corresponding blood glucose,
protein, lactate), and **Cryptococcus Antigen/Indian Ink** in
combination**.**

### Data pre-processing

Most continuous variables were transformed to logarithmic scale and
subsequently centred. In all base models, no scaling was performed as
they could potentially imputed.
<!--# ; CSF bio-markers were standardised only if they were to be fed into a Probabilistic Principal Component Analysis (PPCA) for Dimensionality Reduction-->.
Glasgow Coma Score and its components were not transformed, rather we
translated them to **Loss of GCS** (LoGCS); so that a *G**C**S* = 15
would be equivalent to *L**o**G**C**S* = 0, while *G**C**S* = 3 would be
translated to *L**o**G**C**S* = 12.
<!--# Exceptionally in the imputation step, due to technical difficulties from having different upper bounds, we scaled compartments of $LoGCS$ to unit vectors, so that they all shared the same lower bound at 0 and upper bound at 1.-->

Binary variables were encoded into 0 and 1 and not centred.

### Latent class regression model

We created two-level hierarchical models each of which combines:

-   **Prevalence model**: A logistic regression model trying to estimate
    the prevalence of TBM amongst the study population

-   **Latent class analysis**: estimating the probabilities of having
    positive results from three tests in each class. Similar to previous
    applications (Qu, Tan, and Kutner 1996; Hadgu and Qu 2002;
    Schumacher et al. 2016), we also corrected for individual bacillary
    burden and procedural variance between samples. Latent bacillary
    burden is regressed by individual Gaussian random variables which
    captures noisy fluctuation of test results, and fixed effects coming
    from covariates. We hypothesise that even if two patients were in
    the same TBM positive class, one who had lower bacillary burden
    would be less likely to be tested positive. This lifted the local
    independence assumptions of vanilla Latent Class Analysis
    <!--# citation needed -->.

<!-- -->

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

The inclusion of covariates was performed according to prior knowledge
of potential predictors that share association with the infection risk
and test sensitivity, summarised in table @ref(tab:predictor-tab)). Due
to the lack of understanding about the mechanism of the disease and the
main interest of the analysis, we did not investigate the causal
relationship between different covariates.

<img src="main_files/figure-markdown_strict/skeleton-model-1.png" alt="Model skeleton. Bacillary Burden is only available in model 2+" width="100%" height="5cm" />
<p class="caption">
Model skeleton. Bacillary Burden is only available in model 2+
</p>

The skeleton for all models is shown in figure @ref(fig:skeleton-model).
We use a stepwise approach where we incrementally added up more
flexibility and lifted more constraints. We also include several
extensions to explore many possibilities than can improve performance.
In the main analysis, only CSF neutrophil count had quadratic effect as
suggested by the TBM case definition (Marais et al. 2010). A summary of
all architectures an extensions are shown in table
@ref(tab:model-archs), whereas technical formulations are detailed in
appendix @ref(appendix-details).

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

### Missing data handling

#### Manifest variables

By design, most patients with very high chance of and/or evidently
diagnosed with different diseases were not tested with TBM confirmation
assays (namely ZN Smear, MGIT, and Xpert), unless there were excessive
amount of CSF samples. However, as the tests’ sensitivities were all
firmly believed to be almost perfect (Nhu et al. 2013), we assumed that
patients who had no TBM confirmation tests are all negative. Apart from
one premature death, we have yet to find any patients left un-tested and
un-diagnosed.

#### Predictors

In this analysis, missing predictors’ values were assumed to be Missing
At Random (MAR) and imputed within sampling programmes, together with
the main model. Composite predictor variables, such as *TB-suggested
symptoms* and *Glasgow Coma Score (GCS)* were imputed by compartments;
while potentially correlated variables were grouped and imputed
together. Due to Stan not supporting Multivariate Logistic Regression,
in favour of method consistency, all binary predictors were imputed
using (Multivariate) Probit models. Continuous variables are imputed
using Multivariate Linear
Regression<!--# , with one exception of PPCA-extended models where the imputation of CSF bio-markers was handled by the PPCA directly -->.
HIV status was included in most imputation model as predictors.

As HIV tests were not mandatory in the main study, their chance of
missingness were mostly dependent on whether or not they were to enrol
in a TBM study. Hence, it is safely to assume that HIV status are
Missing at Random (MAR). Accordingly, we imputed HIV using probit
regression, corrected for Blood Lymphocyte and Neutrophil counts.

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

We summaries our rationales and corresponding handling strategies in
table @ref(tab:missing-handling). In the case where missing values were
imputed, figure @ref(fig:impute-model) depicts how those were sampled,
together with their potential hyper-predictors. Note that due to
Hamiltonian Monte Carlo (HMC)’s limitation, imputed LoGCS were treated
as-is in sampling process, but were rounded when estimating model
performance.

<img src="main_files/figure-markdown_strict/impute-model-1.png" alt="Imputation strategy for predictors. Variables in rectangular solid boxes were used in the model, in oval dashed were either compartments of composite ones or only contributed in the imputation model and were not included in the main model. Clustered covariables were imputed together in a multivariate regression. Arrows demonstrate a predictor-response correlation." width="100%" height="100%" />
<p class="caption">
Imputation strategy for predictors. Variables in rectangular solid boxes
were used in the model, in oval dashed were either compartments of
composite ones or only contributed in the imputation model and were not
included in the main model. Clustered covariables were imputed together
in a multivariate regression. Arrows demonstrate a predictor-response
correlation.
</p>

### Prior choices

Following Gelman’s recommendations (Stan Development Team n.d.), we
chose *N**o**r**m**a**l*(0,2.5) as prior distributions for all
intercepts and coefficients in the imputation model, except for LoGCS
(E, V, and M) as their supports were constrained to \[0,1\]; in the
latter case, *U**n**i**f**o**r**m*(0,1) was chosen instead for the
means, and *N**o**r**m**a**l*(0,.5) for the standard deviation.

Cholesky decomposition of the covariance matrix was sampled from a *LKJ
Correlation Cholesky* prior with *s**c**a**l**e* = 4 (Stan Development
Team 2021b):

*L* ∼ *L**K**J**C**o**r**r**C**h**o**l**e**s**k**y*(4);

In every linear sub-model of main one, we used *t*<sub>4</sub>(0,5) for
the intercept on which we impose our weak expectation that its absolute
value cannot be higher than 10 (Boonstra, Barbaro, and Sen 2019). For
the covariates, we considered several sets of prior representing a
spectrum of penalties bestowed upon the model (van Erp, Oberski, and
Mulder 2019):

-   Weakly informative prior: *t*<sub>4</sub>(0,*s*)

-   Ridge-equivalent prior: *N**o**r**m**a**l*(0,*s*)

-   LASSO-equivalent prior:
    *D**o**u**b**l**e**E**x**p**o**n**e**n**t**i**a**l*(0,*s*)

    *s*.*t* : *s* ∼ *N**o**r**m**a**l*(0,2.5) & *s* &gt; 0

Coefficient for known strongly positive risk factors (marked **++** in
table @ref(tab:predictor-tab) were specially imposed a positive Half
Normal distribution. Individual random effects representing unmeasured
bacillary burden was sampled from a *N**o**r**m**a**l*(0,1)
distribution.

For manifest variables, we used highly informative priors for
*specificity (Spc)* basing on previous study (Nhu et al. 2013) and
weakly informative priors for *sensitivity (Sen)* on the logit scale.
These choices were visualised on figure @ref(fig:mv-priors).

$$
\\begin{aligned}
1-Spc\_{Xpert} &\\sim Logistic(logit(0.005), .7)\\\\
1-Spc\_{MGIT}  &\\sim Logistic(logit(0.001), .3)\\\\
1-Spc\_{ZN\\ Smear} &\\sim Logistic(logit(0.001), .3)\\\\
Sen\_{Xpert,\\ MGIT,\\ ZN\\ Smear} &\\sim Logistic(0,.5)
(\\#eq:priors-response)
\\end{aligned}
$$

<img src="main_files/figure-markdown_strict/mv-priors-1.png" alt="Density plots for different priors. A: Logistic(logit(.005),.7), B: Logistic(logit(.001),.3), C: Logistic(0,.5). 1: Logit scale, 2: Linear scale." width="100%" height="10cm" />
<p class="caption">
Density plots for different priors. A: Logistic(logit(.005),.7), B:
Logistic(logit(.001),.3), C: Logistic(0,.5). 1: Logit scale, 2: Linear
scale.
</p>

### Model performances

Models performances were defined and compared by three metrics:

-   Expected log point-wise predictive density (elpd) of hold-out
    observations (Vehtari, Gelman, and Gabry 2016): Models with higher
    elpd arw supposed to have more predictive values for untrained
    observations.

-   Visualised model calibration between estimated predicted and
    observed probabilities of positive confirmation tests (Harrell
    Jr 2021), using a non-parametric *loess* fit. We also used the
    diagnosis at discharge as a pseudo-gold standard to visualise the
    calibration of TBM prevalence. A model with good calibration would
    correctly estimate the observed values.

-   Receiver Operating Characteristic (ROC) curve and corresponding
    Areas Under the Curve (AUC) . Confidence interval for AUC was
    estimated by a 2000-time bootstrapping process ***&lt;I might change
    this to a fully bayesian estimation. shall I? >***. A model with
    good discriminative value would be better to distinguish between two
    class.

-   Additionally, we also visualised class-wise predicted probability
    density plots to visualise how much separable the classes are based
    on the models. This demonstrate how predicted probabilities
    distributed between two classes.

All metrics are based on 5 repetitions of 20-fold cross validated
datasets (i.e. 100 fits, as suggested by Harrell (Harrell Jr 2021)).
Accordingly, we selected three best-performed models and re-estimated
their parameters using the full dataset.

We also test for local independence
assumption<!--# clarification needed -->.

### Exploratory and sensitivity analysis

#### Complete-case analysis and imputation under MNAR assumption

To check the level of impact from our imputation method, we did a
complete-case analysis and a “missing-as-a-category” analysis in which
we considered missing values as a level for binary variables. Suspected
MNAR variables as listed in table @ref(tab:missing-handling) were also
tested for MNAR where we randomly allocated values based on experts’
opinions.

Under MNAR assumption, we did a pattern-mixture method(Mason et al.
2017; White et al. 2007), where we inquired prediction offsets *δ*. *δ*
represents the difference between unobserved part and observed part of
each variables, after correction for all hyper-predictors listed in
@ref(fig:impute-model). The offsets were collected from interviews with
experts working at the Viet Anh Ward at HTD, where they were supposed to
provide an estimation and 95% confidence intervals (95% CI), based on
which *δ*s were then sampled from.

*δ* ∼ *N**o**r**m**a**l*(*μ*<sub>*e**x**p**e**r**t*</sub>,*s**d*<sub>*e**x**p**e**r**t*</sub>)
*where* *μ*<sub>*e**x**p**e**r**t*</sub> *and*
*s**d*<sub>*x**p**e**r**t*</sub> *are estimation and* $\\frac{1}{2}$ *of
95% CI provided by interviewed experts*.

The main model underwent a relief of strict priors for test
specificities in @ref(eq:priors-response), the lifted model used the
same prior for all three tests which cover specificity from at least
90%:

1 − *S**p**c*<sub>*l**i**f**t**e**d*</sub> ∼ *L**o**g**i**s**t**i**c*(*μ*,.7)
*where* *μ* = *l**o**g**i**t*(0.001) for ZN Smear and MGIT,
*μ* = *l**o**g**i**t*(0.005) for Xpert.

Lastly, as recent studies suggested a sup-optimal specificity of Xpert
test on CSF samples(Nhu et al. 2013; Chen et al. 2020), our assumptions
made in table @ref(tab:missing-handling) might not completely valid. To
tackle this, we considered a MAR scenario, where observation chance of
confirmation tests depend on the unknown TBM status and locally
independent to the value of confirmation tests. The observation status
was then included in the model as a separated manifest variables **(to
Ronald: should I left this in the sensitivity or include this in the
main analysis, for Xpert only or for all three?).** The validity of this
method was depicted in a simulation study in
@ref(appendix-simulation-study).

In this analysis, we expected the chance in which at least one
confirmation test was done are 95% for TBM-positive patients, and 50%
for TBM-negative, hence led to two conservative priors:

$$
\\begin{aligned}
obs\_{TBM-} &\\sim Logistic(logit(.5), .7) \\\\
obs\_{TBM+} &\\sim Logistic(logit(.95), 1)
\\end{aligned}
$$

<img src="main_files/figure-markdown_strict/liftes-priors-1.png" alt="Lifted prior for ZN Smear and MGIT specificities, compared to the old one (dark grey). 1: logit scale, 2: linear scale" width="100%" height="8cm" />
<p class="caption">
Lifted prior for ZN Smear and MGIT specificities, compared to the old
one (dark grey). 1: logit scale, 2: linear scale
</p>

#### Non-linearity and Dimensionality Reduction

As part of the exploratory analysis, we considered two models:

-   Non-linearity for all CSF bio-markers: We added quadratic effects
    for all CSF bio-markers in the prevalence model. LASSO-based
    variable selection was implemented for this analysis.

-   Non-linearity for Glasgow coma scores: As GCS is not a continuous
    variable but rather an ordinal one de facto, it is possible that
    there is non-linear correlation between GCS and TBM risk. In this
    analysis, we employed a quadratic effect for GCS to capture this
    potential of non-linearity.

-   Probabilistic Principal Component Analysis (PPCA): Instead of
    performing Variable selection, we performed an implementation of
    Probabilistic Principal Component Analysis for dimensionality
    reduction, especially amongst potentially collinear predictors. In
    this analysis, we only implemented this for collinearity-prone CSF
    bio-markers.

-   Test accuracy with respect to CSF volume: By including the volume of
    sample collected, we can further investigate the effect size of CSF
    volumes on the sensitivity of each confirmation test.

### Simplified approximation of TBM risk

As the full model requires a plethora of predictors and measurements, it
might not be pragmatic in some limited contexts. We hence developed a
simplified version of the prevalence sub-model which exclude all
laboratorial features. The aim of this exploratory analysis is to
provide a decent approximation of TBM risk yet needs only a minimal
amount of information.

In this analysis, we took out the posterior probabilities of TBM from
the best-performed model, on the logit scale, fed them into subsequent
model where the number of predictors were reduced. The error term
follows Logistic distribution with *s**c**a**l**e* = 1.

*z*<sub>*a**p**p**r**o**x*</sub> ∼ *L**o**g**i**s**t**i**c*(*z*,1);

*s.t.* *z* = *l**o**g**i**t*(*P*<sub>*T**B**M*</sub>) and
*z*<sub>*a**p**p**r**o**x*</sub> is an approximation of *z*.

# Results

Venn diagram for confirmation test results are demonstrated in
@ref(fig:venn-test)

<img src="main_files/figure-markdown_strict/venn-test-1.png" alt="Venn diagram for ZN Smear, MGIT, and Xpert" width="60%" />
<p class="caption">
Venn diagram for ZN Smear, MGIT, and Xpert
</p>

Blahblahblah

<!---BLOCK_LANDSCAPE_START--->

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

<!---BLOCK_LANDSCAPE_STOP--->
# Discussion

# Conclusion

# References

Boonstra, Philip S., Ryan P. Barbaro, and Ananda Sen. 2019. “Default
Priors for the Intercept Parameter in Logistic Regressions.”
*Computational Statistics & Data Analysis* 133 (May): 245–56.
<https://doi.org/10.1016/j.csda.2018.10.014>.

Chang, Winston. 2020. *R6: Encapsulated Classes with Reference
Semantics*. <https://CRAN.R-project.org/package=R6>.

Chen, Yuan-Zhi, Li-Chang Sun, Yao-Hong Wen, Zhong-Wei Li, Shu-Jin Fan,
Hong-Kun Tan, Min Qiu, et al. 2020. “Pooled Analysis of the Xpert
MTB/RIF Assay for Diagnosing Tuberculous Meningitis.” *Bioscience
Reports* 40 (1). <https://doi.org/10.1042/bsr20191312>.

Defazio, Aaron, and Huw Campbell. 2020. *Classifierplots: Generates a
Visualization of Classifier Performance as a Grid of Diagnostic Plots*.
<https://CRAN.R-project.org/package=classifierplots>.

Donovan, Joseph, Do Dang Anh Thu, Nguyen Hoan Phu, Vu Thi Mong Dung,
Tran Phu Quang, Ho Dang Trung Nghia, Pham Kieu Nguyet Oanh, et al. 2020.
“Xpert MTB/RIF Ultra Versus Xpert MTB/RIF for the Diagnosis of
Tuberculous Meningitis: A Prospective, Randomised, Diagnostic Accuracy
Study.” *The Lancet Infectious Diseases* 20 (3): 299–307.
<https://doi.org/10.1016/s1473-3099(19)30649-8>.

Gabry, Jonah, Daniel Simpson, Aki Vehtari, Michael Betancourt, and
Andrew Gelman. 2019. “Visualization in Bayesian Workflow.” *J. R. Stat.
Soc. A* 182: 389–402. <https://doi.org/10.1111/rssa.12378>.

Hadgu, A., and Y. Qu. 2002. “A Biomedical Application of Latent Class
Models with Random Effects.” *Journal of the Royal Statistical Society:
Series C (Applied Statistics)* 47 (4): 603–16.
<https://doi.org/10.1111/1467-9876.00131>.

Harrell Jr, Frank E. 2021. *Rms: Regression Modeling Strategies*.
<https://CRAN.R-project.org/package=rms>.

Marais, Suzaan, Guy Thwaites, Johan F Schoeman, M Estée Török, Usha K
Misra, Kameshwar Prasad, Peter R Donald, Robert J Wilkinson, and Ben J
Marais. 2010. “Tuberculous Meningitis: A Uniform Case Definition for Use
in Clinical Research.” *The Lancet Infectious Diseases* 10 (11): 803–12.
<https://doi.org/10.1016/s1473-3099(10)70138-9>.

Mason, Alexina J, Manuel Gomes, Richard Grieve, Pinar Ulug, Janet T
Powell, and James Carpenter. 2017. “Development of a Practical Approach
to Expert Elicitation for Randomised Controlled Trials with Missing
Health Outcomes: Application to the IMPROVE Trial.” *Clinical Trials* 14
(4): 357–67. <https://doi.org/10.1177/1740774517711442>.

Nhu, N. T. Q., D. Heemskerk, D. D. A. Thu, T. T. H. Chau, N. T. H. Mai,
H. D. T. Nghia, P. P. Loc, et al. 2013. “Evaluation of GeneXpert MTB/RIF
for Diagnosis of Tuberculous Meningitis.” *Journal of Clinical
Microbiology* 52 (1): 226–33. <https://doi.org/10.1128/jcm.01834-13>.

Qu, Yinsheng, Ming Tan, and Michael H. Kutner. 1996. “Random Effects
Models in Latent Class Analysis for Evaluating Accuracy of Diagnostic
Tests.” *Biometrics* 52 (3): 797. <https://doi.org/10.2307/2533043>.

R Core Team. 2021. *R: A Language and Environment for Statistical
Computing*. Vienna, Austria: R Foundation for Statistical Computing.
<https://www.R-project.org/>.

Schumacher, Samuel G., Maarten van Smeden, Nandini Dendukuri, Lawrence
Joseph, Mark P. Nicol, Madhukar Pai, and Heather J. Zar. 2016.
“Diagnostic Test Accuracy in Childhood Pulmonary Tuberculosis: A
Bayesian Latent Class Analysis.” *American Journal of Epidemiology* 184
(9): 690–700. <https://doi.org/10.1093/aje/kww094>.

Stan Development Team. 2021a. “RStan: The R Interface to Stan.”
<https://mc-stan.org/>.

———. 2021b. *Stan Modeling Language Users Guide and Reference Manual*.
<https://mc-stan.org/>.

———. n.d. “<span class="nocase">Prior Choice Recommendations
stan-dev/stan Wiki GitHub</span>.” Accessed August 12, 2021.
<https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations>.

Thwaites, GE, TTH Chau, K Stepniewska, NH Phu, LV Chuong, DX Sinh, NJ
White, CM Parry, and JJ Farrar. 2002. “Diagnosis of Adult Tuberculous
Meningitis by Use of Clinical and Laboratory Features.” *The Lancet* 360
(9342): 1287–92. <https://doi.org/10.1016/s0140-6736(02)11318-3>.

van Erp, Sara, Daniel L. Oberski, and Joris Mulder. 2019. “Shrinkage
Priors for Bayesian Penalized Regression.” *Journal of Mathematical
Psychology* 89 (April): 31–50.
<https://doi.org/10.1016/j.jmp.2018.12.004>.

Vehtari, Aki, Andrew Gelman, and Jonah Gabry. 2016. “Practical Bayesian
Model Evaluation Using Leave-One-Out Cross-Validation and WAIC.”
*Statistics and Computing* 27 (5): 1413–32.
<https://doi.org/10.1007/s11222-016-9696-4>.

White, Ian R, James Carpenter, Stephen Evans, and Sara Schroter. 2007.
“Eliciting and Using Expert Opinions about Dropout Bias in Randomized
Controlled Trials.” *Clinical Trials* 4 (2): 125–39.
<https://doi.org/10.1177/1740774507077849>.

Yan, Linlin. 2021. *Ggvenn: Draw Venn Diagram by ’Ggplot2’*.
<https://CRAN.R-project.org/package=ggvenn>.

# (APPENDIX) Appendix

# Model Formulation Details

To be written

# Simulation study

The purpose of this study is to estimate the perfomance an LCA model in
the scenario where manifest variables were partially missing.

In this study, we created *X*, *X*<sub>2</sub>, *X*<sub>3</sub> as three
covariates for latent class *C*. For simplicity,
*Y*<sub>1</sub>, *Y*<sub>2</sub>, and *Y*<sub>3</sub> were three
manifest variables whose data generation processes are locally
independent w.r.t. *C*. *Y*<sub>3</sub>, however were not completely
observed; its observation rate (represented by the binary variable
obs\_rate) are different for *C* = 0 and *C* = 1.

The data generation process was done in version 4.1 (R Core Team 2021)
as below:

    # Misc function
    generate_Y = \(C, probs){
      Cs = sort(unique(C))
      sapply(C, \(c) rbinom(1,1, probs[Cs==c]))
    }

    # Sample size
    N        = 1000 

    # Create data with 3 predictors, X, X2, and X3
    X        =  rnorm(N, 0, 1 )
    X2       =  rnorm(N, 0, 1 )
    X3       = rbinom(N, 1, .3)

    # Create latent class
    probs    = plogis(3*X+X2+5*X3-1)
    C        = sapply(probs, \(p) rbinom(1,1,p))

    # Manifest variables
    Y1       = generate_Y(C, c(.1 ,.3))
    Y2       = generate_Y(C, c(.01,.5))
    Y3       = generate_Y(C, c(.05,.8))

    # Simulate class-aware missing data for Y3. If obs3==1, Y3 is observed
    obs_rate = c(.1, .95)
    obs3     = generate_Y(C, obs_rate)

The likelihood function of the model shall be
$$
P(Y\_1,Y\_2,Y\_3,obs\_3|X,X\_2,X\_3) = \\prod\_{n=1}^N \\sum\_{c=0}^1 (P\_n(Y\_1,Y\_2,Y\_3,obs\_3|C=1)P(C=c|X,X\_2,X\_3)
$$
Consider an individual *n*: If *o**b**s*3 = 1:
$$
\\begin{aligned}
P(y\_1=Y\_1^{(n)},y\_2=Y\_2^{(n)},y\_3=Y\_3^{(n)},obs\_3=1|X\_n,X\_2^{(n)},X\_3^{(n)} &= \\sum\_{c=0}^1 \\prod\_{i=1}^3 P(y\_i=Y\_i^{(n)}|obs\_3=1,C=c)P(obs\_3=1)P(C=c)\\\\
&= \\sum\_{c=0}^1 \\prod\_{i=1}^3 P(y\_i=Y\_i^{(n)}|C=c)P(obs\_3=1)P(C=c)\\\\
\\end{aligned}
$$
(as
*Y*<sub>1</sub>, *Y*<sub>2</sub>, *Y*<sub>3</sub>, *o**b**s*<sub>3</sub>
are conditionally independent on *C*).

Similarly, for those whose *o**b**s*<sub>3</sub> = 0:
$$
\\begin{aligned}
P(y\_1=Y\_1^{(n)},y\_2=Y\_2^{(n)},obs\_3=0|X^{(n)},X\_2^{(n)},X\_3^{(n)} &= \\sum\_{c=0}^1 \\prod\_{i=1}^2 P(y\_i=Y\_i^{(n)}|obs\_3=0,C=c)P(obs\_3=0)P(C=c)\\\\
&= \\sum\_{c=0}^1 \\prod\_{i=1}^2 P(y\_i=Y\_i^{(n)}|C=c)P(obs\_3=0)P(C=c)\\\\
\\end{aligned}
$$

We considered five sets of prior distribution for
*X*, *X*<sub>2</sub>, *X*<sub>3</sub>, *Y*<sub>1</sub>, *Y*<sub>2</sub>, *Y*<sub>3</sub>,
and *o**b**s*<sub>3</sub> corresponding for the level of prior
knowledge. The formulation of priors follow the syntax in
@ref(eq:priors-response).

1.  Weak prior:

$$
  \\begin{aligned}
  (1-Spc\_{X,X\_2,X\_3}) &\\sim Logistic(0,.7) \\\\
  Sen\_{X,X\_2,X\_3} &\\sim Logistic(0,.7) \\\\
  obs\\\_rate|{C=0} &\\sim Logistic(0,.7) \\\\
  obs\\\_rate|{C=1} &\\sim Logistic(0,.7)
  \\end{aligned}
  $$

1.  Good but conservative knowledge for observation rate, weak prior for
    manifest variables:

$$
  \\begin{aligned}
  (1-Spc\_{X,X\_2,X\_3}) &\\sim Logistic(0,.7) \\\\
  Sen\_{X,X\_2,X\_3} &\\sim Logistic(0,.7) \\\\
  obs\\\_rate|{C=0} &\\sim Logistic(logit(.1),.7) \\\\
  obs\\\_rate|{C=1} &\\sim Logistic(0,.7)
  \\end{aligned}
  $$

1.  Good but conservative knowledge for observation rate and manifest
    variables:  
    $$
      \\begin{aligned}
      (1-Spc\_{X}) &\\sim Logistic(logit(.1),.7) \\\\
      (1-Spc\_{X\_2}) &\\sim Logistic(logit(.01),.7) \\\\
      (1-Spc\_{X\_3}) &\\sim Logistic(logit(.05),.7) \\\\
      Sen\_{X,X\_2,X\_3} &\\sim Logistic(0,.7) \\\\
      obs\\\_rate|{C=0} &\\sim Logistic(logit(.1),.7) \\\\
      obs\\\_rate|{C=1} &\\sim Logistic(0,.7)
      \\end{aligned}
      $$

2.  Very strong knowledge for observation rate, weak prior for manifest
    variables:  
    $$
      \\begin{aligned}
      (1-Spc\_{X,X\_2,X\_3}) &\\sim Logistic(0,.7) \\\\
      Sen\_{X,X\_2,X\_3} &\\sim Logistic(0,.7) \\\\
      obs\\\_rate|{C=0} &\\sim Logistic(logit(.1),.25) \\\\
      obs\\\_rate|{C=1} &\\sim Logistic(0,.7)
      \\end{aligned}
      $$

3.  Very strong knowledge for observation rate and manifest variables:  
    $$
      \\begin{aligned}
      (1-Spc\_{X}) &\\sim Logistic(logit(.1),.2) \\\\
      (1-Spc\_{X\_2}) &\\sim Logistic(logit(.01),.2) \\\\
      (1-Spc\_{X\_3}) &\\sim Logistic(logit(.05),.2) \\\\
      Sen\_{X,X\_2,X\_3} &\\sim Logistic(0,.7) \\\\
      obs\\\_rate|{C=0} &\\sim Logistic(logit(.1),.25) \\\\
      obs\\\_rate|{C=1} &\\sim Logistic(0,.7)
      \\end{aligned}
      $$

All models were fitted under the probabilistic language version 2.26
(Stan Development Team 2021b) via inferface (Stan Development Team
2021a).

**Fit results: Hmmm, I remember there is a package showing estimations
for one parameter from different models side-by-side in one ggplot. Do
you remember?**
