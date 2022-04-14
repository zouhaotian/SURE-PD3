Explanation of R code files:

(1) Run surepd_data_query.R for data query and LOWESS curves (Figure 1). The original dataset was not provided but we included the processed dataset in the dataset/sure_pd folder.

(2) Run surepd_unidim_multidim_slope.R for unidimensional and multidimensional IRT models with random slope for pre-levodopa dataset.

(3) Run surepd_post_unidim_multidim_slope.R for unidimensional and multidimensional IRT models with random slope for levodopa therapy dataset.

(4) Run surepd_unidim_multidim.R for unidimensional and multidimensional IRT models without random slope for pre-levodopa dataset.

(5) Run surepd_post_unidim_multidim.R for unidimensional and multidimensional IRT models without random slope for levodopa therapy dataset.

(6) Run surepd_unidim_multidim_spline.R for unidimensional and multidimensional spline IRT models without random slope for pre-levodopa dataset (spline knot at 1 year).

(7) Run surepd_unidim_multidim_cubic.R for unidimensional and multidimensional cubic IRT models without random slope for pre-levodopa dataset (with linear, quadratic, cubic orders of time).

(8) Run surepd_unidim_multidim_slope_spline.R for unidimensional and multidimensional spline IRT models with random slope for pre-levodopa dataset (spline knot at 1 year).

(9) Run surepd_unidim_multidim_slope_cubic.R for unidimensional and multidimensional cubic IRT models with random slope for pre-levodopa dataset (with linear, quadratic, cubic orders of time).

(10) Run surepd_progression_plot.R generate the progression plot.

(11) Run SurePD_diagostic.R and SurePD_post_diagostic.R to generate diagnostic plots in the Supplementary Material (Figures S3a-S10).

(12) Run surepd_summary.R and surepd_summary_pre_post.R to generate summary statistics in Table 1.

(13) Optional: run surepd_lmm.R and surepd_post_lmm.R to fit linear mixed models for sum of scores for pre- and levodopa dataset.

