# Stata-Final-Draft

*SRQM Final Draft - Caroline Perissini Blasque and Nanami Kawashima - Spring 2016
*Title: The Percentage of Women in National Parliament: Social, Political and Economic Determinants 
*Dataset: Quality of Government (QOG) 2013

*Firstly, we load the Quality of Government dataset.
use data/qog2013, clear

* =======================
*  = Dependent Variable =
* =======================
*DEPENDENT VARIABLE: Percentage of Women in National Parliament (lower house) (ipu_w_lower)

*We want to analyze the percentage of women's participation in National Parliament. 
*We can use the variables "% Woman in National Parliament (upper house)" or 
*"% Women in National Parliament (lower house)". 

*Since the latter has more observations (190 against 78 for the former), 
*we choose to analyze the percentage of women in national parliament in the lower house.
lookfor ipu_w_lower

*Describe the dependent variable. 
describe ipu_w_lower 

*Now we want to label our variable. 
la var ipu_w_lower "Percent of Women in National Parliament (Lower House)"

*We rename our variable to make it easier for us. 
ren ipu_w_lower womenpar

*To see if it worked, we use the command "describe". 
d womenpar

*Summary statistics to do a brief overview of our variable. 
su womenpar, d

*Visualisation of distribution.
hist womenpar, percent bin (10)

* We can see here that the mean is 17.2%, but the standard deviation is large (11.1%) because our values
* are a bit spread. Our median (50% of our sample on both sides) is 15.5%. Also, we can see that 
* we have a positive skewness (0.7 - not very symmetrical) because many of the values are low percentages.   


*Recode the variable. We have a lot of values for the percentage of women in 
*parliament. We want to create 13 intervals of percentage for our DV. 
*It can be useful for visualization in later parts. 
recode womenpar ///
(0=0 "0%") ///
(0.1/4.9=1 "0.1%-4.9%") ///
(5/9.9=2 "5%-9.9%") ///
(10/14.9=3 "10%-14.9%") ///
(15/19.9=4 "15%-19.9%") ///
(20/24.9=5 "20%-24.9%") ///
(25/29.9=6 "25%-29.9%") ///
(30/34.9=7 "30%-34.9%") ///
(35/39.9=8 "35%-39.9%") ///
(40/44.9=9 "40%-44.9%") ///
(45/49.9=10 "45%-49.9%") ///
(50/54.9=11 "50%-54.9%") ///
(55/max=12 "55+%"), gen (womenparper) 

*We label that. 
la var womenparper "Percentage Intervals of Women in Parliament (13 groups)"

*We check if it worked. 
fre womenparper 

*We check countries included in the observations. 
fre cname


* ==========================
*  = Independent Variables =
* ==========================
*Now we will move on to our independent variables. We will analyze 8 variables in total, 
*but for the purposes of the dofile, 
*we will describe only 4 of them. 

*The variables we will be looking at will be: women’s political rights (ciri_wopol), 
*women’s economic rights (ciri_wecon), political regime (fh_ipolity2), and average 
*years of education for women (ihme_ayef). 

*We describe the independent variables 
d ciri_wopol ciri_wecon fh_ipolity2 ihme_ayef

*We will rename and label them as the following: 
ren ciri_wopol womenpol
la var womenpol "Women's Political Rights"
ren ciri_wecon womeneco
la var womeneco "Women's Economic Rights"
ren fh_ipolity2 polreg
la var polreg "Political Regime"
ren ihme_ayef yearswomen
la var yearswomen "Years of Education (Women)"

*To check if the re-naming worked, check the frequency for each variable. 
fre womenpol womeneco polreg yearswomen

*Next, we want to visualize and then delete the missing values for the values 
*(both for our DV and IVs). 
misstable pat womenpar womenpol womeneco polreg yearswomen 

*We want to drop the missing observations. 
drop if mi(womenpar, womenpol, womeneco, polreg, yearswomen)
*(23 observation were deleted). 


* ====FIRST INDEPENDENT VARIABLE (IV) ===
*Women's Political Rights (ciri_wopol)
*Here is the frequency for the variable. 
fre womenpol 

*Let's recode to dummies here. 
recode womenpol ///
(0/1 = 0 "No laws or not enforced") ///
(2/3 = 1 "Laws enforced") ///
(else = .), gen(womenpold)
la var womenpold "Women's Political Rights (0/1)"

*Then, we see the result.
fre womenpold

*We want to see our dependent variable for each value of Women's Political Rights
bys womenpold: su womenpar

*The same thing here but we'll crosstabulate the variables now, using womenparper. 
*It's useful for visualization. 
tab womenparper womenpold, col nof 

*A boxplot with our DV and this IV.  
gr hbox womenpar, over(womenpold) ///
name(WPP_POLD, replace)

*Here we see a clear difference. There is a much higher percentage of women in Parliament
*when there are laws that enforce women's political rights. 


*===== SECOND INDEPENDENT VARIABLE (IV) =====
*Women's economic rights (ciri_wecon)

*Checking the frequency.
fre womeneco 

*Let's recode to dummies again. 
recode womeneco ///
(0/1 = 0 " No laws or not enforced") ///
(2/3 = 1 "Laws enforced") ///
(else = .), gen(womenecod)
la var womenecod "Women's Economic Rights (0/1)"

fre womenecod

*Let's summarize it. 
su womenecod


*Here is our independent variable for each value of IV Women's Economic Rights.
bys womenecod: su womenpar

*Crosstabulating now. 
tab womenparper womenecod, col nof

*And a graph. 
gr hbox womenpar, over(womenecod) ///
name(WPP_ECOD, replace)


*There are many countries that do not have laws for women's economic rights. 
*Our outliers for "no laws or not enforced laws" may have a higher level of equality and 
*not be necessary to establish them. However, we do see a small improvement in the percentage
*of women in Parliament when there are such laws. 


*======THIRD INDEPENDENT VARIABLE (IV) =====
*Political Regime (fh_ipolity2)

*Here is the frequency.
fre polreg

*Since we have too many values, we do a recoding once more. 
recode polreg ///
(0/2=0 "Least democratic 0-2") ///
(2.01/4=1 "2-4") ///
(4.01/6=2 "4-6") ///
(6.01/8=3 "6-8") ///
(8.01/max=4 "Most democratic 8-10"), gen (polregime) 

fre polregime 

*Let's summarize the variable.
su polregime

*Now we can observe our dependent variable for each value of our IV Political Regime.
bys polregime: su womenpar

*The same thing here but crosstabulating the variables now, with womenparper for better visualization.
tab womenparper polregime, col nof

*A graph here. 
gr hbox womenpar, over(polregime) ///
name(WPP_PREG, replace)

*For this variable, we actually see that the extreme values are quite similar.
*Both the least democratic and the most democratic have similar values. Therefore, we 
*do not think this is going to be a significant variable for determining our DV. 



*=====FOURTH INDEPENDENT VARIABLE (IV)====
*Average years of Education (female) (ihme_ayef)

*Check frequency. 
fre yearswomen 

*To make it easier for further exploration, we will also recode this variable. 
recode yearswomen ///
(0/3=0 "1 to 3 years") ///
(3.01/6=1 "3 to 6 years") ///
(6.01/9=2 "6 to 9 years") ///
(9.01/12=3 "9 to 12 years") ///
(12.01/max=4 "12 to 15 years"), gen (yearsed) 

*See the independent variable for each value of IV Years of Education.
bys yearsed: su womenpar

*The same thing but we'll crosstabulate the variables.
tab womenparper yearsed, col nof

*Here is a graph. 
gr hbox womenpar, over(yearsed) ///
name(WPP_YED, replace)

*This is an interesting result. We see that with few years of education, 
*the mean is higher(18.6% of women in Parliament) than more educated women (17.9%)
*In the middle (3 to 6 years of education) the mean actually decreases (12% of women 
*in Parliament). 


* =========================
*  = Normality of the DV =
* =========================

*Now that we are done examining the DV and IV’s, we will go back to our DV and check 
*its normality of the distribution by looking at its skewness and kurtosis. 

*First we will look at the summary statistics again. 
su womenpar
tabstat womenpar, s(n mean sd min max)

su womenpar, d
tabstat womenpar, s(p25 median p75 iqr) 

*Visualisation of distribution.
hist womenpar, percent bin (10)
hist womenpar, kdensity

*We can see that our recoded variable is a little closer to normality, but we will 
*use the non recoded one for normality tests. 

*We want to build a histogram to have a visualisation of normality. 
hist womenpar, bin(15) normal kdensity kdenopts(lp(dash) lc(black) bw(1.5)) ///
	note("Normal distribution (solid red) and kernel density (dashed black).") ///
	name(womenpar, replace)

*We also want to have a look at the normal distribution of values of our dependent 
*variable in a histogram.
hist womenpar, percent normal ///
     name (WPNP_hist, replace) 

*We want to see the Kernel density estimate.
kdensity womenpar, normal legend(row(1)) title("") note("") ///
	name(kdenswomen, replace)	

*=== visualization of DV with IV’s ===
*Now we create boxplots to compare our variables. 
*A boxplot with our DV, years of education and women's economic rights. 
gr hbox womenpar, over (yearsed) asyvars over (womenecod) ///
name (WPP_YE_WECO, replace)

*A boxplot with our DV, years of education and women's political rights. 
gr hbox womenpar, over (yearsed) asyvars over (womenpold) ///
name (WPP_YE_WPOL, replace)

*A boxplot with our DV, political regime and women's economic rights. 
gr hbox womenpar, over (polregime) asyvars over (womenecod) ///
name (WPP_PR_WECO, replace)

*A boxplot with our DV, political regime and women's political rights. 
gr hbox womenpar, over (polregime) asyvars over (womenpold) ///
name (WPP_PR_WECO, replace)



*=== SCALARS AND STANDARD DEVIATION ===
*Now we will look at the scalars and standard devision.
*First get summary statistics.
su womenpar, d

*Results of the last command. 
ret li

* Saving the mean and standard deviation of the summarized variable.
sca de mean = r(mean)
sca de sd   = r(sd)


* Saving the 25th and 75th percentiles and computing the interquartile range (IQR), 
*as shown in the dofile from week 4. 
sca de q1  = r(p25)
sca de q3  = r(p75)
sca de iqr = q3 - q1

* Listing all saved scalars. 
sca li


*So let's check if the number of observations in this interval of one standard deviation
*(mean - 1sd and mean +1sd) 
*are close to 68% of the observations.  
count if womenpar > mean - sd & womenpar < mean + sd
di r(N), "observations out of", _N, "(" 100 * round(r(N) / _N, .01) ///
	"% of the sample are within one standard deviation from the mean."

*We have 121 observations out of 171 in this interval, 
*which represents 71% of the sample within one standard deviation. 
*Now we do the same for two standard deviations (mean -2sd, mean +2sd).
* It should come close to 95% of all observations.
count if womenpar > mean - 2 * sd & womenpar < mean + 2 * sd
di r(N), "observations out of", _N, "(" 100 * round(r(N) / _N, .01) ///
	"% of the sample) are within 2 standard deviations from the mean."

	
*We have 165 observations out of 171, which is 96% of the sample in this interval. 


*=== OUTLIERS ====
*Next let's pay attention to the most extreme values in our sample. We summarize here mild or extreme
*outliers (below Q1 and above Q3).
su womenpar if womenpar < q1 - 1.5 * iqr | womenpar > q3 + 1.5 * iqr
su womenpar if womenpar < q1 - 3 * iqr   | womenpar > q3 + 3 * iqr

*Comparing the histograms of our variable before and after recoding
hist womenpar
hist womenparper

tabstat womenpar womenparper, c(s) s(skew kurt)


*We can see that our recoded variable is a little closer to normality, but we will use 
*the non recoded one for normality tests. 

*We want to build a histogram to have a visualisation of normality. 
hist womenpar, bin(15) normal kdensity kdenopts(lp(dash) lc(black) bw(1.5)) ///
	note("Normal distribution (solid red) and kernel density (dashed black).") ///
	name(womenpar, replace)

*We can see it is skewed to the left (too many values below the mean). 
*Let's see how much it deviates from symmetry (red line). 
symplot womenpar, ti("Symmetry plot") ///
	name(womenpar_sym, replace)

*Here is another way to visualize that through the quantiles of the variable. 
*We can see the ends showing an excess of extremes. 
qnorm womenpar, ti("Normal quantile plot") ///
	name(womenpar_qnorm, replace)

*We can check skewness and kurtosis by summarizing the variable 
su womenpar, d


*Variable transformation
*Here we try to get closer to normality. 
*We use the gladder command to see several possible transformations. 
gladder womenpar, ///
name (gladder, replace)

*We checked the logarithm transformation too, which does not appear when we use the gladder command. 
gen logwomenpar = log(womenpar) 
la var logwomenpar "Women in Parliament (log units)"

hist logwomenpar, normal ///
name(logwomenpar, replace)

tabstat womenpar logwomenpar, s(n sk kurtosis min max) c(s)

*We believe the identity is a better solution, so we will not use the log tranformation. There is no 
*need to transform our original variable. 


*=== CONFIDENCE INTERVALS === 
* Here we are looking at the confidence intervals, which are directly connected to standard error, 
* which can be caused by our sample size. 

* Mean womenperpar for our full sample with a 95% Confidence Interval.
ci womenpar

* Mean for the full sample with a 99% Confidence Interval (here we have less precision with a wider interval). 
ci womenpar, level(99)

* Looking at the subsamples of the population. As in the dofile, we see that smaller samples have 
*a bigger confidence interval. 
ci womenpar in 1/10
ci womenpar in 1/100
ci womenpar in 1/170


* Confidence intervals with proportions
* Here we use the categorical variables we recorded as dummies. The binomial distributions shall apply then.
ci womenpold, binomial
ci womenecod, binomial


* These categorical variables reflect proportions (possible range of values of each category in the sample).
prop womenpold
prop womenecod

* As it happened above, the confidence intervals get wider when we have less observations. 
prop womenecod if womenpar < 20
prop womenpold if womenpar < 20


* ===================================
* = SIGNIFICANCE TESTS: DV and IVs = 
* ===================================

*Now that we have checked the normality of our DV, we will start looking at the 
*associations of our IVs with our DV. Among our 4 independent variables, we have 2
*categorical ones that we recoded into dummies (women’s economic rights and women’s 
*political rights), and we have 2 continuous one, considering the original non-recoded 
*variables (average years of education for women and political regime). 

*The significance tests will be conducted for the variables to make sure the associations 
*observed are not due to chance. We will say that if we have a p value of less than 0.05, 
*we will assume that our associations are not random. 


*=== CATEGORICAL IV 1: Women’s Political Rights : womenpold ===
*Because we will be looking at the association between our IV 1 which is categorical with 
*our categorically recoded DV, we will use the spineplot to see their association.

fre womenpold

*Since we have a continuous DV and a dummie for the IV, we use gr dot to graphically represent 
*the relation. 
gr dot womenpar, over(womenpold, sort(1) des) scale(.75) ///
    name(womenparw2, replace)
	
*In this case, we observe that the means of the percentage of women in Parliament are
*very different for when there are laws and when there are not. 

*For the chi2 test, we will recode our DV once again so that it has less categories. We have chosen 17% of women in Parliament 
*to do the recoding because it is both close to the mean (17.68) and to the medium (16.8) of our sample. 
recode womenpar ///
(0/16.9 = 0 "Low") ///
(17/60 = 1 "High") ///
(else = .), gen(womenchi2)

*The chi-squared test. 
tab womenchi2 womenpol, exp chi2  // expected frequencies
tabchi womenchi2 womenpol, noe p  // Pearson residuals
 
*And now the ttest. 
ttest womenpar, by(womenpold)

*Here we can see that we have a difference of 14% based on the existence or not of laws for 
*political rights. Moreover, since our p-value is close to 0, we can affirm that the result 
*is statistically significant. In the chi2 test we see the correlation. 


*=== CATEGORICAL IV 2: Women’s Economic Rights : womenecod ===  

fre womenecod

*Once again we have a continuous DV and a dummie, so here is the gr dot.  
gr dot womenpar, over(womenecod, sort(1) des) scale(.75) ///
    name(womenparw, replace)
	
*For laws concerning economic rights, the difference in means is not very strong. 

* The chi-squared test:
tab womenchi2 womeneco, exp chi2  // expected frequencies
tabchi womenchi2 womeneco, noe p  // Pearson residuals

*Now the ttest. 	
ttest womenpar, by(womenecod)

*Our results are statistically significant because of our p-value (0.0005). The existence or absence of
*laws cause a difference of 6% in value of the DV. From the chi-squared test we do not see any correlation. 


*=== Continuous IV 1: Political regime : polreg ===

*Summarizing our original variable. 
su polreg

*The scatterplot for both continuous variables. 

scatter womenpar polreg, ///
	name(polr_wp, replace)
	
*We do not see any evident linear pattern in this scatterplot, so we 
*do not expect a significant correlation. 

*The pwcorr test to see the significance disregarding non-available data. 
pwcorr polreg womenpar, obs sig

*The results show that there is not a strong correlation (0.16), but it
*is statistically significant with a small risk of error, since the probability level
*of the alternative to the null hypothesis is lower than 0.05, namely 0.028.


*=== Continuous IV 2: Years of education: yearswomen ===
 
*Summarizing the variable. 
su yearswomen

*The scatterplot for our continuous variables. 
scatter womenpar yearswomen, ///
	name(yw_wp, replace)

*The pwcorr test. 
pwcorr yearswomen womenpar, obs sig

*We do not see a linear pattern once again. The correlation is not strong, just 19% (0.188). 
*However it is statistically significant (p-value of 0.0138).


*Scatterplot matrixes
*--------------------

*In this section, do do a series of visualizations of the relation between our DV and
*the continuous variables we have (polreg and yearswomen); 

*The matrix including our DV and IVs. 
gr mat womenpar polreg yearswomen, ///
	name(gr_matrix, replace)
	
*Let's check the correlation again. We do the significance test using pairwise correlation to consider only valid data. 
pwcorr womenpar polreg yearswomen

*We add the star to show the results are statistically significant. 
pwcorr womenpar polreg yearswomen, star(.05)

*Scatterplots with marker labels 
*We want to build informative scatterplots to see how countries score for each IV. 
global ccode "ms(i) mlabpos(0) mlab(cname) legend(off)"

*The scatter plot showing the relation between the DV and political regime with countries. 
sc womenpar polreg, $ccode ///
	name(wp_pr, replace)

*Another scatter plot but relating the DV and years of education with the name of countries. 
sc womenpar yearswomen, $ccode ///
	name(wp_yw, replace)
	
	
*Scatterplots with histograms

*Here is another way to represent the relationship between the variables in quite strange graphs. 
*Firtsly the DV and political regime. 

sc womenpar polreg, ///
	yti("") xti("") ysc(alt) yla(none, angle(v)) xsc(alt) xla(none, grid gmax) ///
	name(plot2, replace) plotregion(style(none))
	
*Just our DV. 
tw hist womenpar, ///
	xsc(alt rev) xla(none) xti("") horiz fxsize(25) ///
	name(plot1, replace) plotregion(style(none))

*The IV political regime. 
tw hist polreg, ///
	ysc(alt rev) yla(none, nogrid angle(v)) yti("") xla(,grid gmax) fysize(25) ///
	name(plot3, replace) plotregion(style(none))

*Let's combine everything. 
gr combine plot1 plot2 plot3, ///
	imargin(0 0 0 0) hole(3) ysiz(5) xsiz(5) ///
	name(wp_polreg, replace)

*As in the dofile, we care more about the final results, so the previous one are dropped. 
gr drop plot1 plot2 plot3
gr di wp_polreg


*Now the same process but for women's years of education. 
sc womenpar yearswomen, ///
	yti("") xti("") ysc(alt) yla(none, angle(v)) xsc(alt) xla(none, grid gmax) ///
	name(plot2, replace) plotregion(style(none))

*Our DV. 
tw hist womenpar, ///
	xsc(alt rev) xla(none) xti("") horiz fxsize(25) ///
	name(plot1, replace) plotregion(style(none))

*the representation of our IV years of education for women. 
tw hist yearswomen, ///
	ysc(alt rev) yla(none, nogrid angle(v)) yti("") xla(,grid gmax) fysize(25) ///
	name(plot3, replace) plotregion(style(none))

*Putting them all together. 
gr combine plot1 plot2 plot3, ///
	imargin(0 0 0 0) hole(3) ysiz(5) xsiz(5) ///
	name(wp_years, replace)

*The results are more important. 
gr drop plot1 plot2 plot3
gr di wp_years

	
*Scatterplots with smoothed lines

*Another type of visualization, but in here we assess the quality of a linear fit. 

*For political regime. 
lowess womenpar polreg, ///
	name(womenpar_pol, replace)

*For years of education. 
lowess womenpar yearswomen, ///
	name(womenpar_yearswomen, replace)
	
*With these smoothed lines, we can see the relation is actually not very linear, 
*especially in the case of years of education. 

*Before going to the regressions, we use the macro to remove the legend and dash
*the regression line of our linear fit. 
global ci "legend(off) lp(dash)"



* =================
* REGRESSION MODELS 
* =================

*Since we have variables with units that are not very large, we do not believe that logistic regression
*is better than the linear one. 

*Simple linear regressions
*--------------------------


*(1) Women in Parliament and Women's Political Rights
*----------------------------------------------------

* Here is the visual inspection of our two variables.
sc womenpar womenpol, $ccode ///
    legend(off) yti("Women in Parliament") ///
    name(wepa_wopol, replace)
	
* With the linear fit.
tw (sc womenpar womenpol, $ccode) (lfit womenpar womenpol, $ci), ///
    yti("Women in Parliament (Political rights)") ///
    name(wepa_wopol2, replace)

* And with the 95% confidence intervals.
tw (sc womenpar womenpol, $ccode) (lfitci womenpar womenpol, $ci), ///
    yti("Women in Parliament (Political rights)") ///
    name(wepa_wopol3, replace)

* Now the regression.
reg womenpar womenpol

* We see here that the F value and the specific p values are good, near zero, 
* showing that the results are statistically significant. The coefficient of determination, r-squared, is
* quite high (0.46). Furthermore, we can observe that the variable womenpol has a significant positive effect
* on womenpar, with an increase of 15.75 for each unit of womenpar.  

*Plotting the results. 

* The simple residuals-versus-fitted plot.
rvfplot, yline(0) ///
    name(rvfplot, replace)

* We see the fitted values.
cap drop yhat
predict yhat

* And the residuals.
cap drop r
predict r, resid

* Plotting residuals against predicted values of the IV.
sc r yhat, yline(0) $ccode ///
    name(rvfplot2, replace)

* Plotting our DV with the observed and predicted values of the IV.
sc womenpar womenpol || conn yhat womenpol, ///
    name(dv_yhat, replace)
	
* We transformed our variable to sqrt since we thought it was more quadratic,
* but we actually achieved a worst fit and a lower r-squared value. Therefore, 
* we keep womenpol as it is.  
	
	
* (2) Women in Parliament and Women's Economic Rights
*-----------------------------------------------------

* We start with the visual inspection of the fit. 
tw (sc womenpar womeneco, $ccode) (lfit womenpar womeneco, $ci), ///
    name(women_eco, replace)
		
* With the linear fit.
tw (sc womenpar womeneco, $ccode) (lfit womenpar womenpol, $ci), ///
    yti("Women in Parliament (Political rights)") ///
    name(wepa_wopol2, replace)

* The 95% confidence intervals.
tw (sc womenpar womeneco, $ccode) (lfitci womenpar womenpol, $ci), ///
    yti("Women in Parliament (Political rights)") ///
    name(wepa_wopol3, replace)

*We do the regression model of the linear form.
reg womenpar womeneco

* Here it is not such a good result as the one before. We do not see much linearity. 
* However, the F and p values are also close to zero. The r-squared is lower (0.18), which means
* womeneco explains womenpar less than womenpol. For each change in unit of womenpar, we see a positive change 
* of 5.28 in womeneco. 

*Plotting the residuals. 
*Simple residuals-versus-fitted plot.
rvfplot, yline(0) ///
    name(rvfplot, replace)
 
* The fitted values.
cap drop yhat
predict yhat
 
* The residuals.
cap drop r
predict r, resid
 
* Plotting residuals against predicted values of IV.
sc r yhat, yline(0) $ccode ///
    name(rvfplot3, replace)
 
* Plotting DV with observed and predicted values of the IV.
sc womenpar womeneco || conn yhat polreg, ///
    name(wom_eco, replace)


* (3) Women in Parliament and Political Regime
* --------------------------------------------
 
* Visual fit.
sc womenpar polreg, $ccode ///
    legend(off) yti("Women in Parliament (percentage)") ///
    name(womenpar_polreg, replace)
 
* Linear fit.
tw (sc womenpar polreg, $ccode) (lfit womenpar polreg, $ci), ///
    yti("Women in Parliament (percentage)") ///
    name(womenpar_polreg2, replace)
 
* Adding 95% confidence interval.
tw (sc womenpar polreg, $ccode) (lfitci womenpar polreg, $ci), ///
    yti("Women in Parliament (percentage)") ///
    name(womenpar_polreg, replace)
 
reg womenpar polreg

* From the regression model, we see that political regime actually explains very little 
* of the variance in women in Parliament (R-squared is 0.09). Nonetheless, the results are significant in
* a statistical perspective, with F and p values close to zero again. 1 unit of womenpar accounts for 1.15 
* variation of polreg. 

*plotting the results here. 
* Simple residuals-versus-fitted plot.
rvfplot, yline(0) ///
    name(rvfplot, replace)
 
* Fitted values.
cap drop yhat
predict yhat
 
* The residuals.
cap drop r
predict r, resid
 
* Plotting residuals against predicted values of IV.
sc r yhat, yline(0) $ccode ///
    name(rvfplot2, replace)
 
* Plotting DV with the observed and predicted values of the IV.
sc womenpar polreg || conn yhat polreg, ///
    name(dv_yhat, replace)
 


* (4) Women in Parliament and Years of Education
* ------------------------------------------------

* Visual fit.
sc womenpar yearswomen, $ccode ///
    legend(off) yti("Women in Parliament (percentage)") ///
    name(womenpar_yearswomen, replace)
 
* Linear fit.
tw (sc womenpar yearswomen, $ccode) (lfit womenpar yearswomen, $ci), ///
    yti("Women in Parliament (percentage)") ///
    name(womenpar_yearswomen2, replace)
 
* 95% CI.
tw (sc womenpar yearswomen, $ccode) (lfitci womenpar yearswomen, $ci), ///
    yti("Women in Parliament (percentage)") ///
    name(womenpar_yearswomen, replace)
 
 
* Regression model to see the predicted effect. 
reg womenpar yearswomen

* Here our F and p values are 0.016. The r-square is low (0.054), so we are not able to explain much of
* the DV with this variable. For a change in unit of DV we have a variance of 0.725. 
 
 
* Plotting regression results
* ---------------------------
 
* Simple residuals-versus-fitted plot.
rvfplot, yline(0) ///
    name(rvfplot, replace)
 
* Get fitted values.
cap drop yhat
predict yhat
 
* Get residuals.
cap drop r
predict r, resid
 
* Plot residuals against predicted values of IV.
sc r yhat, yline(0) $ccode ///
    name(rvfplot2, replace)
 
* Plot DV with observed and predicted values of IV.
sc womenpar polreg || conn yhat polreg, ///
    name(dv_yhat, replace)



* Multiple linear regression
* --------------------------

* Here we add our other three variables: unemployment, gender equality and GDP.
* Renaming variables.
renvars wdi_ue pwt_rgdpch wef_gend \ unemployment gdp gender

* We check the missing values and delete them. 
misstable pat unemployment gdp gender
drop if mi(unemployment, gdp, gender)

* Since we are dealing with GDP, we change it to log.
gen log_gdp = ln(unna_gdp / unna_pop)
la var log_gdp "Real GDP/capita (constant USD, logged)"

* Now let's add all variables to our multiple regression model. 
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

reg womenpar log_gdp unemployment

* With a cleaner and better output.
leanout: 

* First of all, we see that even though the F value is close to zero, the p values for womeneco, polreg,
* GDP and unemployment are too high and, thus, we cannot draw significant considerations from the observed coefficients because they
* can be random. We do see, however, that we have a high r-squared value (0.65).  
* Moreover, womenpol and yearswomen have a statistically significant result. Womenpol has a high coefficient (13.09) and 
* is related to the value of our DV by these units. On the other side, years of education actually have a modest negative relation to our DV (-0.73). 

* In the multiple regression model, only with womenpol and gender the p value remains close to zero. From the moment we add womeneco, 
* the p value becomes too high. The same happens with polreg, yearswomen, log_gdp and unemployment.  

reg womenpar womenpol gender 

* We did not run the beta coefficients because the do file advised against it, due to 
* the fact they are controversial. 



* Dummy (for control). 
*----------------------

* Here we shall use a variable called Electoral Process to create a dummy. 
fre fh_ep

* Creating the dummy. 
recode fh_ep ///
(0/6 = 0 "Not free nor fair elections") ///
(7/12 = 1 "Free and fair elections") ///
(else = .), gen(election)
la var election "Electoral Process (0/1)"

* Regression line for our dummy.
tw (sc yhat election) (lfit yhat election), xlab(0 "Low" 1 "High") ///
	name(reg_elec, replace)

* T-test and regression results for dummy.
ttest yhat, by(election)
reg yhat i.election


* ===========================
* = REGRESSION DIAGNOSTICS =
* ===========================

* (1) Standardized residuals
* --------------------------

* We want to store the residuals. 
cap drop r
predict r, resid

* And look at their normality. 
kdensity r, norm legend(off) ti("") ///
    name(diag_kdens, replace)

* Here we check the homoskedasticity (constant variance) of the residuals x fitted values for our DV. 
rvfplot, yline(0) ms(i) mlab(cname) name(diag_rvf, replace)

* We see that there is some variance with some outliers. 

* Storing the residuals and making them more homoskedastic than before.
cap drop rsta
predict rsta, rsta

* We have the outliers for two standard deviations (95% CI) in a more homoskedastic pattern.
sc rsta yhat, yline(-2 2) || sc rsta yhat if abs(rsta) > 2, ///
    ylab(-3(1)3) mlab(cname) legend(lab(2 "Outliers")) ///
    name(diag_rsta, replace)
	
* We still have values that are fall out of the standardization, namely Nepal, Macedonia, 
* Botswana, Belize and Sri Lanka. They are regionally diverse and small countries, and may have 
* specific variables influencing the number of women in Parliament.    


* Heteroskedasticity
* ------------------

* Let's have a look at each of the IVs now to see how they impact the model. 

* For the 1st IV: womenpol
sc r womenpol, ///
	yline(0) mlab(cname) legend(lab(2 "Outliers")) ///
	name(diag_edu1, replace)

* The 2nd IV: womeneco
sc r womeneco, ///
	yline(0) mlab(cname) legend(lab(2 "Outliers")) ///
	name(diag_edu2, replace)
	
* The 3rd IV: polreg
sc r polreg, ///
	yline(0) mlab(cname) legend(lab(2 "Outliers")) ///
	name(diag_edu3, replace)
	
*The 4th IV: yearswomen
sc r yearswomen, ///
	yline(0) mlab(cname) legend(lab(2 "Outliers")) ///
	name(diag_edu4, replace)
	
* The 5th IV: log_gdp
sc r log_gdp, ///
	yline(0) mlab(cname) legend(lab(2 "Outliers")) ///
	name(diag_edu5, replace)
	
*The 6th IV: unemployment
sc r unemployment, ///
	yline(0) mlab(cname) legend(lab(2 "Outliers")) ///
	name(diag_edu6, replace)

*The 7th IV: gender
sc r gender, ///
	yline(0) mlab(cname) legend(lab(2 "Outliers")) ///
	name(diag_edu7, replace)

* We can see that womenpol has a big impact for the configuration of the values we have seen for our DV in the graph. 

* With the LOWESS curve to compare the clouds, seing how variance occurs as a result of the IV
* by looking at the curve deviating from the axis. 
* For the 1st IV: womenpol
lowess rsta womenpol, bw(.5) yline(0) ///
	name(diag_edu1, replace)
	

* The 2nd IV: womeneco
lowess rsta womeneco, bw(.5) yline(0) ///
	name(diag_edu2, replace)
	
	
* The 3rd IV: polreg
lowess rsta womenpol, bw(.5) yline(0) ///
	name(diag_edu3,	replace)


*The 4th IV: yearswomen
lowess rsta yearswomen, bw(.5) yline(0) ///
	name(diag_edu4, replace)
	
* The 5th IV: log_gdp
lowess rsta log_gdp, bw(.5) yline(0) ///
	name(diag_edu5, replace)
	
*The 6th IV: unemployment
lowess rsta unemployment, bw(.5) yline(0) ///
	name(diag_edu6, replace)
	
*The 7th IV: gender
lowess rsta gender, bw(.5) yline(0) ///
	name(diag_edu7, replace)
	



* (3) Variance inflation and interaction terms
* --------------------------------------------

* Here we check for multicollinearity with the Variance Inflation Factor (VIF), in order to make sure we 
* are not considering twice the same factors. 
vif 

* We see they are all below 10, which is tolerable, except for womeneco that is incredibly high. Moreover, 
* for the interaction womenecoXgender the value is also very high.  

* Here we account for most of the possible interactions. We do not consider all of them due to the small
* influence they have on the regression model and because of the number of statistically irrelevant values. 
  
* (1)For womenpol and womeneco.
gen womenpolXwomeneco = womenpol * womeneco
la var womenpolXwomeneco "Women's political rights * Economic rights"

* Original regression model.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* Regression model with the interaction term.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender womenpolXwomeneco

* The f value is close to zero, which means statistical significance but all variables but womenpol and gender
* are not significant at all, with a p value of more than 0.05. Therefore, we cannot draw any meaningful conclusions 
* about them, but we do see that the coefficient for womenpol goes down when we account for the interaction with womeneco (in about 3 units).  

* (2)For womenpol and polreg
gen womenpolXpolreg = womenpol * polreg
la var womenpolXpolreg "Women's political rights * Political regime"

* Original regression model.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* Regression model with the interaction term.
reg  womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender womenpolXpolreg 

* When we account for this interaction, we can see that womenpol decreases a bit, as well as gender equality. 
* Years of education has a minor reduction too. These are the only statistically significant values. 

* (3) For womeneco and polreg
gen womenecoXpolreg = womeneco * polreg
la var womenecoXpolreg "Women's economic rights * political regime"

* Original regression model.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* Regression model with the interaction term.
reg  womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender womenecoXpolreg

* In this regression, we do not see any change in the variable that have a p value below 0.05. 

* (4) For womeneco and yearswomen.
gen womenecoXyearswomen = womeneco * yearswomen
la var womenecoXyearswomen "Women's economic rights * Education"

* Original regression model.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* Regression model with the interaction term.
reg  womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender womenecoXyearswomen

*When we take into consideration this interaction, we observe that womenpol and yearswomen suffer a very small reduction, while 
* gender actually goes up by 3 points. F value is close to zero, but these are the only significant p values. 

* (5) For polreg and yearswomen.
gen polregXyearswomen = polreg * yearswomen
la var polregXyearswomen "Political regime * Education"

* Original regression model.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* Regression model with the interaction term.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender polregXyearswomen

* Once again we have few p values below 0.05. We see that womenpol barely changed, while the value for polreg became significant
* in statistical terms, with a negative relation to the DV (-1.68). Moreover, yearswomen decreased in its negative relation to the DV.
* Gender went up by 4 points, on the other hand.  

* (6) Womenpol and gender equality. 
gen womenpolXgender = womenpol * gender
la var womenpolXgender "Women's political rights * Gender Equality"

* Original regression model.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* Regression model with the interaction term.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender womenpolXgender

* With this interaction, womenpol and gender become statistically insignificant. Yearswomen goes down slightly in the negative relation. 

* (7) Womeneco and gender equality. 
gen womenecoXgender = womeneco * gender
la var womenecoXgender "Women's economic rights * Gender Equality"

* Original regression model.
reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* Regression model with the interaction term.
reg  womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender womenecoXgender

* Womenpol and yearswomen fall a bit. Gender drops a lot (from 80 to 60) here, but the rest is all 
* statistically insignificant  for us.


* = EXPORT MODEL RESULTS =
* ========================

* Erasing the previous ones.
eststo clear

* Model 1: Baseline. 
eststo M1: qui reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender

* In a cleaner way. 
leanout:

* Model 2: With the dummy for control.
eststo M2: qui reg womenpar womenpol womeneco polreg yearswomen log_gdp gender unemployment election

* Simplified form again. 
leanout:

* Model 3: Adding the interaction between womeneco and gender equality, the one that had such a big vif. 
eststo M3: qui reg womenpar womenpol womeneco polreg yearswomen log_gdp unemployment gender election womenecoXgender

* Re-read, in simplified form.
leanout: 

* Compare all models on screen.
esttab M1 M2 M3, lab b(1) se(1) sca(rmse) ///
    mti("Baseline" "Control" "Interaction")

* Export all models for comparison and reporting.
esttab M1 M2 M3 using finaldraft.txt, replace /// 
	lab b(1) se(1) sca(rmse) ///
    mti("Baseline" "Controls" "Interactions")
	
* Here is clear that lack of statistical significance for the values we have found is a problem to draw sound 
* conclusions from this model. We only have good p values for women's political rights, gender equality and 
* women's years of education, so we cannot say much about the other variables and the model as a whole. 
* If we consider women's political rights, we can see that there is an impact of this IV on the DV. 
* We expected this due to the very nature of the variable, which accounts for women being able to vote,
* run for office and organize in society, basic requisites for becoming a representative at the National 
* Parliament. As for gender equality, the variable involves economic participation, opporunity, political 
* empowerment, educational attainment and health and well-being of women in society. In our regression, we 
* could notice that it had a very significant effect on the DV, with a large coefficient. However, it introduced 
* a problem of multicollinearity, since it took into consideration factors that were in our other variables. 
* Due to that fact, we tried to check for the interaction on our Model 3, and we got a more trustworthy coefficient there. 
* It can be seen that the values of the coefficient are always high and, therefore, this variable should have a large impact 
* on the value of women in Parliament.
* The results for years of education actually show a negative predicted impact of the years of education a women has to the 
* percentage of women in Parliament. From our previous hypothesis, it is the contrary to what we expected to see. 
* Thus, considering the results, we see that women's political rights and gender equality are important variables that 
* explain quite a lot of the DV, while the years of education has a negative impact, but it is a milder value when compared to the other two. 
* We admit that our model does not explain much about the percentage of women in Parliament, but we could see some interesting results. 
* We were constrained by the data we had and, most importantly, the fact that many important variables we wanted to use in our model 
* did not have observations for the year 2009, the only one we had observations for our DV, limiting our choices.
* We believe that better indicators for the variables we saw in the theory could have had resulted in more satisfying conclusions. 
* For instance, a variable considering quotas for women in politics could have proved to be very explanatory, having in mind
* the case of Rwanda, which has a high number of women in Parliament mainly due to the quotas introduced after the genocide in the 1990s, 
* according to literature.  
 
 *End of final draft. Thank you very much for the semester! :)
