# CO2-modeling-with-Stan
Bayesian modeling for atmospheric CO<sub>2</sub> forecasts

Data viz-replete report is `co2_writeup.pdf`. 

Using weekly atmospheric CO<sub>2</sub> measurements from the Mauna Loa Observatory in Hawaii 
(public dataset available through the Scripps Institute of Oceanography), I predict 
atmospheric carbon dioxide levels through 2058, which is 100 years since the Scripps
CO<sub>2</sub> program first began taking measurements at Mauna Loa. I also infer the year by 
which we can expect with high probability that atmospheric CO<sub>2</sub> will exceed 450 ppm — 
the threshold over which we critically reduce our chance to stabilize the average 
global temperature. I evaluate and compare a quadratic (likelihood) model and an 
exponential (likelihood) model for the task. I propose priors for their parameters, 
explaining my reasoning, and then use RStan, an imperative probabilistic programming 
language, to arrive at the parameters’ posterior predictive distributions. Having 
arrived at appropriate parameters, I use RStan to generate future predicted values for 
CO<sub>2</sub> in ppm. 

Ultimately, I predict that On January 5th, 2058, the 95% confidence interval for atmospheric 
CO<sub>2</sub> has a lower bound of **516.0349** ppm, and an upper bound of **519.9429** ppm. 
Furthermore, the 95% confidence interval bounds 450 ppm at the dates February 18th, 2034, 
and March 17th, 2035. This means that in 100 experiments with the same amount of 
samples, we could expect 95 to contain the true value in their confidence intervals. 
Moreover, it means that the year of 2034 is looking pretty dangerous: once we hit 
450 ppm, we critically diminish our chance at limiting global warming to 2 degrees Celsius. 

![450 ppm](/mauna_loa/zoom-in.png)
