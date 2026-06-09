We assume that shipments have approximately a 0.1% infestation (mu = 0.001). We include a little variance with a high kappa (e.g. 500) because there will be some randomness in the visual sampling protocol outcome (i.e., at a 95% confidence level)



# What is the actual distribution of disease in these shipments?
You see about 1 detection every 2 years with ∼4000 shipments per year:

Detections per year ≈0.5.
Shipments per year ≈4000.
So the per-shipment detection probability is
r = 0.5/4000 = 1.25x10^-4


We therefore want (μ,κ) such that
1−B(μκ,(1−μ)κ+n) / B(μκ,(1−μ)κ)  =  r  = 1.25x10^-4


## Low variability scenario
Low variability (large κ): prevalence nearly constant across shipments. Then D(p)≈nμ for tiny μ, so μ≈r/n≈4.17×10−8 (vanishingly small). This contradicts your belief that some shipments are near 1%1\%1%.
## high variabiliity scenario
High variability (small κ): most shipments near zero prevalence, a tiny tail with moderate/large p. Then μ can be much larger than r/n while still producing a tiny average detection probability, because whenever p is moderate, D(p) saturates close to 1.

I believe there is high variability -- e.g., some shipments have potatoes from all clean fields while other shipments come have many infected potatoes from an infected field


#See table for solution, where mu and kappa are solved for to produce the required alpha, beta. It also solves the likelihood of that prevalence is greater than or equal to 0.1% (the lever it is tested to), as well as 1%.


k       | u (mu)        | a (alpha) = uk        | B (beta) = (1 - u)k       | P(p >= 0.1%)  | P(p >= 1%)
0.01       0.000115137      1.151*10^-6             0.00999885                  1.23x10^-4      1.204*10^-4
0.1        0.0000678272     6.783×10⁻⁶              0.0999932                   1.1363×10⁻⁴     9.796×10⁻⁵
1          0.0000145633     1.456×10⁻⁵              0.9999854                   1.0059×10⁻⁴     6.7065×10⁻⁵
10.0       0.00000217111    2.171×10⁻⁵              9.999978                    8.8747×10⁻⁵     4.0478×10⁻⁵

Where k = 0.01, P(p >= 0.1%) is approx. 1 in 8130
Where k = 0.1, P(p >= 0.1%) is approx. 1 in 8800
Where k = 1, P(0.1%) is approx. 1 in 9941
Where k = 10, P(0.1%) is ~ 1 in 11268


# Thus, we can use k = 0.01, mu = 0.00011517 as an alternative infective scenario.