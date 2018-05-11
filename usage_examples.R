# Note: When plotting to PDF, use the function cairo_pdf() instead of pdf()
# to get greek characters in the results.

# 50% credible intervals after seeing 6 successes out of 9 trials,
# with flat prior:
binomial.hpdi(n.successes=6, n.trials=9, prior.alpha=1, prior.beta=1, prob=0.5)
plot.binomial.hpdi(6, 9, 1, 1, 0.5)
binomial.pi(6, 9, 1, 1, 0.5)
plot.binomial.pi(6, 9, 1, 1, 0.5)

# 50% credible intervals after seeing 6 successes out of 9 trials,
# with prior assuming one earlier success and one earlier failure:
binomial.hpdi(6, 9, 2, 2, 0.5)
plot.binomial.hpdi(6, 9, 2, 2, 0.5)
binomial.pi(6, 9, 2, 2, 0.5)
plot.binomial.pi(6, 9, 2, 2, 0.5)

# 50% credible intervals after seeing 7 successes our of 11 trials
# (equivalent to 6/9 with prior alpha=2, beta=2):
binomial.hpdi(7, 11, 1, 1, 0.5)
plot.binomial.hpdi(7, 11, 1, 1, 0.5)
binomial.pi(7, 11, 1, 1, 0.5)
plot.binomial.pi(7, 11, 1, 1, 0.5)

# 50% credible intervals after seeing 1 successes out of 2 trials,
# with flat prior:
binomial.hpdi(1, 1, 1, 1, 0.5)
plot.binomial.hpdi(1, 1, 1, 1, 0.5)
binomial.pi(1, 1, 1, 1, 0.5)
plot.binomial.pi(1, 1, 1, 1, 0.5)

# Probability of parameter being larger than 0.5 after seeing 6
# successes out of 9 trials with flat prior:
binomial.pgt(6, 9, 1, 1, 0.5)

# Probability of parameter being larger than 0.5 after seeing 6
# successes out of 9 trials with prior assuming one earlier success
# and one earlier failure:
binomial.pgt(6, 9, 2, 2, 0.5)

# Probability of parameter being larger than 0.5 after seeing 7
# successes our of 11 trials (equivalent to 6/9 with prior alpha=2,
# beta=2):
binomial.pgt(7, 11, 1, 1, 0.5)