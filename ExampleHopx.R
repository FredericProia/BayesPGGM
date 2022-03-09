source("GibbsPGGM.R")

# Example on the Hopx dataset
load("data.X.rda")
load("data.Y.Hopx.rda")
load("MAP.file.rda")
Y = as.matrix(data.Y.Hopx)
colnames(Y)[1] = "Adrenal"
X = as.matrix(data.X)
q = ncol(Y)
p = ncol(X)
n = nrow(X)

# Groups naturally arise from chromosomes
NGrp = c(74, 67, 63, 60, 39, 45, 52, 43, 31, 51, 21, 26, 33, 22, 15, 27, 18, 30, 34, 19)

# Group-sparsity with adaptative shrinkage
Est = GibbsPGGM(Y, X, ListHyp = list(a=1, b=20), type = "gs", shr="a", NGrp = NGrp, initl = 1000)
Est$EstD
Est$EstOy
Est$ProbInc # Probability of inclusion of each group

# Sparse-group-sparsity with adaptative shrinkage
Est = GibbsPGGM(Y, X, ListHyp = list(a1=3, b1=1, a2=1, b2=1), type = "sgs", shr="a", NGrp = NGrp, initl = 1000, N = 5000, Nd = 2500)
Est$EstD
Est$EstOy

# Remark: the small sample size relative to the number of covariates (29/770) is problematic.
# It would probably be better to run the algorithms on randomly chosen subsets of observations and to aggregate the results...
