source("GibbsPGGM.R")

# Example with :
# - q = 2 outcomes
# - p = 100 predictors (10 groups of size 10) with only groups 3 and 8 half filled
# - n = 500
q = 2
p = 100
n = 500
D = matrix(0, q, p)
D[, seq(21, 30, 2)] = 2
D[, seq(71, 80, 2)] = -3
Oy = matrix(c(2, -1, -1, 2), q, q)
iOy = solve(Oy)
B = -t(D)%*%iOy
X = matrix(rnorm(n*p), n, p)
E = matrix(0, n, q)
for (i in 1:n){
  E[i, ] = rmvnorm(1, rep(0, q), sigma = iOy)
}
Y = X%*%B + E

# Sparsity with adaptative shrinkage
PGGM = GibbsPGGM(Y, X, ListHyp = list(a=10, b=1), type = "s", shr = "a")
PGGM$EstD
PGGM$EstOy

# Group-sparsity with global shrinkage
# Warning : the shrinkage hyperparameter may have to be controlled with a high initial value (to avoid numerical pitfalls)
PGGM = GibbsPGGM(Y, X, ListHyp = list(a=3, b=1), type = "gs", shr = "g", NGrp = rep(10, 10), initl = 100)
PGGM$EstD
PGGM$EstOy
PGGM$ProbInc # Probability of inclusion of each group

# Sparse-group-sparsity with adaptative shrinkage
# Warning : the shrinkage hyperparameter may have to be controlled with a high initial value (to avoid numerical pitfalls)
PGGM = GibbsPGGM(Y, X, ListHyp = list(a1=5, b1=1, a2=10, b2=1), type = "sgs", shr = "a", NGrp = rep(10, 10), initl = 100)
PGGM$EstD
PGGM$EstOy

# No sparsity (not recommended in this example !)
PGGM = GibbsPGGM(Y, X, type = "ns")
PGGM$EstD
PGGM$EstOy
