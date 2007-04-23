"MCMCpoissonChangepoint" <- 
        function (data, m = 1, burnin = 1000, mcmc = 1000, thin = 1, verbose = 0,
        seed = NA, c0, d0, a = NULL, b = NULL, marginal.likelihood = c("none", "Chib95"), ...)
{
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    cl <- match.call()
    if (!is.na(seed))
        set.seed(seed)
    y <- data
    n <- length(y)
    A0 <- trans.mat.prior(m = m, n = n, a = a, b = b)
    marginal.likelihood <- match.arg(marginal.likelihood)
    lambda.store <- matrix(NA, mcmc/thin, m + 1)
    P.store <- matrix(NA, mcmc/thin, (m + 1)^2)
    ps.store <- matrix(0, n, m + 1)
    s1.store <- matrix(NA, mcmc/thin, n)
    py <- rep(0, m + 1)
    pdf.P.store <- matrix(NA, mcmc/thin, m + 1)
    lambda1 <- rep(NA, m + 1)
    P1 <- matrix(NA, m + 1, m + 1)
    lambda0 <- runif(m + 1)
    P0 <- trans.mat.prior(m = m, n = n, a = 0.9, b = 0.1)
    for (iter in 1:totiter) {
        state.out <- Poisson.state.sampler(m = m, y = y, lambda = lambda0,
            P = P0)
        s1 <- state.out$s1
        ps1 <- state.out$ps1
        for (j in 1:(m + 1)) {
            ej <- as.numeric(s1 == j)
            yj <- y[ej == 1]
            nj <- length(yj)
            c1 <- sum(yj) + c0
            d1 <- nj + d0
            lambda1[j] <- rgamma(1, c1, d1)
        }
        switch <- switchg(s1)
        for (j in 1:(m + 1)) {
            switch1 <- A0[j, ] + switch[j, ]
            pj <- rdirichlet.cp(1, switch1)
            P1[j, ] <- pj
        }
        lambda0 <- lambda1
        P0 <- P1
        if (iter > burnin && (iter%%thin == 0)) {
            lambda.store[(iter - burnin)/thin, ] <- lambda1
            P.store[(iter - burnin)/thin, ] <- as.vector(t(P1))
            s1.store[(iter - burnin)/thin, ] <- s1
            ps.store <- ps.store + ps1
        }
        if (verbose > 0 && iter%%verbose == 0) {
            cat("----------------------------------------------",
                "\n")
            cat("iteration = ", iter, "\n")
            cat("lambda = ", lambda1, "\n")
            cat("Transition Matrix", "\n")
            for (i in 1:(m + 1)) cat(paste("", P1[i, ]), fill = TRUE,
                labels = paste("{", i, "}:", sep = ""), sep = ",")
        }
    }
    if (marginal.likelihood == "Chib95") {
        lambda.st <- apply(lambda.store, 2, mean)
        P.vec.st <- apply(P.store, 2, mean)
        P.st <- t(matrix(P.vec.st, m + 1, m + 1))
        density.lambda <- matrix(NA, (mcmc/thin), m + 1)
        for (i in 1:(mcmc/thin)) {
            for (j in 1:(m + 1)) {
                ej <- as.numeric(s1.store[i, ] == j)
                yj <- y[ej == 1]
                nj <- length(yj)
                c1 <- sum(yj) + c0
                d1 <- nj + d0
                density.lambda[i, j] <- dgamma(lambda.st[j],
                  c1, d1)
            }
        }
        pdf.lambda <- log(prod(apply(density.lambda, 2, mean)))
        for (g in 1:(mcmc/thin)) {
            state.out <- Poisson.state.sampler(m = m, y = y,
                lambda = lambda.st, P = P0)
            s1 <- state.out$s1
            ps1 <- state.out$ps1
            switch <- switchg(s1)
            for (j in 1:(m + 1)) {
                switch1 <- A0[j, ] + switch[j, ]
                pj <- rdirichlet.cp(1, switch1)
                P1[j, ] <- pj
                pdf.P.store[g, j] <- ddirichlet.cp(P.st[j, ],
                  switch1)
            }
            P0 <- P1
        }
        pdf.P <- log(prod(apply(pdf.P.store, 2, mean)))
        F <- matrix(NA, n, m + 1)
        like <- rep(NA, n)
        pr1 <- c(1, rep(0, m))
        for (t in 1:n) {
            py <- sapply(c(1:(m + 1)), function(i) {
                poisson.pdf(y[t], lambda.st[i])
            })
            if (t == 1) {
                pstyt1 = pr1
            }
            else {
                pstyt1 <- F[t - 1, ] %*% P.st
            }
            unnorm.pstyt <- pstyt1 * py
            pstyt <- unnorm.pstyt/sum(unnorm.pstyt)
            F[t, ] <- pstyt
            like[t] <- sum(unnorm.pstyt)
        }
        loglik <- sum(log(like))
        nprior.lambda <- nprior.P <- rep(NA, m + 1)
        nprior.lambda <- sapply(c(1:(m + 1)), function(i) {
            dgamma(lambda.st[i], c0, d0, log = TRUE)
        })
        nprior.P <- sapply(c(1:(m + 1)), function(i) {
            log(ddirichlet.cp(P.st[i, ], A0[i, ]))
        })
        prior.lambda <- sum(nprior.lambda)
        prior.P <- sum(nprior.P)
        numerator <- loglik + prior.lambda + prior.P
        denominator <- pdf.lambda + pdf.P
        logmarglike <- numerator - denominator
        if (verbose > 0) {
            cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
            cat("Log Marginal Likelihood\n")
            cat("-------------------------------------------------",
                "\n")
            cat("log(marglike)= ", logmarglike, "\n")
            cat("log(likelihood)= ", loglik, "\n")
            cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
        }
    }
    else {
        logmarglike <- marginal <- loglik <- NULL
    }
    output <- as.mcmc(lambda.store)
    varnames(output) <- paste("lambda.", 1:(m + 1), sep = "")
    attr(output, "title") <- "MCMCpoissonChangepoint Posterior Sample"
    attr(output, "y") <- data
    attr(output, "call") <- cl
    attr(output, "logmarglike") <- logmarglike
    attr(output, "loglik") <- loglik
    attr(output, "prob.state") <- ps.store/(mcmc/thin)
    return(output)
}

