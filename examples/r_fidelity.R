# r_fidelity.R
# Fidelity check: compare the Stata nnsreg adaptation with the reference R
# implementation (NNS.reg from the NNS package) on the motorcycle data.
#
# nnsreg is a deliberate univariate subset of NNS.reg (see manuscript section
# 2.3): fixed user-chosen order, mean/median partition centers, and its own
# endpoint and central-point handling.  The two are therefore expected to
# agree in shape, not term for term.  This script quantifies that agreement.
#
# Input : moto_for_r.csv with columns time, accel, connect_fit, ols_fit,
#         where connect_fit and ols_fit are the Stata nnsreg fitted values
#         at order 3 (method(connect) and method(ols)).  Produce it with:
#           webuse motorcycle, clear
#           nnsreg accel time, order(3) method(connect) noplot
#           predict double connect_fit
#           nnsreg accel time, order(3) method(ols) noplot
#           predict double ols_fit
#           export delimited time accel connect_fit ols_fit using moto_for_r.csv, replace
# Output: prints a small table and writes moto_rfidelity.csv.
#
# Reference implementation: NNS package (Viole), version 13.1.

suppressMessages(library(NNS))

args <- commandArgs(trailingOnly = TRUE)
csv <- if (length(args) >= 1) args[1] else "moto_for_r.csv"
d <- read.csv(csv)

# Stata nnsreg connect at order 3 places 13 interior knots on these data.
stata_connect_knots <- 13

get_r_fit <- function(ord) {
    f <- NNS.reg(d$time, d$accel, order = ord, plot = FALSE)
    fx <- f$Fitted.xy
    ag <- aggregate(y.hat ~ x, data = fx, FUN = mean)  # collapse x-ties
    yhat <- ag$y.hat[match(d$time, ag$x)]
    list(yhat = yhat,
         interior = nrow(f$regression.points) - 2,
         r2 = f$R2)
}

rows <- list()
for (o in 3:6) {
    g <- get_r_fit(o)
    rows[[length(rows) + 1]] <- data.frame(
        r_order        = o,
        r_interior_knots = g$interior,
        r_R2           = round(g$r2, 3),
        cor_connect    = round(cor(g$yhat, d$connect_fit, use = "complete.obs"), 4),
        cor_ols        = round(cor(g$yhat, d$ols_fit, use = "complete.obs"), 4)
    )
}
out <- do.call(rbind, rows)

cat(sprintf("Stata nnsreg connect, order(3): %d interior knots\n\n", stata_connect_knots))
cat("R NNS.reg fitted values vs Stata nnsreg fitted values (correlation):\n")
print(out, row.names = FALSE)

write.csv(out, "moto_rfidelity.csv", row.names = FALSE)
cat("\nWrote moto_rfidelity.csv\n")
