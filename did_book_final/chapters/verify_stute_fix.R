# Verification script for R stute_test fix.
# BUG: StuteTest::stute_test() never calls set.seed(seed) internally.
# FIX: call set.seed() BEFORE each stute_test() call.
#
# Run this to verify reproducible p-values close to Stata.

library(StuteTest)
library(dplyr)

load(url("https://raw.githubusercontent.com/Credible-Answers/did_book/main/cc_xd_didtextbook_2025_9_30/Data%20sets/Pierce%20and%20Schott%202016/pierce_schott_didtextbook.RData"))
df <- as.data.frame(df)

cat("=================================================================\n")
cat("GQ4: Stute Tests of Linearity (order=1)\n")
cat("=================================================================\n")
cat("Stata ref: delta2001=0.042, delta2002=0.134, delta2004=0.478, delta2005=0.716\n\n")
for (yr in c(2001, 2002, 2004, 2005)) {
    set.seed(1)  # <-- FIX: set seed BEFORE each call
    res <- stute_test(df = df, Y = paste0("delta", yr), D = "ntrgap", brep = 500)
    print(res)
}

# Panel joint test
cat("\nPanel Joint Test (Stata ref: 0.0089, p=0.5400):\n")
df_long <- reshape(df,
    varying = list(
        delta = paste0("delta", c(1997, 1998, 2001, 2002, 2004, 2005)),
        deltalintrend = paste0("deltalintrend", c(1997, 1998, 2001, 2002, 2004, 2005))
    ),
    v.names = c("delta", "deltalintrend"),
    timevar = "year",
    times = c(1997, 1998, 2001, 2002, 2004, 2005),
    direction = "long"
)
df_panel_gq4 <- df_long %>% filter(year >= 2001)
set.seed(1)
res <- stute_test(df = as.data.frame(df_panel_gq4), Y = "delta", D = "ntrgap",
                  group = "indusid", time = "year", brep = 500)
print(res)

cat("\n=================================================================\n")
cat("GQ7: Stute Tests of Mean Independence (order=0)\n")
cat("=================================================================\n")
cat("Stata ref: deltalintrend1998=0.302, deltalintrend1997=0.508\n\n")
for (yr in c(1998, 1997)) {
    set.seed(1)
    res <- stute_test(df = df, Y = paste0("deltalintrend", yr), D = "ntrgap",
                      order = 0, brep = 500)
    print(res)
}

# Panel joint test
cat("\nPanel Joint Test (Stata ref: 0.0012, p=0.4720):\n")
df_pre_lt <- df_long %>% filter(year %in% c(1997, 1998))
set.seed(1)
res <- stute_test(df = as.data.frame(df_pre_lt), Y = "deltalintrend", D = "ntrgap",
                  group = "indusid", time = "year", order = 0, brep = 500)
print(res)

cat("\n=================================================================\n")
cat("GQ9: Stute Tests of Linearity with Linear Trends (order=1)\n")
cat("=================================================================\n")
cat("Stata ref: 2001=0.378, 2002=0.062, 2004=0.384, 2005=0.530\n\n")
for (yr in c(2001, 2002, 2004, 2005)) {
    set.seed(1)
    res <- stute_test(df = df, Y = paste0("deltalintrend", yr), D = "ntrgap",
                      brep = 500)
    print(res)
}

# Panel joint test
cat("\nPanel Joint Test (Stata ref: 0.0119, p=0.4020):\n")
df_post_lt <- df_long %>% filter(year >= 2001)
set.seed(1)
res <- stute_test(df = as.data.frame(df_post_lt), Y = "deltalintrend", D = "ntrgap",
                  group = "indusid", time = "year", brep = 500)
print(res)
