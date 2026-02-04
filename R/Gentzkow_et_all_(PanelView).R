library(haven)
datos <- read_dta("~/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Gentzkow et al 2011/gentzkowetal_didtextbook.dta")
install.packages("panelView")
library(panelView)
names(datos)
length(unique(datos$cnty90))
# All of counties
panelview(prestout ~ same_treat_after_first_change, data = datos,
          index = c("cnty90", "year"),
          xlab = "Year", ylab = "County",
          display.all = TRUE, gridOff = TRUE, by.timing = TRUE,
          axis.lab.gap = c(4, 100),
          axis.lab = "both",
          axis.lab.angle = 45,
          background = "white",
          color = c("#7BB3E0", "#2A6496", "#D5D5D5"),
          legend.labs = c("Treated (Pre)", "Treated (Post)", "Control"),
          main = "Treatment Status: Gentzkow et al. (2011)",
          cex.main = 16, cex.axis.x = 8, cex.axis.y = 5,
          cex.lab = 16, cex.legend = 12)

# With 100 units
subs <- c(head(treated_ids, 60), head(control_ids, 40))
datos100 <- datos[datos$cnty90 %in% subs, ]

panelview(prestout ~ same_treat_after_first_change, data = datos100,
          index = c("cnty90", "year"),
          xlab = "Year", ylab = "County",
          display.all = TRUE, gridOff = TRUE, by.timing = TRUE,
          axis.lab.gap = c(4, 10),
          axis.lab = "both",
          axis.lab.angle = 45,
          background = "white",
          color = c("#7BB3E0", "#2A6496", "#D5D5D5"),
          legend.labs = c("Treated (Pre)", "Treated (Post)", "Control"),
          main = "Treatment Status: Gentzkow et al. (2011)",
          cex.main = 16, cex.axis.x = 8, cex.axis.y = 5,
          cex.lab = 16, cex.legend = 12)