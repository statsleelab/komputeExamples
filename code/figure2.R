library(gridExtra)

load("docs/figure/figures.Rmd/sim_BC_v16.rdata")
load("docs/figure/figures.Rmd/sim_CC_v16.rdata")
load("docs/figure/figures.Rmd/sim_OF_v16.rdata")

fig2 <- grid.arrange(fig2.1, fig2.2, fig2.3, widths = c(1, 1, 1.23))
fig2
ggsave(file="docs/figure/figures.Rmd/Figure2.pdf", fig2, width=12, height=4)
ggsave(file="docs/figure/figures.Rmd/Figure2.png", fig2, width=12, height=4)



