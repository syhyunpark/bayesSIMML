##  12/13/2022 
# HP. summary of the simulation results 
# the results from "bsim_HTE_simulation_06182022.R"

## need to run the following chunk first to make plots.. 
################
#############
library(proto)
# Detect and prevent collisions.
# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
  # Determine height
  if (!is.null(height)) {
    # height set manually
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y - height / 2
      data$ymax <- data$y + height / 2
    }
  } else {
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y
      data$ymax <- data$y
    }
    
    # height determined from data, must be floating point constant
    heights <- unique(data$ymax - data$ymin)
    heights <- heights[!is.na(heights)]
    
    #   # Suppress warning message since it's not reliable
    #     if (!zero_range(range(heights))) {
    #       warning(name, " requires constant height: output may be incorrect",
    #         call. = FALSE)
    #     }
    height <- heights[1]
  }
  
  # Reorder by x position, relying on stable sort to preserve existing
  # ordering, which may be by group or order.
  data <- data[order(data$ymin), ]
  
  # Check for overlap
  intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
  intervals <- intervals[!is.na(intervals)]
  
  if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
    warning(name, " requires non-overlapping y intervals", call. = FALSE)
    # This is where the algorithm from [L. Wilkinson. Dot plots.
    # The American Statistician, 1999.] should be used
  }
  
  if (!is.null(data$xmax)) {
    plyr::ddply(data, "ymin", strategy, height = height)
  } else if (!is.null(data$x)) {
    data$xmax <- data$x
    data <- plyr::ddply(data, "ymin", strategy, height = height)
    data$x <- data$xmax
    data
  } else {
    stop("Neither x nor xmax defined")
  }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
  if (nrow(df) == 1) return(df)
  
  n <- nrow(df) + 1
  x <- ifelse(is.na(df$x), 0, df$x)
  if (all(is.na(df$y))) {
    heights <- rep(NA, n)
  } else {
    heights <- c(0, cumsum(x))
  }
  
  df$xmin <- heights[-n]
  df$xmax <- heights[-1]
  df$x <- df$xmax
  df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
  stacked <- pos_stackv(df, height)
  stacked$xmin <- stacked$xmin / max(stacked$xmax)
  stacked$xmax <- stacked$xmax / max(stacked$xmax)
  stacked$x <- stacked$xmax
  stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
  n <- length(unique(df$group))
  if (n == 1) return(df)
  
  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }
  
  d_height <- max(df$ymax - df$ymin)
  
  # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
  # ggplot(df, aes(n, div)) + geom_point()
  
  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidy <- match(df$group, sort(unique(df$group), decreasing = TRUE))
  
  # Find the center for each group, then use that to calculate xmin and xmax
  df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
  df$ymin <- df$y - d_height / n / 2
  df$ymax <- df$y + d_height / n / 2
  
  df
}


#' Adjust position by dodging overlaps to the side.
#'
#' @inheritParams ggplot2::position_identity
#' @param height Dodging height, when different to the height of the individual
#'   elements. This is useful when you want to align narrow geoms with wider
#'   geoms. See the examples for a use case.
#' @family position adjustments
#' @export
#' @examples
#' ggplot(mtcars, aes(factor(cyl), fill = factor(vs))) +
#'   geom_bar(position = "dodge")
#' \donttest{
#' ggplot(diamonds, aes(price, fill = cut)) +
#'   geom_histogram(position="dodge")
#' # see ?geom_boxplot and ?geom_bar for more examples
#'
#' # To dodge items with different heights, you need to be explicit
#' df <- data.frame(x=c("a","a","b","b"), y=2:5, g = rep(1:2, 2))
#' p <- ggplot(df, aes(x, y, group = g)) +
#'   geom_bar(
#'     stat = "identity", position = "dodge",
#'     fill = "grey50", colour = "black"
#'   )
#' p
#'
#' # A line range has no height:
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1), position = "dodge")
#' # You need to explicitly specify the height for dodging
#' p + geom_linerange(aes(ymin = y-1, ymax = y+1),
#'   position = position_dodge(width = 0.9))
#'
#' # Similarly with error bars:
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1), width = 0.2,
#'   position = "dodge")
#' p + geom_errorbar(aes(ymin = y-1, ymax = y+1, height = 0.2),
#'   position = position_dodge(width = 0.90))
#' }
position_dodgev <- function(height = NULL) {
  ggproto(NULL, PositionDodgeV, height = height)
}


PositionDodgeV <- ggproto("PositionDodgeV", Position,
                          required_aes = "y",
                          height = NULL,
                          setup_params = function(self, data) {
                            if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
                              warning("height not defined. Set with `position_dodgev(height = ?)`",
                                      call. = FALSE)
                            }
                            list(height = self$height)
                          },
                          
                          compute_panel = function(data, params, scales) {
                            collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
                          }
)




######

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] +
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}



#####################

library("collapse")
library("dplyr")
library("data.table")
library("tidyverse")
library("kableExtra")
library("table1")
#load("C:/Users/Danni/OneDrive - NYU Langone Health/BayesSIMML/results_testlinearcode/v7brms/all_24_scenarios_brms.rda")
load("all_24_scenarios_brms.rda")  # read-in the saved data 

ls(results.aggregated2)

results <- unlist2d(results.aggregated2, idcols = "replicate",DT = TRUE)#the first columns: scenario's id;the second column:simulation's id for each secnario
head(results)

pri_scen <- function(x){
  paste0(x["replicate"],". ","n=",x["n"],","," p=",x["p"],","," g.choice=",x["g.choice"],",",
         " m.choice=",x["m.choice"])
}
results$scenario <- apply(results,1,pri_scen)
kable(results, "html") %>% kable_styling("striped") %>% scroll_box(height = "300px")


## The expected outcome under the treatment regime
e_y <- results[,c("replicate","bsim.value","sng.bayes.value","opt.value")]
e_y <- reshape2::melt(e_y, id.var=c("replicate"))
e_y$replicate <- as.factor(e_y$replicate)

ggplot(e_y, aes(x=replicate, y=value,fill=variable)) +
  geom_boxplot(outlier.size=0.5)+ 
  labs(title="Boxplot of the expected outcome under the treatment regime",
       y = "The expected outcome under the treatment regime", x= "Scenario")+ theme(legend.position="bottom")+
  scale_fill_discrete(name = "Model", labels = c("BSIM", "BLM","True optimal value using true value of parameters"))








###########
## plots ##
###########

library(ggplot2)
library(gridExtra)
library(grid)
library(latex2exp)
library(ggthemes)


#load("famTEMsel_sim_ITR_08082020.RData")
results.aggregated2
for(i in 1:24){
  results.aggregated2[[i]] <- as.data.frame(results.aggregated2[[i]])
}

scenarios
dodge <- position_dodge(width = 0.6)

n.rep <- 100


# for simulation set "A" (m linear)
# make 4 different datasets- for each p in (5, 10) and g in (nonlinear or linear)

# 1) for p=5 and g = nonlinear  # 10, 11, 12
Value  <- c(results.aggregated2[[10]][1:n.rep, "bsim.value"],
            results.aggregated2[[10]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[10]][1:n.rep, "opt.value"],
            results.aggregated2[[11]][1:n.rep, "bsim.value"],
            results.aggregated2[[11]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[11]][1:n.rep, "opt.value"],
            results.aggregated2[[12]][1:n.rep, "bsim.value"],
            results.aggregated2[[12]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[12]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data1 <- data.frame(Value, Method, n)


# 2) for p=10 and g = nonlinear # 7, 8 ,9
Value  <- c(results.aggregated2[[7]][1:n.rep, "bsim.value"],
            results.aggregated2[[7]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[7]][1:n.rep, "opt.value"],
            results.aggregated2[[8]][1:n.rep, "bsim.value"],
            results.aggregated2[[8]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[8]][1:n.rep, "opt.value"],
            results.aggregated2[[9]][1:n.rep, "bsim.value"],
            results.aggregated2[[9]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[9]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data2 <- data.frame(Value, Method, n)

# 3) for p=5 and g = linear # 1, 2, 3
Value  <- c(results.aggregated2[[1]][1:n.rep, "bsim.value"],
            results.aggregated2[[1]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[1]][1:n.rep, "opt.value"],
            results.aggregated2[[2]][1:n.rep, "bsim.value"],
            results.aggregated2[[2]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[2]][1:n.rep, "opt.value"],
            results.aggregated2[[3]][1:n.rep, "bsim.value"],
            results.aggregated2[[3]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[3]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data3 <- data.frame(Value, Method, n)


# 4) for p=10 and g = linear # 4, 5, 6
Value  <- c(results.aggregated2[[4]][1:n.rep, "bsim.value"],
            results.aggregated2[[4]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[4]][1:n.rep, "opt.value"],
            results.aggregated2[[5]][1:n.rep, "bsim.value"],
            results.aggregated2[[5]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[5]][1:n.rep, "opt.value"],
            results.aggregated2[[6]][1:n.rep, "bsim.value"],
            results.aggregated2[[6]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[6]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data4 <- data.frame(Value, Method, n)





P1  <- ggplot(data1, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P1


P2  <- ggplot(data2, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P2


P3  <- ggplot(data3, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P3


P4  <- ggplot(data4, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P4




#theme_gray()
# A1 <- grid.arrange(P1+ theme_gray() + theme(plot.title = element_text(size=12)) + theme(legend.position="none")+
#                      ylab(TeX("$\\hat{V}(\\hat{D}^{opt}) - \\hat{V}(D^{opt})$"))
#                    ,
#                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12)) + theme(legend.position="none"),
#                    nrow=1, ncol=2)
# 
# A2 <- grid.arrange(P3+ ylab(" ") + theme(plot.title = element_text(size=12)),
#                    P4+ ylab(" ") +  theme(plot.title = element_text(size=12)),
#                    nrow=1, ncol=2)

plot1 <- grid_arrange_shared_legend(P1+ theme_gray() + theme(plot.title = element_text(size=12))+ylim(0.42, 0.51) +
                                      #ylab(TeX("$V(\\hat{D}^{opt})$")), # - V(D^{opt})$")),
                                      ylab("Expected outcome"), 
                                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12))+ylim(0.42, 0.51),
                                    P3+ ylab(" ") + theme(plot.title = element_text(size=12))+ylim(0.42, 0.51),
                                    P4+ ylab(" ") + theme(plot.title = element_text(size=12)) +ylim(0.42, 0.51), 
                                    ncol = 4, nrow = 1, position = "right")
plot1 <- grid.arrange(plot1, top = "Expected outcome under regime")
plot1




 


### PCD 

# 1) for p=5 and g = nonlinear  # 10, 11, 12
Value  <- c(results.aggregated2[[10]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[10]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[11]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[11]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[12]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[12]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data1 <- data.frame(Value, Method, n)


# 2) for p=10 and g = nonlinear # 7, 8 ,9
Value  <- c(results.aggregated2[[7]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[7]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[8]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[8]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[9]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[9]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data2 <- data.frame(Value, Method, n)

# 3) for p=5 and g = linear # 1, 2, 3
Value  <- c(results.aggregated2[[1]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[1]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[2]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[2]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[3]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[3]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data3 <- data.frame(Value, Method, n)


# 4) for p=10 and g = linear # 4, 5, 6
Value  <- c(results.aggregated2[[4]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[4]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[5]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[5]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[6]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[6]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data4 <- data.frame(Value, Method, n)





P1  <- ggplot(data1, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P1


P2  <- ggplot(data2, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P2


P3  <- ggplot(data3, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P3


P4  <- ggplot(data4, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P4

 
plot2 <- grid_arrange_shared_legend(P1+ theme_gray() + theme(plot.title = element_text(size=12))+ ylim(0.4, 1) +
                                      #ylab(TeX("$V(\\hat{D}^{opt})$")), # - V(D^{opt})$")),
                                      ylab("PCD"), 
                                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12))+  ylim(0.4, 1),
                                    P3+ ylab(" ") + theme(plot.title = element_text(size=12)) + ylim(0.4, 1),
                                    P4+ ylab(" ") + theme(plot.title = element_text(size=12)) + ylim(0.4, 1), 
                                    ncol = 4, nrow = 1, position = "right")
plot2 <- grid.arrange(plot2, top = "Proportion of correct decision (PCD)")
plot2





### Deviance
 
# 1) for p=5 and g = nonlinear  # 10, 11, 12
Value  <- c(results.aggregated2[[10]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[10]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[11]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[11]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[12]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[12]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data1 <- data.frame(Value, Method, n)


# 2) for p=10 and g = nonlinear # 7, 8 ,9
Value  <- c(results.aggregated2[[7]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[7]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[8]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[8]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[9]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[9]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data2 <- data.frame(Value, Method, n)

# 3) for p=5 and g = linear # 1, 2, 3
Value  <- c(results.aggregated2[[1]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[1]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[2]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[2]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[3]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[3]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data3 <- data.frame(Value, Method, n)


# 4) for p=10 and g = linear # 4, 5, 6
Value  <- c(results.aggregated2[[4]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[4]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[5]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[5]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[6]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[6]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data4 <- data.frame(Value, Method, n)





P1  <- ggplot(data1, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P1


P2  <- ggplot(data2, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P2


P3  <- ggplot(data3, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P3


P4  <- ggplot(data4, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P4


plot3 <- grid_arrange_shared_legend(P1+ theme_gray() + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.42) +
                                      #ylab(TeX("$V(\\hat{D}^{opt})$")), # - V(D^{opt})$")),
                                      ylab("Deviance"), 
                                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.42),
                                    P3+ ylab(" ") + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.42),
                                    P4+ ylab(" ") + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.42),
                                    ncol = 4, nrow = 1, position = "right")
plot3 <- grid.arrange(plot3, top = "Deviance")
plot3

grid.arrange(plot3, plot2, plot1)




#########################################
#########################################
# for simulation set "B" (m nonlinear)
# make 4 different datasets- for each p in (5, 10) and g in (nonlinear or linear)

# 1) for p=5 and g = nonlinear  # 10, 11, 12
Value  <- c(results.aggregated2[[22]][1:n.rep, "bsim.value"],
            results.aggregated2[[22]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[22]][1:n.rep, "opt.value"],
            results.aggregated2[[23]][1:n.rep, "bsim.value"],
            results.aggregated2[[23]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[23]][1:n.rep, "opt.value"],
            results.aggregated2[[24]][1:n.rep, "bsim.value"],
            results.aggregated2[[24]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[24]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data1 <- data.frame(Value, Method, n)


# 2) for p=10 and g = nonlinear # 7, 8 ,9
Value  <- c(results.aggregated2[[19]][1:n.rep, "bsim.value"],
            results.aggregated2[[19]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[19]][1:n.rep, "opt.value"],
            results.aggregated2[[20]][1:n.rep, "bsim.value"],
            results.aggregated2[[20]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[20]][1:n.rep, "opt.value"],
            results.aggregated2[[21]][1:n.rep, "bsim.value"],
            results.aggregated2[[21]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[21]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data2 <- data.frame(Value, Method, n)

# 3) for p=5 and g = linear # 1, 2, 3
Value  <- c(results.aggregated2[[13]][1:n.rep, "bsim.value"],
            results.aggregated2[[13]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[13]][1:n.rep, "opt.value"],
            results.aggregated2[[14]][1:n.rep, "bsim.value"],
            results.aggregated2[[14]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[14]][1:n.rep, "opt.value"],
            results.aggregated2[[15]][1:n.rep, "bsim.value"],
            results.aggregated2[[15]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[15]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data3 <- data.frame(Value, Method, n)


# 4) for p=10 and g = linear # 4, 5, 6
Value  <- c(results.aggregated2[[16]][1:n.rep, "bsim.value"],
            results.aggregated2[[16]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[16]][1:n.rep, "opt.value"],
            results.aggregated2[[17]][1:n.rep, "bsim.value"],
            results.aggregated2[[17]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[17]][1:n.rep, "opt.value"],
            results.aggregated2[[18]][1:n.rep, "bsim.value"],
            results.aggregated2[[18]][1:n.rep, "sng.bayes.value"],
            results.aggregated2[[18]][1:n.rep, "opt.value"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep), rep("True optimal", n.rep)),
                     levels = c("Index model","Linear regression", "True optimal")), 3)
n <- factor(c(rep("500", 3*n.rep), rep("1000", 3*n.rep), rep("2000", 3*n.rep)),
            levels = c("500", "1000", "2000"))
data4 <- data.frame(Value, Method, n)

 

P1  <- ggplot(data1, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P1


P2  <- ggplot(data2, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P2


P3  <- ggplot(data3, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P3


P4  <- ggplot(data4, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P4


#theme_gray()
# A1 <- grid.arrange(P1+ theme_gray() + theme(plot.title = element_text(size=12)) + theme(legend.position="none")+
#                      ylab(TeX("$\\hat{V}(\\hat{D}^{opt}) - \\hat{V}(D^{opt})$"))
#                    ,
#                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12)) + theme(legend.position="none"),
#                    nrow=1, ncol=2)
# 
# A2 <- grid.arrange(P3+ ylab(" ") + theme(plot.title = element_text(size=12)),
#                    P4+ ylab(" ") +  theme(plot.title = element_text(size=12)),
#                    nrow=1, ncol=2)

plot1 <- grid_arrange_shared_legend(P1+ theme_gray() + theme(plot.title = element_text(size=12))+ylim(0.415, 0.51) +
                                      #ylab(TeX("$V(\\hat{D}^{opt})$")), # - V(D^{opt})$")),
                                      ylab("Expected outcome"), 
                                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12))+ylim(0.415, 0.51),
                                    P3+ ylab(" ") + theme(plot.title = element_text(size=12))+ylim(0.415, 0.51),
                                    P4+ ylab(" ") + theme(plot.title = element_text(size=12)) +ylim(0.415, 0.51), 
                                    ncol = 4, nrow = 1, position = "right")
plot1 <- grid.arrange(plot1, top = "Expected outcome under regime")
plot1







### PCD 

# 1) for p=5 and g = nonlinear  # 10, 11, 12
Value  <- c(results.aggregated2[[22]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[22]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[23]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[23]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[24]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[24]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data1 <- data.frame(Value, Method, n)


# 2) for p=10 and g = nonlinear # 7, 8 ,9
Value  <- c(results.aggregated2[[19]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[19]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[20]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[20]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[21]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[21]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data2 <- data.frame(Value, Method, n)


# 3) for p=5 and g = linear # 1, 2, 3
Value  <- c(results.aggregated2[[13]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[13]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[14]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[14]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[15]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[15]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data3 <- data.frame(Value, Method, n)


# 4) for p=10 and g = linear # 4, 5, 6
Value  <- c(results.aggregated2[[16]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[16]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[17]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[17]][1:n.rep, "sng.bayes.accuracy"], 
            results.aggregated2[[18]][1:n.rep, "bsim.accuracy"],
            results.aggregated2[[18]][1:n.rep, "sng.bayes.accuracy"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data4 <- data.frame(Value, Method, n)





P1  <- ggplot(data1, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P1


P2  <- ggplot(data2, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P2


P3  <- ggplot(data3, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P3


P4  <- ggplot(data4, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P4


plot2 <- grid_arrange_shared_legend(P1+ theme_gray() + theme(plot.title = element_text(size=12))+ ylim(0.4, 1) +
                                      #ylab(TeX("$V(\\hat{D}^{opt})$")), # - V(D^{opt})$")),
                                      ylab("PCD"), 
                                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12))+  ylim(0.4, 1),
                                    P3+ ylab(" ") + theme(plot.title = element_text(size=12)) + ylim(0.4, 1),
                                    P4+ ylab(" ") + theme(plot.title = element_text(size=12)) + ylim(0.4, 1), 
                                    ncol = 4, nrow = 1, position = "right")
plot2 <- grid.arrange(plot2, top = "Proportion of correct decision (PCD)")
plot2





### Deviance

# 1) for p=5 and g = nonlinear  # 10, 11, 12
Value  <- c(results.aggregated2[[22]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[22]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[23]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[23]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[24]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[24]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data1 <- data.frame(Value, Method, n)


# 2) for p=10 and g = nonlinear # 7, 8 ,9
Value  <- c(results.aggregated2[[19]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[19]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[20]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[20]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[21]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[21]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data2 <- data.frame(Value, Method, n)


# 3) for p=5 and g = linear # 1, 2, 3
Value  <- c(results.aggregated2[[13]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[13]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[14]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[14]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[15]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[15]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data3 <- data.frame(Value, Method, n)


# 4) for p=10 and g = linear # 4, 5, 6
Value  <- c(results.aggregated2[[16]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[16]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[17]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[17]][1:n.rep, "sng.bayes.deviance"], 
            results.aggregated2[[18]][1:n.rep, "bsim.deviance"],
            results.aggregated2[[18]][1:n.rep, "sng.bayes.deviance"])
Method <- rep(factor(c(rep("Index model", n.rep), rep("Linear regression", n.rep)),
                     levels = c("Index model","Linear regression")), 3)
n <- factor(c(rep("500", 2*n.rep), rep("1000", 2*n.rep), rep("2000", 2*n.rep)),
            levels = c("500", "1000", "2000"))
data4 <- data.frame(Value, Method, n)





P1  <- ggplot(data1, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P1


P2  <- ggplot(data2, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & nonlinear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P2


P3  <- ggplot(data3, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 5~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P3


P4  <- ggplot(data4, aes(x=n, y=Value, fill=Method, color = Method)) +
  labs(x=" ", y = "Value")  + #ylim(-0.6, 0.0) +
  #geom_hline(yintercept=0, linetype="dotted", color = "gray", size=0.7) +
  geom_boxplot(width=.25, outlier.colour=NA, position = dodge,
               aes(x=n, y=Value, fill=Method, color = Method) ) +
  theme_classic() + theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="BuPu") + theme_minimal() + theme(legend.position="none")+
  ggtitle(bquote(p== 10~" & linear g")) + xlab(bquote(n)) + ylab("Value")  + scale_colour_wsj("colors6")
P4


plot3 <- grid_arrange_shared_legend(P1+ theme_gray() + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.47) +
                                      #ylab(TeX("$V(\\hat{D}^{opt})$")), # - V(D^{opt})$")),
                                      ylab("Deviance"), 
                                    P2+ theme_gray() + ylab(" ") + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.47),
                                    P3+ ylab(" ") + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.47),
                                    P4+ ylab(" ") + theme(plot.title = element_text(size=12))+ ylim(1.3, 1.47),
                                    ncol = 4, nrow = 1, position = "right")
plot3 <- grid.arrange(plot3, top = "Deviance")
plot3



grid.arrange(plot3, plot2, plot1)






 


######################################################################
## END OF THE FILE
######################################################################