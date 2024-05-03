#' Produces RMSE boxplots based on the list containinng rmses and their sum for the number of dimensions of interest
#'
#' @param obj List containing rmses
#' @param plot_name Name to be given to the plot file
#' @param d Number of input dimensions
#' @param nexps Number of experiments
#' @param start Start dimension
#' @param end End dimension
#' @param width Width of the plot
#' @param height Height of the plot
#' @param is_mch If \code{TRUE}, plot MCH (this is the worst method for all real-world optimization problems)
#' @noRd
#'
#' @return saves the plot to the file named \code{obj_name}
#'
#' @examples
#' obj_name <- "ptest_unif"
#' obj <- readRDS("./rds/normalized_outputs/ptest_unif.rds")
#' rmse_boxplot(obj, obj_name, d=3, nexps=10, start=1, end=3, width=3.5, height=2, is_mch=TRUE)
#'
rmse_boxplot <- function(obj, plot_name, d, nexps, start, end, width, height, is_mch=FALSE) {
  # obj <- readRDS("18:44:41.rds")
  # plot_name <- "zezr"
  # d <- 6
  # nexps <- 30
  # start <- 3
  # end <- 6

  # obj <- Result
  color_map <- list(AG="magenta", LP="orange", MCH="cyan", SSPD="green", SEE="yellow", FG="red", Zahm="blue")
  method_abbr <- list(AG="A", LP="L", MCH="M", SSPD="P", SEE="S", FG="F", Zahm="Z")
  legend_abbr_exp <- list(AG="AG - A", LP="LP - L", MCH="MCH - M", SSPD="SSPD - P", SEE="SEE - S", FG="FG - F", Zahm="Zahm - Z")

  if(!is_mch) obj <- obj[names(obj)!="MCH"]
  methods <- names(obj)

  methods_expl <- as.vector(unlist(legend_abbr_exp[(names(legend_abbr_exp) %in% methods)]))
  Q <- lapply(obj, function(x) do.call(rbind,x))
  Q <- do.call(cbind,lapply(Q, function(x) x[,ncol(x)]))
  Q <- cbind(rep(1:d, times=rep(nexps, d)), Q)

  colnames(Q) <- c("d", methods)
  Q <- Q[((start-1)*nexps+1):((end-1)*nexps+nexps),]

  df <- as.data.frame(Q)

  melted_df <-pivot_longer(df, cols = -d, names_to = "Methods", values_to = "Value")

  melted_df$Methods <- factor(melted_df$Methods, levels = methods)

  gp <- ggplot(data=melted_df, mapping=aes(x = factor(d), y = Value)) +
    geom_boxplot(aes(fill =Methods), lwd=0.1, outlier.size = 0.01, width=0.85) +
    labs(x = "d", y = "Sum of RMSE")+
    scale_fill_manual(values = color_map[names(color_map) %in% methods], labels=methods_expl)


  gp <- gp+theme_minimal()+theme(panel.grid.major = element_line(size = 0.01), panel.grid.minor = element_line(size = 0.01), legend.position = c(0.85, 0.6), legend.key.size = unit(0.5, 'cm'),
                                 legend.text = element_text(size=4.5), axis.title=element_text(size=7), axis.text = element_text(size = 6), legend.title = element_blank())

  # Create a summary data frame to position text inside boxplots
  text_df <- melted_df %>% group_by(d, Methods) %>% summarise(text_pos = max(Value))
  dl <- length(unique(melted_df$d))
  abbr_arr <- rep(as.vector(unlist(method_abbr[names(method_abbr) %in% melted_df$Methods])), dl)

  text_df$abbr <- abbr_arr

  if(d %in% melted_df$d) text_df[text_df$d==d,]$abbr <- ""

  # Add text inside boxplots
  gp <- gp +
    geom_text(
      data = text_df,
      aes(x = as.factor(d), y = text_pos, label = abbr, group = Methods),
      vjust = -0.5,
      size = 1.5,
      show.legend = FALSE,
      position = position_dodge(width = 0.8)
    )


  vline_rep <- seq(1.5, length.out = dl-1, by = 1)

  gp <- gp + geom_vline(xintercept=vline_rep, linetype=3, linewidth=0.2, col="blue")
  print(gp)

  ggsave(paste0("fig_",plot_name, ".pdf"), plot = gp, width = width, height = height)

    }

