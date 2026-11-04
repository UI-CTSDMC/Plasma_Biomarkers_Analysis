# install.packages("png")
# install.packages("grid")
# install.packages("patchwork")
# install.packages("ggplotify")
# install.packages("ggplot2")

library(png)
library(grid)
library(patchwork)
library(ggplotify)
library(ggplot2)

#------------------------------------------------------------
# DIRECTORIES
#------------------------------------------------------------

lme_dir <- "file_path"  # Directory containing the LME PNG files

km_dir <- "file_path"  # Directory containing the KM PNG files

out_dir <- "file_path" # Directory for saving the PNG files


#------------------------------------------------------------
# COMMON HELPERS
#------------------------------------------------------------

img_to_plot <- function(file_path) {
  
  img <- readPNG(file_path)
  
  g <- rasterGrob(
    img,
    interpolate = FALSE,
    width = unit(1, "npc"),
    height = unit(1, "npc")
  )
  
  as.ggplot(g) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

make_header <- function(text, size = 11.0) {
  
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.52,
      label = text,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = size,
      lineheight = 0.90,
      color = "black"
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(3, 3, 3, 3),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

make_label <- function(text, tertiles, biomarker_size, tertile_size) {
  
  ggplot() +
    
    # Biomarker name
    annotate(
      "text",
      x = 0.50,
      y = 0.62,
      label = text,
      hjust = 0.5,
      vjust = 0.5,
      fontface = "bold",
      size = biomarker_size,
      lineheight = 0.92,
      color = "black"
    ) +
    
    # T1
    annotate(
      "text",
      x = 0.50,
      y = 0.46,
      label = tertiles[1],
      hjust = 0.5,
      vjust = 0.5,
      fontface = "plain",
      size = tertile_size,
      color = "black"
    ) +
    
    # T2
    annotate(
      "text",
      x = 0.50,
      y = 0.36,
      label = tertiles[2],
      hjust = 0.5,
      vjust = 0.5,
      fontface = "plain",
      size = tertile_size,
      color = "black"
    ) +
    
    # T3
    annotate(
      "text",
      x = 0.50,
      y = 0.26,
      label = tertiles[3],
      hjust = 0.5,
      vjust = 0.5,
      fontface = "plain",
      size = tertile_size,
      color = "black"
    ) +
    
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(0, 5, 0, 5),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

make_blank_cell <- function() {
  
  ggplot() +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

#------------------------------------------------------------
# FIGURE 1 FILES
# Population-Averaged MOCA Over Time
# + Time to Progression Milestone (Cognition)
#------------------------------------------------------------

files_fig1 <- list(
  
  Lilly = list(
    biomarker = "pTau217\nplasma",
    tertiles = c(
      "T1: < 0.131",
      "T2: 0.131-0.174",
      "T3: > 0.174"
    ),
    lme = file.path(lme_dir, "P-Tau217 MOCA1.png"),
    km  = file.path(km_dir, "pTau217plasma_Cognitivedomain_KM.png")
  ),
  
  NFL = list(
    biomarker = "NFL",
    tertiles = c(
      "T1: < 9.877",
      "T2: 9.877-14.538",
      "T3: > 14.538"
    ),
    lme = file.path(lme_dir, "NFL MOCA1.png"),
    km  = file.path(km_dir, "NFLplasma_Cognitivedomain_KM.png")
  ),
  
  GFAP = list(
    biomarker = "GFAP",
    tertiles = c(
      "T1: < 64.478",
      "T2: 64.478-105.341",
      "T3: > 105.341"
    ),
    lme = file.path(lme_dir, "GFAP MOCA1.png"),
    km  = file.path(km_dir, "GFAPplasma_Cognitivedomain_KM.png")
  ),
  
  AB4240 = list(
    biomarker = "A\u03b242/A\u03b240\nplasma",
    tertiles = c(
      "T1: < 0.104",
      "T2: 0.104-0.111",
      "T3: > 0.111"
    ),
    lme = file.path(lme_dir, "Aß42_40 MOCA1.png"),
    km  = file.path(km_dir, "AB4240plasma_Cognitivedomain_KM.png")
  )
)

#------------------------------------------------------------
# FIGURE 2 FILES
# Population-Averaged MDS-UPDRS I Over Time
# + Autonomic Dysfunction
# + Functional Dependence
#------------------------------------------------------------

files_fig2 <- list(
  
  Lilly = list(
    biomarker = "pTau217\nplasma",
    tertiles = c(
      "T1: < 0.131",
      "T2: 0.131-0.174",
      "T3: > 0.174"
    ),
    lme     = file.path(lme_dir, "P-Tau217 MDS1.png"),
    km_auto = file.path(km_dir, "pTau217plasma_Autonomicdysfunctiondomain_KM.png"),
    km_fd   = file.path(km_dir, "pTau217plasma_Functionaldependencedomain_KM.png")
  ),
  
  NFL = list(
    biomarker = "NFL",
    tertiles = c(
      "T1: < 9.877",
      "T2: 9.877-14.538",
      "T3: > 14.538"
    ),
    lme     = file.path(lme_dir, "NFL MDS1.png"),
    km_auto = file.path(km_dir, "NFLplasma_Autonomicdysfunctiondomain_KM.png"),
    km_fd   = file.path(km_dir, "NFLplasma_Functionaldependencedomain_KM.png")
  ),
  
  GFAP = list(
    biomarker = "GFAP",
    tertiles = c(
      "T1: < 64.478",
      "T2: 64.478-105.341",
      "T3: > 105.341"
    ),
    lme     = file.path(lme_dir, "GFAP MDS1.png"),
    km_auto = file.path(km_dir, "GFAPplasma_Autonomicdysfunctiondomain_KM.png"),
    km_fd   = file.path(km_dir, "GFAPplasma_Functionaldependencedomain_KM.png")
  ),
  
  AB4240 = list(
    biomarker = "A\u03b242/A\u03b240\nplasma",
    tertiles = c(
      "T1: < 0.104",
      "T2: 0.104-0.111",
      "T3: > 0.111"
    ),
    lme     = file.path(lme_dir, "Aß42_40 MDS1.png"),
    km_auto = file.path(km_dir, "AB4240plasma_Autonomicdysfunctiondomain_KM.png"),
    km_fd   = file.path(km_dir, "AB4240plasma_Functionaldependencedomain_KM.png")
  )
)

#------------------------------------------------------------
# FIGURE 3 FILES
# Population-Averaged MDS-UPDRS III-OFF Over Time
# + Motor Complications
#------------------------------------------------------------

files_fig3 <- list(
  
  Lilly = list(
    biomarker = "pTau217\nplasma",
    tertiles = c(
      "T1: < 0.131",
      "T2: 0.131-0.174",
      "T3: > 0.174"
    ),
    lme      = file.path(lme_dir, "P-Tau217 MDS OFF1.png"),
    km_motor = file.path(km_dir, "pTau217plasma_Motorcomplicationsdomain_KM.png")
  ),
  
  NFL = list(
    biomarker = "NFL",
    tertiles = c(
      "T1: < 9.877",
      "T2: 9.877-14.538",
      "T3: > 14.538"
    ),
    lme      = file.path(lme_dir, "NFL MDS OFF1.png"),
    km_motor = file.path(km_dir, "NFLplasma_Motorcomplicationsdomain_KM.png")
  ),
  
  GFAP = list(
    biomarker = "GFAP",
    tertiles = c(
      "T1: < 64.478",
      "T2: 64.478-105.341",
      "T3: > 105.341"
    ),
    lme      = file.path(lme_dir, "GFAP MDS OFF1.png"),
    km_motor = file.path(km_dir, "GFAPplasma_Motorcomplicationsdomain_KM.png")
  ),
  
  AB4240 = list(
    biomarker = "A\u03b242/A\u03b240\nplasma",
    tertiles = c(
      "T1: < 0.104",
      "T2: 0.104-0.111",
      "T3: > 0.111"
    ),
    lme      = file.path(lme_dir, "Aß42_40 MDS OFF1.png"),
    km_motor = file.path(km_dir, "AB4240plasma_Motorcomplicationsdomain_KM.png")
  )
)

#------------------------------------------------------------
# CHECK ALL FILES
#------------------------------------------------------------

all_paths <- c(
  unlist(lapply(files_fig1, function(x) c(x$lme, x$km))),
  unlist(lapply(files_fig2, function(x) c(x$lme, x$km_auto, x$km_fd))),
  unlist(lapply(files_fig3, function(x) c(x$lme, x$km_motor)))
)

missing_files <- all_paths[!file.exists(all_paths)]

if (length(missing_files) > 0) {
  cat("Missing files:\n")
  print(missing_files)
  stop("Please fix the missing file names before running.")
}

#------------------------------------------------------------
# FIGURE 1 PANEL
#------------------------------------------------------------

make_row_fig1 <- function(item) {
  
  label_plot <- make_label(
    item$biomarker,
    item$tertiles,
    biomarker_size = 10.0,
    tertile_size   = 8.0
  )
  
  lme_plot <- img_to_plot(item$lme)
  km_plot  <- img_to_plot(item$km)
  
  label_plot + lme_plot + km_plot +
    plot_layout(widths = c(0.38, 1, 1))
}

header_row_fig1 <- make_blank_cell() +
  make_header("Population-Averaged MOCA Over Time", size = 11.0) +
  make_header("Time to Progression Milestone\n(Cognition)", size = 11.0) +
  plot_layout(widths = c(0.38, 1, 1))

fig1_panel <- (
  header_row_fig1 /
    make_row_fig1(files_fig1$Lilly) /
    make_row_fig1(files_fig1$NFL) /
    make_row_fig1(files_fig1$GFAP) /
    make_row_fig1(files_fig1$AB4240)
) +
  plot_layout(
    heights = c(0.22, 1, 1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(
        t = 4, r = 72, b = 4, l = 4, unit = "pt"
      )
    )
  )

ggsave(
  filename = file.path(
    out_dir,
    "MOCA_LME_KM_CognitiveDomain_Panel_Final_Fig1_v6.png"
  ),
  plot = fig1_panel,
  width = 29,
  height = 25,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

#------------------------------------------------------------
# FIGURE 2 PANEL
#------------------------------------------------------------

make_row_fig2 <- function(item) {
  
  label_plot <- make_label(
    item$biomarker,
    item$tertiles,
    biomarker_size = 9.0,
    tertile_size   = 7.8
  )
  
  lme_plot  <- img_to_plot(item$lme)
  auto_plot <- img_to_plot(item$km_auto)
  fd_plot   <- img_to_plot(item$km_fd)
  
  label_plot + lme_plot + auto_plot + fd_plot +
    plot_layout(widths = c(0.36, 1, 1, 1))
}

header_row_fig2 <- make_blank_cell() +
  make_header("Population-Averaged\nMDS-UPDRS I Over Time", size = 11.0) +
  make_header("Time to Progression Milestone\n(Autonomic Dysfunction)", size = 11.0) +
  make_header("Time to Progression Milestone\n(Functional Dependence)", size = 11.0) +
  plot_layout(widths = c(0.36, 1, 1, 1))

fig2_panel <- (
  header_row_fig2 /
    make_row_fig2(files_fig2$Lilly) /
    make_row_fig2(files_fig2$NFL) /
    make_row_fig2(files_fig2$GFAP) /
    make_row_fig2(files_fig2$AB4240)
) +
  plot_layout(
    heights = c(0.20, 1, 1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(
        t = 4, r = 72, b = 4, l = 4, unit = "pt"
      )
    )
  )

ggsave(
  filename = file.path(
    out_dir,
    "MDS1_LME_KM_Autonomic_FunctionalDependence_Panel_Final_Fig2_v6.png"
  ),
  plot = fig2_panel,
  width = 29,
  height = 25,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

#------------------------------------------------------------
# FIGURE 3 PANEL
#------------------------------------------------------------

make_row_fig3 <- function(item) {
  
  label_plot <- make_label(
    item$biomarker,
    item$tertiles,
    biomarker_size = 10.0,
    tertile_size   = 8.0
  )
  
  lme_plot <- img_to_plot(item$lme)
  km_plot  <- img_to_plot(item$km_motor)
  
  label_plot + lme_plot + km_plot +
    plot_layout(widths = c(0.38, 1, 1))
}

header_row_fig3 <- make_blank_cell() +
  make_header("Population-Averaged MDS-UPDRS III-OFF\nOver Time", size = 11.0) +
  make_header("Time to Progression Milestone\n(Motor Complications)", size = 11.0) +
  plot_layout(widths = c(0.38, 1, 1))

fig3_panel <- (
  header_row_fig3 /
    make_row_fig3(files_fig3$Lilly) /
    make_row_fig3(files_fig3$NFL) /
    make_row_fig3(files_fig3$GFAP) /
    make_row_fig3(files_fig3$AB4240)
) +
  plot_layout(
    heights = c(0.22, 1, 1, 1, 1)
  ) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(
        t = 4, r = 72, b = 4, l = 4, unit = "pt"
      )
    )
  )

ggsave(
  filename = file.path(
    out_dir,
    "MDSOFF_LME_KM_MotorComplications_Panel_Final_fig3_v6.png"
  ),
  plot = fig3_panel,
  width = 29,
  height = 25,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

#------------------------------------------------------------
# CONFIRM OUTPUT
#------------------------------------------------------------

cat(
  "Saved:\n",
  file.path(out_dir, "MOCA_LME_KM_CognitiveDomain_Panel_Final_Fig1_v6.png"), "\n",
  file.path(out_dir, "MDS1_LME_KM_Autonomic_FunctionalDependence_Panel_Final_Fig2_v6.png"), "\n",
  file.path(out_dir, "MDSOFF_LME_KM_MotorComplications_Panel_Final_fig3_v6.png"), "\n"
)