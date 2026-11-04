# ================================
# 1. Load libraries
# ================================
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# Run once if patchwork is not installed:
# install.packages("patchwork")


# ================================
# 2. Import SAS dataset
# ================================
filepath <- "file_path"  # Path where the SAS dataset is located

data <- read_sas(
  file.path(filepath, "combined_biomarker_with_clinical.sas7bdat")
)

# Check imported variable names
names(data)


# ================================
# 3. Define subgroup order
# ================================
group_order <- c(
  "Sporadic PD",
  "Healthy Control"
)

group_display <- c(
  "Sporadic PD" = "sPD NSD+",
  "Healthy Control" = "HC NSD-"
)


# ================================
# 4. Main cohort filter
#
# Keep:
#   - Early visits only
#   - sPD participants with NSD = 1
#   - HC participants with NSD = 0
#
# Biomarker-specific date filtering
# is performed inside the plotting
# function.
# ================================
plot_data <- data %>%
  mutate(
    subgroup_definition =
      str_squish(as.character(subgroup_definition)),
    
    visit_range =
      str_squish(as.character(visit_range)),
    
    COHORT_num =
      suppressWarnings(as.numeric(as.character(COHORT))),
    
    NSD_num =
      suppressWarnings(as.numeric(as.character(NSD)))
  ) %>%
  filter(
    visit_range == "Early (BL-V02)"
  ) %>%
  filter(
    (
      COHORT_num == 1 &
        subgroup_definition == "Sporadic PD" &
        NSD_num == 1
    ) |
      (
        COHORT_num == 2 &
          subgroup_definition == "Healthy Control" &
          NSD_num == 0
      )
  ) %>%
  mutate(
    subgroup_definition = factor(
      subgroup_definition,
      levels = group_order
    )
  )


# ================================
# 5. Biomarkers and corresponding
#    assay date variables
# ================================
biomarker_info <- tibble(
  var = c(
    "Lilly_ptau217p",
    "IU_GFAP_plasma_mean",
    "IU_NFL_plasma_mean",
    "C2N_AB40",
    "C2N_AB42",
    "C2N_AB4240"
  ),
  
  date_var = c(
    "Lilly_ptau217p_date",
    "IU_Plasma_date",
    "IU_Plasma_date",
    "C2N_AB_date",
    "C2N_AB_date",
    "C2N_AB_date"
  ),
  
  short_label = c(
    "pTau217 plasma",
    "GFAP",
    "NFL",
    "Aβ40 plasma",
    "Aβ42 plasma",
    "Aβ42/Aβ40 plasma"
  ),
  
  ylab = c(
    "pTau217",
    "GFAP",
    "NFL",
    "Aβ40",
    "Aβ42",
    "Aβ42/Aβ40"
  )
)


# ================================
# 6. Box-and-whisker plotting
#    function
#
# Each biomarker requires:
#   1. Non-missing corresponding date
#   2. Non-missing biomarker value
#
# Box:
#   Median and IQR
#
# Whiskers:
#   Most extreme values within
#   1.5 × IQR
# ================================
make_boxplot <- function(
    df,
    biomarker_var,
    biomarker_date_var,
    biomarker_label,
    ylab_text) {
  
  # ----------------------------
  # Biomarker-specific analysis
  # data
  # ----------------------------
  df_plot <- df %>%
    filter(
      !is.na(.data[[biomarker_date_var]]),
      !is.na(.data[[biomarker_var]])
    )
  
  # ----------------------------
  # Check for multiple records
  # per participant
  # ----------------------------
  duplicate_check <- df_plot %>%
    count(
      PATNO,
      subgroup_definition,
      name = "n_records"
    ) %>%
    filter(n_records > 1)
  
  if (nrow(duplicate_check) > 0) {
    warning(
      paste0(
        biomarker_label,
        ": Some participants have more than one eligible record. ",
        "Review the duplicate check before finalizing the figure."
      )
    )
  }
  
  # ----------------------------
  # Participant Ns by subgroup
  # ----------------------------
  biomarker_counts <- df_plot %>%
    group_by(subgroup_definition) %>%
    summarise(
      n_patno = n_distinct(PATNO),
      .groups = "drop"
    )
  
  count_frame <- tibble(
    subgroup_definition = factor(
      group_order,
      levels = group_order
    )
  ) %>%
    left_join(
      biomarker_counts,
      by = "subgroup_definition"
    ) %>%
    mutate(
      n_patno = replace_na(n_patno, 0),
      
      group_lbl = paste0(
        unname(
          group_display[
            as.character(subgroup_definition)
          ]
        ),
        "\n(n = ",
        n_patno,
        ")"
      )
    )
  
  # Add display labels to plotting data
  df_plot <- df_plot %>%
    left_join(
      count_frame %>%
        select(
          subgroup_definition,
          group_lbl
        ),
      by = "subgroup_definition"
    ) %>%
    mutate(
      group_lbl = factor(
        group_lbl,
        levels = count_frame$group_lbl
      )
    )
  
  # ----------------------------
  # Missing biomarker counts
  # among participants with a
  # corresponding assay date
  # ----------------------------
  missing_counts <- df %>%
    filter(
      !is.na(.data[[biomarker_date_var]]),
      is.na(.data[[biomarker_var]])
    ) %>%
    group_by(subgroup_definition) %>%
    summarise(
      n_missing = n_distinct(PATNO),
      .groups = "drop"
    )
  
  missing_frame <- tibble(
    subgroup_definition = factor(
      group_order,
      levels = group_order
    )
  ) %>%
    left_join(
      missing_counts,
      by = "subgroup_definition"
    ) %>%
    mutate(
      n_missing = replace_na(n_missing, 0)
    )
  
  # ----------------------------
  # QC summary
  # ----------------------------
  qc_summary <- count_frame %>%
    select(
      subgroup_definition,
      n_plotted = n_patno
    ) %>%
    left_join(
      missing_frame,
      by = "subgroup_definition"
    ) %>%
    mutate(
      biomarker = biomarker_label,
      subgroup = unname(
        group_display[
          as.character(subgroup_definition)
        ]
      )
    ) %>%
    select(
      biomarker,
      subgroup,
      n_plotted,
      n_missing
    )
  
  # ----------------------------
  # Create box-and-whisker plot
  # ----------------------------
  p <- ggplot(
    df_plot,
    aes(
      x = group_lbl,
      y = .data[[biomarker_var]],
      fill = subgroup_definition
    )
  ) +
    geom_boxplot(
      width = 0.58,
      linewidth = 0.55,
      outlier.shape = 16,
      outlier.size = 1.8,
      outlier.alpha = 0.65,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "Sporadic PD" = "#7570b3",
        "Healthy Control" = "#a6761d"
      ),
      drop = FALSE
    ) +
    scale_x_discrete(
      drop = FALSE
    ) +
    scale_y_continuous(
      name = ylab_text,
      limits = c(0, NA),
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +
    labs(
      title = biomarker_label,
      x = NULL
    ) +
    theme_minimal(
      base_size = 11
    ) +
    theme(
      legend.position = "none",
      
      panel.grid.major.x =
        element_blank(),
      
      panel.grid.minor =
        element_blank(),
      
      axis.text.x = element_text(
        size = 10,
        face = "bold",
        margin = margin(t = 5)
      ),
      
      axis.text.y = element_text(
        size = 9
      ),
      
      axis.title.y = element_text(
        size = 10
      ),
      
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.margin = margin(
        t = 8,
        r = 10,
        b = 8,
        l = 8
      )
    )
  
  return(
    list(
      plot = p,
      qc = qc_summary,
      duplicates = duplicate_check
    )
  )
}


# ================================
# 7. Create all six plots
# ================================
plot_list <- vector(
  mode = "list",
  length = nrow(biomarker_info)
)

qc_list <- vector(
  mode = "list",
  length = nrow(biomarker_info)
)

duplicate_list <- vector(
  mode = "list",
  length = nrow(biomarker_info)
)

for (i in seq_len(nrow(biomarker_info))) {
  
  result <- make_boxplot(
    df = plot_data,
    
    biomarker_var =
      biomarker_info$var[i],
    
    biomarker_date_var =
      biomarker_info$date_var[i],
    
    biomarker_label =
      biomarker_info$short_label[i],
    
    ylab_text =
      biomarker_info$ylab[i]
  )
  
  plot_list[[i]] <- result$plot
  
  qc_list[[i]] <- result$qc
  
  duplicate_list[[i]] <-
    result$duplicates %>%
    mutate(
      biomarker =
        biomarker_info$short_label[i]
    )
}


# ================================
# 8. Combine into one panel
#
# Layout:
#   2 columns × 3 rows
# ================================
panel_plot <- wrap_plots(
  plot_list,
  ncol = 2
) +
  plot_annotation(
    title =
      "Early-Visit Plasma Biomarker Distributions",
    
    subtitle =
      "Sporadic PD NSD+ compared with Healthy Control NSD-; visit range: Early (BL-V02)",
    
    
    tag_levels = "A",
    
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 11,
        hjust = 0.5,
        margin = margin(b = 12)
      ),
      
      plot.caption = element_text(
        size = 9,
        hjust = 0,
        margin = margin(t = 12)
      ),
      
      plot.tag = element_text(
        size = 12,
        face = "bold"
      )
    )
  )


# Display panel in R
print(panel_plot)


# ================================
# 10. QC summaries
# ================================
qc_summary <- bind_rows(qc_list)

cat(
  "\nParticipant and missing-value QC:\n"
)

print(qc_summary)


# Participants with multiple eligible
# records for the same biomarker
duplicate_summary <- bind_rows(
  duplicate_list
)

cat(
  "\nParticipants with more than one eligible record:\n"
)

if (nrow(duplicate_summary) == 0) {
  
  cat(
    "No duplicate participant records detected.\n"
  )
  
} else {
  
  print(
    duplicate_summary %>%
      select(
        biomarker,
        PATNO,
        subgroup_definition,
        n_records
      ) %>%
      arrange(
        biomarker,
        subgroup_definition,
        PATNO
      )
  )
}