setwd("~/Documents/C.elegans_sync/Final_files/02-ANALYSES/")

###### Load the script with the functions I made ######
source("~/Documents/Scripts/R_FUNCTIONS.R")

##### Compute correlation between nb of syn and nsyn muts #####
## Load mutation table
MEGA <- read.table("R_txt_files/Mutation_List.txt",header=T,sep="\t")

#### Table for mutation per position ####
Count_Pos <- as.data.frame(table(MEGA$POS_REF))
colnames(Count_Pos) <- c("POS_REF","COUNT")
write.table(Count_Pos, file = "R_txt_files/Mut_per_pos.txt")
Count_Pos <- read.table("R_txt_files/Mut_per_pos.txt")

#### Plot ####
plot(Count_Pos$POS_REF, Count_Pos$COUNT, xlab = "Distance to Origin", ylab = "Number of Mutations", main = "Mutations vs. Distance")

lm_model <- lm(Count_Pos$COUNT ~ Count_Pos$POS_REF)
summary(lm_model)

#### windowed ####
{
  # Define your window size (adjust as needed)
  window_size <- 100  # Size of the window (e.g., 100 base pairs)
  
  # Create a new column for the window (floor the position divided by window size)
  Count_Pos$window <- floor(Count_Pos$POS_REF / window_size)
  
  # Aggregate the mutation counts by window and calculate the average distance to origin
  agg_df <- Count_Pos %>%
    group_by(window) %>%
    summarise(
      avg_mutations = mean(COUNT),  # Average number of mutations in the window
      avg_distance = mean(POS_REF)  # Average distance to origin for the window
    )
  
  # Load ggplot2 package for plotting
  library(ggplot2)
  
  # Function to plot the windowed approach regression
  plot_windowed_regression <- function(agg_df, lm_windowed) {
    # Create a scatter plot of average mutations vs. average distance
    ggplot(agg_df, aes(x = avg_distance, y = avg_mutations)) +
      geom_point(color = "blue", alpha = 0.7) +   # Scatter points
      geom_smooth(method = "lm", color = "red", se = FALSE) +  # Regression line
      labs(
        title = "Average Mutations vs. Average Distance (Windowed)",
        x = "Average Distance to Origin",
        y = "Average Number of Mutations"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12)
      )
  }
  
  # Plot the data
  plot_windowed_regression(agg_df, lm_windowed)
  
  # Fit a linear model to the aggregated data
  lm_windowed <- lm(avg_mutations ~ avg_distance, data = agg_df)
  
  # View the summary of the model
  summary(lm_windowed)
  
}