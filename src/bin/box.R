library(dplyr)
library(ggplot2)

# Load the data
summary_df <- read.csv("Summary.csv")

# Create the plots
p1 <- ggplot(summary_df, aes(x = Group, y = RT_PedMeters, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "RT_PedMeters across Groups",
       y = "Cumulative PedMeters (first 4h)")

p2 <- ggplot(summary_df, aes(x = Group, y = RT_Water, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "RT_Water across Groups",
       y = "Cumulative Water (first 4h)")

# Save them to files
ggsave("RT_PedMeters_boxplot.jpg", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("RT_Water_boxplot.jpg", plot = p2, width = 8, height = 6, dpi = 300)


library(dplyr)
library(ggplot2)

# Load the data
summary_df <- read.csv("Summary.csv")

# --- F_WT vs F_S1RKO ---
f_data <- filter(summary_df, Group %in% c("F_WT", "F_S1RKO"))
wilcox.test(RT_PedMeters ~ Group, data = f_data)
wilcox.test(RT_Water ~ Group, data = f_data)

# --- M_WT vs M_S1RKO ---
m_data <- filter(summary_df, Group %in% c("M_WT", "M_S1RKO"))
wilcox.test(RT_PedMeters ~ Group, data = m_data)
wilcox.test(RT_Water ~ Group, data = m_data)

# Boxplots
p_f_ped <- ggplot(f_data, aes(x = Group, y = RT_PedMeters, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Female WT vs S1RKO: RT_PedMeters")

p_f_water <- ggplot(f_data, aes(x = Group, y = RT_Water, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Female WT vs S1RKO: RT_Water")

p_m_ped <- ggplot(m_data, aes(x = Group, y = RT_PedMeters, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Male WT vs S1RKO: RT_PedMeters")

p_m_water <- ggplot(m_data, aes(x = Group, y = RT_Water, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Male WT vs S1RKO: RT_Water")

# Save plots to files
ggsave("Female_PedMeters_boxplot.jpg", plot = p_f_ped, width = 8, height = 6, dpi = 300)
ggsave("Female_Water_boxplot.jpg", plot = p_f_water, width = 8, height = 6, dpi = 300)
ggsave("Male_PedMeters_boxplot.jpg", plot = p_m_ped, width = 8, height = 6, dpi = 300)
ggsave("Male_Water_boxplot.jpg", plot = p_m_water, width = 8, height = 6, dpi = 300)
