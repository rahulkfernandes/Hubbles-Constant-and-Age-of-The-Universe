#### Part I:
# “I certify that I indeed finished reading Ch. 2 from 
# An Introduction to Statistical Learning, by James
# Gareth, Daniela Witten, Trevor Hastie, Robert Tibshirani”

library(tidyverse)

#### Part II:
### Q1
load("./galaxies.RData")

### Q2
View(galaxies)
str(galaxies)

# Check for missing values
sum(is.na(galaxies))

# Summary statistics
summary(galaxies)

attach(galaxies)

# Spread Statistics
range(velocity)
range(distance)
IQR(velocity)
IQR(distance)

var(velocity)
sd(velocity)

var(distance)
sd(distance)

# To check for any repeating galaxy names
unique(Galaxy)
duplicated(galaxies)

## Univariate Analysis
# Histogram for velocity
ggplot(galaxies, aes(x = velocity)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 4, alpha = 0.7) +  # Histogram
  geom_vline(aes(xintercept = mean(velocity, na.rm = TRUE)), 
             color = "black", linetype = "solid", linewidth = 1.2) +  # Solid black mean line
  geom_vline(aes(xintercept = median(velocity, na.rm = TRUE)), 
             color = "red", linetype = "dashed", linewidth = 1.2) +  # Dashed red median line
  labs(title = "Histogram of Velocity", x = "Velocity", y = "Count") +
  theme_minimal() +
  annotate("text", x = mean(velocity, na.rm = TRUE), y = 10, 
           label = "Mean", color = "black", angle = 90, vjust = -1, hjust=0.8) +
  annotate("text", x = median(velocity, na.rm = TRUE), y = 10, 
           label = "Median", color = "red", angle = 90, vjust = -1, hjust=0.8)


# Histogram for distance
ggplot(galaxies, aes(x = distance)) +
  geom_histogram(fill = "lavender", color = "black", bins = 5, alpha = 0.7) +  # Histogram
  geom_vline(aes(xintercept = mean(distance)), 
             color = "black", linetype = "solid", linewidth = 1.2) +  # Solid black mean line
  geom_vline(aes(xintercept = median(distance)), 
             color = "red", linetype = "dashed", linewidth = 1.2) +  # Dashed red median line
  labs(title = "Histogram of Distance", x = "Velocity", y = "Count") +
  theme_minimal() +
  annotate("text", x = mean(distance), y = 10, 
           label = "Mean", color = "black", angle = 90, vjust = -1, hjust=0.8) +
  annotate("text", x = median(distance), y = 10, 
           label = "Median", color = "red", angle = 90, vjust = -1, hjust=0.8)

# Density Plot for velocity
ggplot(galaxies , aes(x = velocity)) +
  geom_density(alpha = 0.6, fill='lightblue') +
  labs(title = "Density Plot of Velocity") +
  theme_minimal()

# Density Plot for distance
ggplot(galaxies , aes(x = distance)) +
  geom_density(alpha = 0.6, fill='lavender') +
  labs(title = "Density Plot of Distance") +
  theme_minimal()

# Boxplot for velocity
ggplot(galaxies, aes(y = velocity)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "Velocity Boxplot") +
  theme_minimal()

# Boxplot for distance
ggplot(galaxies, aes(y = distance)) +
  geom_boxplot(fill = "salmon") +
  labs(title = "Distance Boxplot") +
  theme_minimal()

#############
# Check for outliers in velocity
Q1 <- quantile(velocity, 0.25)
Q3 <- quantile(velocity, 0.75)
IQR_value <- IQR(velocity)

lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

outliers <- velocity[velocity < lower_bound | velocity > velocity]
print(paste("Outliers in velocity: ", outliers))

# Check for outliers in distance
Q1 <- quantile(distance, 0.25)  # First quartile (25th percentile)
Q3 <- quantile(distance, 0.75)  # Third quartile (75th percentile)
IQR_value <- IQR(distance)      # Compute IQR

lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

outliers <- distance[distance < lower_bound | distance > upper_bound]
print(paste("Outliers in distance: ", outliers))


## Bivariate Analysis
# Correlation matrix
numerical_data <- galaxies %>%
  select(velocity, distance)

cor_matrix <- cor(numerical_data, method = "pearson")

print(cor_matrix)

# Reshape matrix using base R
cor_matrix_melt <- as.data.frame.table(cor_matrix, responseName = "value")
colnames(cor_matrix_melt) <- c("Var1", "Var2", "value")

# Create heatmap with ggplot2
ggplot(cor_matrix_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limit = c(-1, 1)
  ) +
  labs(title = "Correlation Heatmap", x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Pearson correlation test
cor_test <- cor.test(distance, velocity, method = "pearson")
print(cor_test)

# Scatter plot with regression line
ggplot(galaxies, aes(x = distance, y = velocity)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Velocity vs. Distance", x = "Distance (Mpc)", y = "Velocity (km/s)") +
  theme_minimal()

### Q3
## Fitting Model
# No intercept linear model
model <- lm(velocity ~ distance + 0 , data = galaxies)

# View model summary
summary(model)

ggplot(galaxies, aes(x = distance, y = velocity)) +
  geom_point(alpha = 0.6) +
  geom_smooth(
    method = "lm",
    formula = y ~ x + 0,  # No-intercept line
    color = "red",
    se = TRUE  # Show confidence interval
  ) +
  labs(
    title = "Hubble's Law Fit: Velocity vs Distance",
    x = "Distance (Mpc)",
    y = "Velocity (km/s)"
  ) +
  theme_minimal()

### Q4
# Residuals vs. fitted values plot
ggplot(data = model, aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Residuals vs. Fitted Values",
    x = "Fitted Velocity",
    y = "Residuals"
  ) +
  theme_minimal()

# Residual distribution
ggplot(data = model, aes(x = .resid)) +
  geom_histogram(fill = "skyblue", bins = 5) +
  labs(title = "Distribution of Residuals") +
  theme_minimal()

# Q-Q plot for normality
qqnorm(residuals(model))
qqline(residuals(model), col = "red")

# Model fit quality
sigma(model)  # RSE (in km/s)
sqrt(mean(residuals(model)^2))  # RMSE (in km/s)
summary(model)$coefficients["distance", c("Estimate", "Std. Error", "Pr(>|t|)")]
summary(model)$r.squared # R-squared

### Q5
## Calculating Hubble's constant
hubble_constant <- coef(model)["distance"]
print(paste("Hubble Constant (H_0):", round(hubble_constant, 2), "km/s/Mpc"))

### Q6
## Approximating age of the universe
# Constants
mpc_to_km <- 3.086e19   # 1 Mpc = 3.086e19 km
seconds_per_year <- 3.154e7

# Calculate age in seconds
age_seconds <- (1 / hubble_constant) * mpc_to_km

# Convert to years
age_years <- age_seconds / seconds_per_year

# Print results
print(paste("Age of the Universe (approx):", 
            round(age_years / 1e9, 2), "billion years"))
