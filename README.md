# Hubbles Constant and Age of The Universe
Rahul Kenneth Fernandes
02/13/2025



## Objectives:

1.  Read in the data from the object, which is saved as galaxies.RData.
2.  Perform EDA of the data.
3.  Fit a linear no-intercept model (called Hubble’s Law).
4.  Assess the quality of the model fit.
5.  Estimate Hubble’s constant, including units.
6.  Approximate the age of the universe.

## Setting Up the Environment


``` r
library(tidyverse)
```

## Load Data


``` r
load("./galaxies.RData")
```

## Exploratory Data Analysis

### Dataset


``` r
head(galaxies)
```

```
##     Galaxy velocity distance
## 1  NGC0300      133     2.00
## 2  NGC0925      664     9.16
## 3 NGC1326A     1794    16.14
## 4  NGC1365     1594    17.95
## 5  NGC1425     1473    21.88
## 6  NGC2403      278     3.22
```

``` r
glimpse(galaxies)
```

```
## Rows: 24
## Columns: 3
## $ Galaxy   <fct> NGC0300, NGC0925, NGC1326A, NGC1365, NGC1425, NGC2403, NGC2541, NGC20…
## $ velocity <int> 133, 664, 1794, 1594, 1473, 278, 714, 882, 80, 772, 642, 768, 609, 14…
## $ distance <dbl> 2.00, 9.16, 16.14, 17.95, 21.88, 3.22, 11.22, 11.75, 3.63, 13.80, 10.…
```

The dataset contains three columns: Galaxy, velocity and distance. Where, the Galaxy columns contains the catalog names of galaxies, the velocity column contains the velocity at which the respective galaxy is moving and the distance column contains the distance of the respective galaxy from Earth. The velocity is measured in km/s and the distance is measured in Mpc or megaparsec, where 1Mpc = 3.086e19 km.

### Summary Statistics


``` r
summary(galaxies)
```

```
##       Galaxy      velocity         distance    
##  IC4182  : 1   Min.   :  80.0   Min.   : 2.00  
##  NGC0300 : 1   1st Qu.: 616.5   1st Qu.: 8.53  
##  NGC0925 : 1   Median : 827.0   Median :13.08  
##  NGC1326A: 1   Mean   : 924.4   Mean   :12.05  
##  NGC1365 : 1   3rd Qu.:1423.2   3rd Qu.:15.87  
##  NGC1425 : 1   Max.   :1794.0   Max.   :21.98  
##  (Other) :18
```


``` r
attach(galaxies)
```

\newpage
### Univariate Analysis
In this part, we look at velocity and distance individually.

Range of Velocity:

``` r
range <-range(velocity)
print(paste("Min: ", range[1], "; Max: ", range[2]))
```

```
## [1] "Min:  80 ; Max:  1794"
```

Variance and Standard Deviation of Velocity:

``` r
print(paste("Variance:", var(velocity)))
```

```
## [1] "Variance: 262978.157608696"
```

``` r
print(paste("Standard Deviation:", sd(velocity)))
```

```
## [1] "Standard Deviation: 512.813960036869"
```
IQR of Velocity:

``` r
print(paste("IQR of velocity:", IQR(velocity)))
```

```
## [1] "IQR of velocity: 806.75"
```

Outliers in velocity (using 1.5 * IQR):

``` r
Q1 <- quantile(velocity, 0.25)
Q3 <- quantile(velocity, 0.75)
IQR_value <- IQR(velocity)

lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

outliers <- velocity[velocity < lower_bound | velocity > velocity]
print(outliers)
```

```
## integer(0)
```


``` r
ggplot(galaxies, aes(x = velocity)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 5, alpha = 0.7) +  # Histogram
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
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-9-1.png" alt="plot of chunk unnamed-chunk-9" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-9</p>
</div>


``` r
ggplot(galaxies , aes(x = velocity)) +
  geom_density(alpha = 0.6, fill='lightblue') +
  labs(title = "Density Plot of Velocity") +
  theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-10-1.png" alt="plot of chunk unnamed-chunk-10" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-10</p>
</div>


``` r
ggplot(galaxies, aes(y = velocity)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "Velocity Boxplot") +
  theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-11-1.png" alt="plot of chunk unnamed-chunk-11" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-11</p>
</div>

Velocity: The median velocity (827.0 km/s) is lesser than the mean velocity (924.4 km/s), hence indicating a right skew in the data. This can be observed in the histogram of velocity. This can also be observed in the box plot, where the median line is not at the center of the box plot. The box plot also shows no outliers. Considering the minimum velocity (80.0 km/s), maximum velocity (1794.0 km/s), the IQR (806.75 km/s) and the density plot, we can observe that the data is not closely packed near the center and has a wide spread. We can also see this when we calculate the Standard Deviation of the velocity data (512.814 km/s). Hence, there is a deviation from the normal expected bell curve of a normal distribution.


Range of Distance:

``` r
range <- range(distance)
print(paste("Min: ", range[1], "; Max: ", range[2]))
```

```
## [1] "Min:  2 ; Max:  21.98"
```
Variance and Standard Deviation of Distance:

``` r
print(paste("Variance:", var(distance)))
```

```
## [1] "Variance: 33.8101389492754"
```

``` r
print(paste("Standard Deviation:", sd(distance)))
```

```
## [1] "Standard Deviation: 5.81464865226398"
```
IQR of Distance:

``` r
print(paste("IQR of distance:", IQR(distance)))
```

```
## [1] "IQR of distance: 7.34"
```

Outliers in distance (using 1.5 * IQR):

``` r
Q1 <- quantile(distance, 0.25)
Q3 <- quantile(distance, 0.75)
IQR_value <- IQR(distance)

lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

outliers <- distance[distance < lower_bound | distance > upper_bound]
print(outliers)
```

```
## numeric(0)
```


``` r
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
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-16-1.png" alt="plot of chunk unnamed-chunk-16" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-16</p>
</div>


``` r
ggplot(galaxies , aes(x = distance)) +
  geom_density(alpha = 0.6, fill='lavender') +
  labs(title = "Density Plot of Distance") +
  theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-17-1.png" alt="plot of chunk unnamed-chunk-17" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-17</p>
</div>


``` r
ggplot(galaxies, aes(y = distance)) +
  geom_boxplot(fill = "lavender") +
  labs(title = "Distance Boxplot") +
  theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-18-1.png" alt="plot of chunk unnamed-chunk-18" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-18</p>
</div>

Distance: the median distance (13.08 Mpc) is greater than the mean distance (12.05 Mpc), hence indicating a left skew in the data. This can be observed in the histogram of distances. This can also be observed in the box plot, where the median line is not at the center of the box plot. The box plot also shows no outliers. Considering the minimum distance (2.00 Mpc), maximum distance (21.98 Mpc), the IQR (7.34 Mpc) and the density plot, we can observe that the data is closely packed near the center and has a wide spread. We can also see this when we calculate the Standard Deviation of the distances (5.814649 Mpc). Hence, there is a deviation from the normal expected bell curve of a normal distribution.


### Bivariate Analysis
Correlation matrix:


``` r
numerical_data <- galaxies %>%
  select(velocity, distance)

cor_matrix <- cor(numerical_data, method = "pearson")
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
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-19-1.png" alt="plot of chunk unnamed-chunk-19" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-19</p>
</div>

Pearson Correlation Test:

``` r
cor_test <- cor.test(distance, velocity, method = "pearson")
print(cor_test)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  distance and velocity
## t = 8.0189, df = 22, p-value = 5.677e-08
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.7054127 0.9394555
## sample estimates:
##       cor 
## 0.8631815
```

Linear Relationship:


``` r
ggplot(galaxies, aes(x = distance, y = velocity)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Velocity vs. Distance", x = "Distance (Mpc)", y = "Velocity (km/s)") +
  theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-21-1.png" alt="plot of chunk unnamed-chunk-21" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-21</p>
</div>
The correlation matrix heat map and the Pearson Correlation test show us that there is high correlation between velocity and distance with a Pearson's Coefficient of 0.86. The linear relationship is confirmed by the scatter plot of velocity and distance with a regression line plotted through it. From this information, we can proceed to fit a linear model.

## Fitting a No-Intercept Model
No-Intercept Linear Model (also called as Hubble's Law):

``` r
model <- lm(velocity ~ distance + 0, data = galaxies)
```

Model Summary:

``` r
summary(model)
```

```
## 
## Call:
## lm(formula = velocity ~ distance + 0, data = galaxies)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -736.5 -132.5  -19.0  172.2  558.0 
## 
## Coefficients:
##          Estimate Std. Error t value Pr(>|t|)    
## distance   76.581      3.965   19.32 1.03e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 258.9 on 23 degrees of freedom
## Multiple R-squared:  0.9419,	Adjusted R-squared:  0.9394 
## F-statistic: 373.1 on 1 and 23 DF,  p-value: 1.032e-15
```

Plotting the Linear Model:


``` r
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
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-24-1.png" alt="plot of chunk unnamed-chunk-24" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-24</p>
</div>

## Assessing the quality of the model fit
Residuals vs. fitted values plot:


``` r
ggplot(data = model, aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Residuals vs. Fitted Values",
    x = "Fitted Velocity",
    y = "Residuals"
  ) +
  theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-25-1.png" alt="plot of chunk unnamed-chunk-25" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-25</p>
</div>
From Residuals vs. fitted values plot, we can see that the points are randomly distributed around 0 without any pattern. This shows that the linear model holds good and the relationship between velocity and distance is linear.

Residual distribution:


``` r
ggplot(data = model, aes(x = .resid)) +
  geom_histogram(fill = "skyblue", bins = 5) +
  labs(title = "Distribution of Residuals") +
  theme_minimal()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-26-1.png" alt="plot of chunk unnamed-chunk-26" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-26</p>
</div>

Q-Q plot for normality:


``` r
qqnorm(residuals(model))
qqline(residuals(model), col = "red")
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-27-1.png" alt="plot of chunk unnamed-chunk-27" width="80%" />
<p class="caption">plot of chunk unnamed-chunk-27</p>
</div>
In the Q-Q plot, majority of the points fall near the red line, indicating that the residuals are approximately normal.The points at the left and right tail, deviate from the red line which indicates the presence of heavy tails in the distributions. This can also be observed from the residual histogram which follows an approximate bell shaped curve.  

Model fit quality:

```
## [1] "RSE (in km/s): 258.933064035632"
```

```
## [1] "RMSE (in km/s)) 253.481231058161"
```

```
## [1] "R-squared:  0.941931023545281"
```
The RSE (258.93) and RSME (253.48) are close to each other, which suggests homoscedasticity, i.e, residuals have constant variance. The model typical prediction error is around 253-259 km/s. A high R-squared (0.9419) shows that the model fits the data very well.

Distance Coefficient Significance:

``` r
summary(model)$coefficients["distance", c("Estimate", "Std. Error", "Pr(>|t|)")]
```

```
##     Estimate   Std. Error     Pr(>|t|) 
## 7.658117e+01 3.964794e+00 1.031907e-15
```
The Estimate of the coefficient is calculated as 76.58. The standard error is small, indicating a more precise estimation. Since p < 0.05, "distance" is a statistically significant predictor. The near-zero p value shows that distance has a meaningful effect on velocity.

## Calculating Hubble's Constant
We can calculate Hubble's Constant by calculating the distance coefficient from the no-intercept linear model. 


``` r
hubble_constant <- coef(model)["distance"]
print(paste("Hubble Constant (H_0):", round(hubble_constant, 2), "km/s/Mpc"))
```

```
## [1] "Hubble Constant (H_0): 76.58 km/s/Mpc"
```

## Approximating the Age of the Universe

``` r
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
```

```
## [1] "Age of the Universe (approx): 12.78 billion years"
```

