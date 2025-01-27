% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxOR.R
\name{fitLogisticRegression}
\alias{fitLogisticRegression}
\title{Fit a logistic regression model and return the OR and 95\% confidence interval}
\usage{
fitLogisticRegression(data)
}
\arguments{
\item{data}{A data frame with columns: y, age_p}
}
\value{
A matrix of OR and 95\% confidence intervals
}
\description{
Fit a logistic regression model and return the OR and 95\% confidence interval
}
\examples{
# Generate simulated data dataC
dataC <- createData(200)

# Fit the semi-parametric model
spm.fit <- fitSemiParamModel(dataC)

# Find two cut-off points
dataC <- findCutoffs(spm.fit, dataC)

# Calculate the two cut-off points after limiting sensitivity se and specificity sp
result <- calculateCutoffs(dataC)
cutoffs <- result$cutoffs
dataC2 <- result$filteredData

# Discretize age variable based on the two cutoffs
dataC2 <- discretizeAge(dataC2, cutoffs)

# Fitting logistic regression models and obtaining OR values and 95\% confidence intervals
OR_Results <- fitLogisticRegression(dataC2)
print(round(OR_Results, 3))
}
