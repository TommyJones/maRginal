
# marginal <img src="figures/marginal_v1.svg" align="right" alt="logo" width="200"/>

marginal will be an R package for calculating marginal effects for arbitrary prediction models. As of now, code and research are in early stages. 

Open issues include:

* Finding a valid way to calculate confidence intervals
* Handling categorical variables more elegantly
* Finding the best default behaviors for calculations
* And more!

Some basic usage
```
devtools::install_github("tommyjones/marginal")

library(marginal)
library(randomForest)

data(mtcars)

fit <- randomForest(mpg ~., data = mtcars)

mfx <- CalcMfx(object = fit, X = mtcars)

plot(mfx)
```
