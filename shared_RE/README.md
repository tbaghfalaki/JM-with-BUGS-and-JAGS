
# Simulated data
The data set is generated based on a shared random effects. The longitudinal profiles and KM plot are shown below: 
![](/Figures/srm.png) 

# Difference between JAGS and BUGS in handling censoring
In OpenBUGS (specifically R2OpenBUGS), censoring can be handled using the commands $`\texttt{I (a, )}`$, $`\texttt{I ( , b)}`$, and $`\texttt{I (a, b)}`$ for right, left, and interval censoring, respectively. For example, when dealing with right censoring, we can use the following BUGS code 
  ```{r setup, include=FALSE}
 st_n[i] ~ dweib(kappa,mut[i])I(st_cen[i],)
 log(mut[i])<-inprod(alpha[],W[i,])+inprod(b[i,],gamma[])
```
where $`st\_n[i]`$ represents the actual survival times for uncensored or missing data (using NA) for censored data. Also, $`st\_cen[i]`$ takes on zero values when the survival times are observed and represent the censoring time for others.

On the other hand, JAGS uses the function $`\texttt{dinterval}(,)`$ to represent censored data as a distribution. The function takes two values, 0 and 1. If we write $`\texttt{is.censored}[i]\sim \texttt{dinterval}(t[i], c[i])`$, then for censored data $`\texttt{is.censored}[i]=1`$ and for the observed data $`\texttt{is.censored}[i]=0`$.
In survival analysis, it is also necessary to specify pairs $`t[i]`$ and $`c[i]`$. In the statement, $`t[i]`$ is either an actual time or a missing value (=NA), depending on whether the observation represents the time of the actual event or a censoring time, respectively. Also, $`c[i]`$ represents the censoring time. If the ith individual is not censored, we assign a very large number. The large number might represent the longest possible follow-up time in the study under analysis. In Listing following, the R syntax is provided for modifying the data to accommodate the distinctions between observed and censored data in JAGS.
```{r setup, include=FALSE}
is.censored <- 1 - death
t <- c <- rep(0, n)
for (i in 1:n) {
  if (death[i] == 0) ((t[i] <- NA) & (c[i] <- st[i]))
  if (death[i] == 1) ((t[i] <- st[i]) & (c[i] <- max(st)))
}
```
For more details about interval censoring or left censoring, refer to Rosner et al. (2021) and Alvares et al. (2021). After preparing the data, the BUGS code for this purpose is provided in Listing \ref{jagscode}.
```{r setup, include=FALSE}
t[i] ~ dweib(kappa,mut[i]) 
    log(mut[i])<-inprod(alpha[],W[i,])+inprod(b[i,],gamma[])
    is.censored[i]~dinterval(t[i],c[i])
```
As a comparison of the input for the R2jags and R2OpenBUGS, the following code shows the format of input required for the R2OpenBUGS:
```{r setup, include=FALSE}
st_n <- st_cen <- rep(0, n)
for (i in 1:n) {
  if (death[i] == 0) ((st_n[i] <- NA) & (st_cen[i] <- st[i]))
  if (death[i] == 1) ((st_n[i] <- st[i]) & (st_cen[i] <- 0))
}
```
Also, the following code shows the corresponding BUGS code in this package. 
```{r setup, include=FALSE}
st_n[i] ~ dweib(nu,mut[i])I(st_cen[i],)
log(mut[i])<-inprod(alpha[],W[i,])+inprod(b[i,],gamma[])
```

