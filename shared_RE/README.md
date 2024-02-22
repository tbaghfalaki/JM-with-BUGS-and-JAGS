
# Simulated data
The data set is generated based on a shared random effects. The longitudinal profiles and KM plot are shown below: 
![](/Figures/srm.png) 

# Difference between JAGS and BUGS in handling censoring
In OpenBUGS (specifically R2OpenBUGS), censoring can be handled using the commands $`\texttt{I (a, )}`$, $`\texttt{I ( , b)}`$, and $`\texttt{I (a, b)}`$ for right, left, and interval censoring, respectively. For example, when dealing with right censoring, we can use the following BUGS code 
  ```{r setup, include=FALSE}
 st_n[i] ~ dweib(kappa,mut[i])I(st_cen[i],)
 log(mut[i])<-inprod(alpha[],W[i,])+inprod(b[i,],gamma[])
```
where $st\_n[i]$ represents the actual survival times for uncensored or missing data (using NA) for censored data. Also, $st\_cen[i]$ takes on zero values when the survival times are observed and represent the censoring time for others.\\ 
On the other hand, JAGS uses the function $\texttt{dinterval}(,)$ to represent censored data as a distribution. The function takes two values, 0 and 1. If we write $\texttt{is.censored}[i]\sim \texttt{dinterval}(t[i], c[i])$, then for censored data $\texttt{is.censored}[i]=1$ and for the observed data $\texttt{is.censored}[i]=0$.
In survival analysis, it is also necessary to specify pairs $t[i]$ and $c[i]$. In the statement, $t[i]$ is either an actual time or a missing value (=NA), depending on whether the observation represents the time of the actual event or a censoring time, respectively. Also, $c[i]$ represents the censoring time. If the $i$th individual is not censored, we assign a very large number. The large number might represent the longest possible follow-up time in the study under analysis\cite{rosner2021bayesian}. In Listing \ref{censjags}, the R syntax is provided for modifying the data to accommodate the distinctions between observed and censored data in JAGS.

\begin{lstlisting}[caption=The right censored data in R syntax for JAGS.,label=censjags]
is.censored <- 1 - death
t <- c <- rep(0, n)
for (i in 1:n) {
  if (death[i] == 0) ((t[i] <- NA) & (c[i] <- st[i]))
  if (death[i] == 1) ((t[i] <- st[i]) & (c[i] <- max(st)))
}
\end{lstlisting}
For more details about interval censoring or left censoring, refer to Rosner et al.\cite{rosner2021bayesian} and Alvares et al.\cite{alvares2021bayesian}. After preparing the data, the BUGS code for this purpose is provided in Listing \ref{jagscode}.
\begin{lstlisting}[caption=The BUGS code for the AFT model for JAGS.,label=jagscode]
t[i] ~ dweib(kappa,mut[i]) 
    log(mut[i])<-inprod(alpha[],W[i,])+inprod(b[i,],gamma[])
    is.censored[i]~dinterval(t[i],c[i])
\end{lstlisting}
As a comparison of the input for the \texttt{R2jags} and \texttt{R2OpenBUGS}, Listing \ref{cenopen} shows the format of input required for the \texttt{R2OpenBUGS}. Also, Listing \ref{codeopen} shows the corresponding BUGS code in this package. It should be mentioned that, overall, \texttt{R2jags} has a shorter computational time than \texttt{R2OpenBUGS}. 
All the R codes including a set of generated data under the assumption of the shared random effects model can be accessed at \url{https://github.com/tbaghfalaki/JM-with-BUGS-and-JAGS/tree/main/shared_RE}. \\
This section involves simulating real data and analysing them from a single Gaussian longitudinal marker and time-to-event outcome. To see an illustrative example of such a dataset, please refer to \cite{alvares2021bayesian}. In the following sections, we will expand on joint modeling and provide real datasets for each extension.
