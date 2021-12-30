\name{bayes.capacity}
\alias{bayes.capacity}
\title{
A Bayesian Model of Capacity Across Trials
}
\description{
Uses a Bayesian updating method to calculate the Capacity Coefficient for First-Terminating (OR) Processing across trials (time).
}
\usage{
bayes.capacity(n=1000, s, x, plist, a=20, ss, sm=40, fun=sq.exp, cap=OR, chans=c(CH12))
}
\arguments{
  \item{n}{The number of Monte Carlo samples. Default = 1000}
  \item{s}{Step size. This assumes equal intervals between $Step.}
  \item{x}{The data frame to estimate over. This should already be subset to the subject/block, contain only correct data, and have columns $Step, $RT, $Channel1, $Channel2}
  \item{plist}{A list of parameter values to use for the prior of each single Channel and their combination. The list should be structured such that plist[[1]][1] = shape Channel 1, plist[[1]][2] = scale Channel 1, etc.  NOTE: The scale value used in this model is defined differently than the base R function provides. You must transform the scale BEFORE running this model. Use the scale_bayes() function should you fit and create the params list yourself. Alternatively, you may use the Weibull.paramfit() function to estimate these using your data; the result will have shape and scale values already tranformed for using this model.}
  \item{a}{A weight value, instantiated as number of samples, given to the prior. Default = 20}
  \item{ss}{A vector of the dicrete steps in time (i.e., trial number, time stamps) to estimate capacity.}
  \item{sm}{The drop-off rate by which data from previous steps influence the current estimate. This is a parameter used for creating the weighted smoothing function. Default = 40.}
  \item{fun}{The weighting function to smooth over the data from t=0 up to the current estimate. Default = sq.exp (squared exponential)}
  \item{cap}{The stopping rule for the capacity estimate. Default = OR}
  \item{chans}{The channels (CH1, CH2, CH12) assumed to be nonstationary. Default = CH12}
}
\details{
The capacity coefificent models efficient performance based on assumptions about how humans may process multiple channels of information. There are multiple metrics for measuring processing efficiency using a UCIP model comparison depending on the construct of interest. For instance, if the UCIP model responds as a first-terminating decision process about multiple targets (OR stopping rule) and \eqn{H(t)} is the cumulative hazard function, the capacity coefficient is the ratio of observed performance when all channels are present, \eqn{H_1...n}, to the UCIP model prediction, 
\deqn{C_{\rm OR}(t)=\frac{H_{\rm OR}(t)}{\sum_i H_i(t)}.}{Cor(t)=[Hor(t)]/sum[Hi(t)],} 
where \eqn{n} is the total number of channels. A nonparametric estimator of \eqn{C_{OR}(t)} can be calculated using the capacity.OR() function, and the appropriate interpretation of results are explained in existing documentation. 

A parameteric estimator of the capacity coefficient uses a Bayesian approach to achieve high temporal specificity. A Weibull distribution with two parameters, shape, \beta, and scale, \theta, can be used to describe response times such that the hazard function of the Weibull is given by: 
\deqn{h_{\rm t|\beta,\theta}=\frac{\beta}{\theta}\left[\frac{t}{\theta}\right]^{\beta-1}}{h(t|beta,theta)=[beta/theta][t/theta]^(beta-1),}

and the cumulative hazard function by:
\deqn{H_{\rm t|\beta,\theta}=\left(\frac{t}{\theta}\right)^\beta}{H(t|beta,theta)=(t/theta)^beta}

Assuming both single channels and their combination are nonstationary, we can estimate the capacity coefficient by estimating the shape and scale parameters for multi-channel performance \eqn{(beta_{1...n}, theta_{1...n})} and performance with each single channel \eqn{(beta_{i}, theta_{i})}, respectively. The parametric OR capacity coefficient assuming a Weibull distribution is, 

\deqn{C(t) = \frac{H \left(t|\beta_{1...n}, \theta{1...n}\right)} {\sum_{i=1}^{n}H\left(t|\beta_{i}, \theta_{i}\right)} = \frac{\left(\frac{t}{\theta_{1...n}\right)^\beta_{1...n}}}{\sum_{i=1}^{n}\left(\frac{t}{\theta_{i}}\right)^\beta_{i}}.}{C(t) = H(t/beta_{1...n}, theta_{1...n})/sum_{i=1}^n(H(t|beta_{i}, theta_{i})) = (t/theta_{1...n})^{beta_{1...n}}/sum_{i=1}^{n}(t/theta_{i})^beta_{i}.}

To estaimte a trial(time)-varying capacity coefficient, we can rely on Bayesian updating using a conjugate prior of the Weibull function, which simplies the posterior estimation. The Weibull function with a known shape, \eqn{\beta} belongs to the exponential family, therefore a conjugate prior exists. The conjugate prior of the \eqn{\theta} parameter of the Weibull function with a known \beta is an inverse-gamma distribution, \pi, with two parameters, a and b. 

\deqn{\pi\left(\theta \mid a,b\right) = \frac{b^{a-1}e^{\frac{-b}{\theta}}}{\Gamma\left(a-1\right)\theta^{a}}.}{pi(theta|a,b)=(b^{a-1}e{(-b)/theta})/(Gamma(a-1)theta^{a}).}

The posterior distribution of the scale (\theta) parameter, given a fixed shape (\beta), hyperparameters a and b and m observations x_{i}...x_{m}, is given by,

\deqn{\left(\theta \mid X\right) = \pi\bigleft( \theta \mid a + m, b + \sum_{i=1}^{m} x^{\beta}_{i} \bigright),}{P(theta|X) = pi(theta|a + m, b + sum_{i=1}^{m}(x_{i}^{beta})),}

and the mode of the maximum a posteriori estimate of the scale can easily be calculated by,

\deqn{\theta_{\mid X} = \frac{b+\sum_{x \in X}X^{\beta}}{a+m+1}.}{theta_{|X} = (b+sum_{x in X}X^beta)/(a+m+1).}

Nonstationary performance in trials containing information from a single channel, multiple channels, or both could result in a varying (or stable given the same rate of nonstationarity in each term) capacity estimate across time. One may choose to update estimates of either or both the UCIP model and multichannel performance after the observation of each respective term used to calculate the capacity coefficient. 

A simplifying assumption is that there is a single shape parameter across conditions and throughout the experiment. Therefore, capacity on a given trial, s, is given by,

\deqn{C\left(t,s\right) = \frac{\theta^{-\beta}_{s,1...n}}{\sum_{i=1}^{n}{\theta_{s,i}^{-\beta}}}.}{C(t,s) = (theta_{s,1...n}^{-beta})/(sum_{i=1}^{n}(theta^{-beta}_{s,i})).}

If the prior a is high relative to m, then the posterior will change less relative to the prior. A parameter w_{j}, for each X_{j} in X, can explicitly set the relative value of each new observed data point given by,

\deqn{P\left(\theta \mid X\right) = \pi\bigleft(\theta \mid a+ \sum_{j=1}^{m}{w_{j}}, b+\sum_{j=1}^{m}{X^{\beta}_{j}w_{j}}\bigright).}{P(theta|X)=pi(theta|a+sum_{j=1}^{m}(w_{j}), b+sum_{j=1}^{m}(X^{beta}_{j}w_{j})).}

There are numerous potential options for the form of the weights. For example, a squared exponential function: if a response time, \eqn{X_{j}} was observed at time \eqn{S_{j}}, then the weight given at trial S is, 

\deqn{w_{j}(t) = exp \bigleft[ \frac{-\left(S_{j}-s\right)^2}{d}\bigright],}{w_{j}(t) = exp [-(S_{j}-S)^2/d],}

where d is the smoothing parameter (sm) for the distribution of weights across trials (time). 


}
\value{
  \item{hdi}{A list with [[1]] = the mode capacity estimate at each step, [[2]] = the corresponding upper (97.5\%) HDI and [[3]] = the lower (2.5\%) HDI at each step, respectively.}
}

\references{
Townsend, J.T. & Nozawa, G. (1995). Spatio-temporal properties of elementary perception: An investigation of parallel, serial and coactive theories. \emph{Journal of Mathematical Psychology, 39}, 321-360.

Houpt, J.W. & Townsend, J.T. (2012). Statistical Measures for Workload Capacity Analysis.  \emph{Journal of Mathematical Psychology, 56}, 341-355.

Houpt, J.W., Blaha, L.M., McIntire, J.P., Havig, P.R. and Townsend, J.T. (2013). Systems Factorial Technology with R. \emph{Behavior Research Methods, 46}, 307-330.

Fox, E.L. and Houpt, J.W. (2021). A Bayesian model of capacity across trials. \emph{Journal of Mathematical Psychology, 105}, 102604.
}
\author{
Betsy Fox <elizabeth.fox.9@us.af.mil>
}
\seealso{
\code{\link{ucip.test}},
\code{\link{capacity.or}}
}

\examples{
#fixed shape (a), prior scale
shape = 1.6; scale=1.6

up = c(shape, 1.5) #unlimited
lp = c(shape, 2) #limited
sp = c(shape, 1) #super

#unlimited capacity 
u_ch1b = abs(rnorm(1, mean = up[2], sd = .1))
u_ch2b = abs(rnorm(1, mean = up[2], sd = .1))

#lower capacity --- slow downs in processing each source
l_ch1b = abs(rnorm(1, mean = lp[2], sd = .1)) 
l_ch2b = abs(rnorm(1, mean = lp[2], sd = .1))

#higher capacity --- speed ups in processing each source
s_ch1b = abs(rnorm(1, mean = sp[2], sd = .1))
s_ch2b = abs(rnorm(1, mean = sp[2], sd = .1))

nsamps = 100


#limited capacity
l_Channel1 = rweibull(nsamps, shape=shape, scale=u_ch1b)
l_Channel2 = rweibull(nsamps, shape=shape, scale=u_ch2b)
l_Channel12 = pmin(rweibull(nsamps, shape=shape, scale=l_ch1b), rweibull(nsamps, shape=shape, scale=l_ch2b))
l.dat = data.frame(cap=I(rep("limited", length.out=nsamps*3)), Step=I(rep(1:nsamps, 3)), Channel1=I(c(rep(1, length.out=nsamps), rep(0, length.out=nsamps),rep(1, length.out=nsamps))), Channel2=I(c(rep(0, length.out=nsamps), rep(1, length.out=nsamps),rep(1, length.out=nsamps))), RT = I(c(l_Channel1, l_Channel2, l_Channel12)))
l.rows = sample(nrow(l.dat)); l.dat = l.dat[l.rows,] # shuffle the order of Channel1, Channel2, Channel12 in sequence
l.dat['Step'] = seq(1, nsamps*3)


#unlimited capacity
u_Channel1 = rweibull(nsamps, shape=shape, scale=u_ch1b)
u_Channel2 = rweibull(nsamps, shape=shape, scale=u_ch2b)
u_Channel12 = pmin(rweibull(nsamps, shape=shape, scale=u_ch1b), rweibull(nsamps, shape=shape, scale=u_ch2b))
u.dat = data.frame(cap=I(rep("unlimited", length.out=nsamps*3)), Step=I(rep(1:nsamps, 3)), Channel1=I(c(rep(1, length.out=nsamps), rep(0, length.out=nsamps),rep(1, length.out=nsamps))), Channel2=I(c(rep(0, length.out=nsamps), rep(1, length.out=nsamps),rep(1, length.out=nsamps))), RT = I(c(u_Channel1, u_Channel2, u_Channel12)))
u.rows = sample(nrow(u.dat)); u.dat = u.dat[u.rows,] # shuffle the order of Channel1, Channel2, Channel12 in sequence
u.dat['Step'] =  seq(nsamps*3+1, nsamps*3*2)

#super capacity
s_Channel1 = rweibull(nsamps, shape=shape, scale=u_ch1b)
s_Channel2 = rweibull(nsamps, shape=shape, scale=u_ch2b)
s_Channel12 = pmin(rweibull(nsamps, shape=shape, scale=s_ch1b), rweibull(nsamps, shape=shape, scale=s_ch2b))
s.dat = data.frame(cap=I(rep("super", length.out=nsamps*3)), Step=I(rep(1:nsamps, 3)), Channel1=I(c(rep(1, length.out=nsamps), rep(0, length.out=nsamps),rep(1, length.out=nsamps))), Channel2=I(c(rep(0, length.out=nsamps), rep(1, length.out=nsamps),rep(1, length.out=nsamps))), RT = I(c(s_Channel1, s_Channel2, s_Channel12)))
s.rows = sample(nrow(s.dat)); s.dat = s.dat[s.rows,] # shuffle the order of Channel1, Channel2, Channel12 in sequence
s.dat['Step'] =  seq(nsamps*3*2+1, nsamps*3*3)

dat = rbind(l.dat, u.dat, s.dat)
dat$RT = as.numeric(dat$RT)

plist = Weibull.paramfit(cap='capacity.or', sh=shape, sc=scale, ns=nsamps)

res = bayes.capacity(n=1000, s=1, x=dat, plist=plist, a=40, ss=dat$Step, sm=20, fun='sq.exp', cap='OR', chans=CH12)

windows()
plot(0:1, 0:1, type = 'n', xlim = c(0, max(dat$Step)), ylim=c(0,max(res[[2]])), axes=FALSE, ann=FALSE)
lines(res[[1]], lwd=2, col='black')
axis(side=1, at=seq(from=0, to=length(res[[1]]), length.out=20), labels=round(seq(from=1, to=length(res[[1]]), length.out=20),0))
axis(side=2, at=seq(from=0, to=3.5, by=.5), labels=round(seq(from=0, to=3.5, by=.5),2), cex=.8)
title(ylab="Capacity (Ratio of Hazards)(s)", xlab="Sample #")
title(main="Example of Bayes Ct Across Trials", cex.main=.9)
lines(res[[2]], lty=2, col='black'); lines(res[[3]], lty=2, col='black')
abline(h=1, lty=1) #unlimited cap
abline(v=c(nsamps*3, nsamps*6), lty=3)
mtext(c('limited','unlimited','super'),side=1, line=-2.5, at=c(125, 440, 755))
legend('topleft', y=c(1.65,2.0), legend=c('Bayes C(t) Across Time', '95\% HDI', 'Unlimited Capacity'), col='black', lty=c(1,2,1), lwd=c(2,1,1), cex=.8)
}
\keyword{ sft }