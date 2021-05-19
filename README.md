# nimble-virtual-2021

Materials for the virtual NIMBLE workshop, May 26-28 2021.

To prepare for the workshop:

 - Install NIMBLE (see below)
 - Install additional packages (see below and check back for updates)
 - (UNDER CONSTRUCTION) Join our Slack workspace (highly recommended but not required; see information in email)
 - (UNDER CONSTRUCTION) Download these materials (check back before the workshop on Wednesday for updates)
 - (UNDER CONSTRUCTION) See email for Zoom invitation.

All materials for the workshop will be in this GitHub repository. If you're familiar with Git/GitHub, you already know how to get all the materials on your computer. If you're not, simply click [here](https://github.com/nimble-training/nimble-virtual-2021/archive/master.zip).

There is some overview information [here (UNDER CONSTRUCTION)](https://htmlpreview.github.io/?https://github.com/nimble-training/nimble-virtual-2021/blob/master/overview.html), including links to the content modules in order.

Time: 8 am - 1 pm California time (GMT-8); 11 am - 4 pm Eastern US time.

Location: Zoom (see email) and [Slack (UNDER CONSTRUCTION)](https://2021nimbleworkshop.slack.com).

## Slack and Zoom instructions and tips/tricks

Please see [this document (UNDER CONSTRUCTION)](https://docs.google.com/document/d/1hhm6Eco0KevM30aDGdpo0n-gYDcT7IPbw7_vTrb-BVQ/edit?usp=sharing) for discussion of how to use Zoom and Slack during the workshop.

## Tentative Schedule

All times are California time. If we need to modify the start times of
units, we'll try to announce in advance on the Slack #general channel.

Day 1 (Wednesday May 26):

1. (8 am - 8:45 am) Introduction to NIMBLE: Basic concepts and workflows
2. (9 am - 10 am) Working with NIMBLE models and converting from WinBUGS/JAGS
3. (10:30 am - 11:45 pm) Comparing and customizing MCMC methods in NIMBLE
4. (12 pm - 1 pm) Strategies for improving MCMC

Day 2 (Thursday May 27):

5. (8 am - 9 am) Writing your own functions and distributions 
6. (9:15 am - 10:00 am) The basics of programming algorithms in NIMBLE
7. (10:30 pm - 1 pm with a 15 minute break) Spatial modeling

Day 3 (Friday May 28):

8. (8 am - 9:30 am) Model selection and Bayesian nonparametrics
9. (9:45 am - 1 pm with breaks) Special topics breakouts
   - (a) Ecological models
   - (b) Programming in NIMBLE, including writing your own MCMC sampler and calling out to R and C++

## Help with NIMBLE

Our user manual is [here](https://r-nimble.org/html_manual/cha-welcome-nimble.html).

We have a 'cheatsheet' [here](https://r-nimble.org/documentation).

## Installing NIMBLE

NIMBLE is an R package on CRAN, so in general it will be straightforward to install as with any R package, but you do need a compiler and related tools on your system.  

In summary, here are the steps.

1. Install compiler tools on your system. [https://r-nimble.org/download](https://r-nimble.org/download) has more details on how to install *Rtools* on Windows and how to install the command line tools of *Xcode* on a Mac. Note that if you have packages requiring a compiler (e.g., *Rcpp*) on your computer, you should already have the compiler tools installed.

2. Install the *nimble* package from CRAN in the usual fashion for an R package. More details (including troubleshooting tips) can also be found in Section 4 of the [NIMBLE manual](https://r-nimble.org/html_manual/cha-installing-nimble.html).

3) To test that things are working please run the following code in R:

```
library(nimble)
code <- nimbleCode({
  y ~ dnorm(0,1)
})
model <- nimbleModel(code)
cModel <- compileNimble(model)
```


If that runs without error, you're all set. If not, please see the troubleshooting tips and email nimble.stats@gmail.com directly if you can't get things going.  

In general we encourage you to update to the most recent version of NIMBLE, 0.11.0.


#### (Not required) Development version(s) of NIMBLE

Sometimes we make an update or new feature available on a github branch before it is released.  In the event a need arises to install from a branch, you can do so as follows (for branch "devel"):

```
library(remotes)
install_github("nimble-dev/nimble", ref = "devel", subdir = "packages/nimble")
```

## Installing additional packages

Some of the packages we will use (beyond those automatically installed with `nimble`) can be installed as follows:

```
install.packages(c("mcmcplots", "CARBayesdata", "sp", "spdep", "classInt"))
```

`compareMCMCs` is a package in development that is not yet on CRAN:

```
library(remotes)
install_github("nimble-dev/compareMCMCs", subdir = "compareMCMCs")
```

Windows users will probably need to use this invocation:

```
library(remotes)
install_github("nimble-dev/compareMCMCs", subdir = "compareMCMCs", INSTALL_opts = "--no-multiarch")
```


