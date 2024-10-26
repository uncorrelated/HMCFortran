\name{hmc.ologit}
\alias{hmc.ologit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
hmc.ologit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hmc.ologit(frml, beta.mu, beta.Sigma, gamma.mu, gamma.Sigma, data = NULL, N = 3000, BI = as.integer(N * 0.1), adjustSpan = BI, L = 20, epsilons = NULL, nchains = 2, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frml}{
%%     ~~Describe \code{frml} here~~
}
  \item{beta.mu}{
%%     ~~Describe \code{beta.mu} here~~
}
  \item{beta.Sigma}{
%%     ~~Describe \code{beta.Sigma} here~~
}
  \item{gamma.mu}{
%%     ~~Describe \code{gamma.mu} here~~
}
  \item{gamma.Sigma}{
%%     ~~Describe \code{gamma.Sigma} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{BI}{
%%     ~~Describe \code{BI} here~~
}
  \item{adjustSpan}{
%%     ~~Describe \code{adjustSpan} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{epsilons}{
%%     ~~Describe \code{epsilons} here~~
}
  \item{nchains}{
%%     ~~Describe \code{nchains} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
