#' The MIT Cellphone Data
#'
#' The Massachusetts Institute of Technology (MIT) cellphone data consists of the human interactions via cellphone activity, among 96 participants for a duration of 232 days. The data were taken from 15-Sept-2004 to 4-May-2005 inclusive, which covers the winter and spring vacations in the MIT 2004-2005 academic calendar.
#'
#' @format A list with 232 matrices of dimension of 96 by 96.
#'
#' @references Eagle, Nathan. "Reality mining: sensing complex social systems." Personal and ubiquitous computing 10.4 (2006): 255-268.
#'
#' @examples
#' data("MITphone")
#' result <- CPD_STERGM(MITphone, directed=FALSE, network_stats=c("edges", "isolates", "triangles"))
#'
#' theta_change <- result$theta_change; threshold <- result$threshold
#' plot(1:length(theta_change), theta_change, type='l',ylab="", xlab="", xaxt="n")
#' abline(h = threshold, col='red',lwd=2)
#'
#' xtick <- result$est_CP; seq_date <- seq(as.Date("2004-09-15"), as.Date("2005-05-04"), by="days")
#' axis(side=1, at=xtick-2, labels = FALSE, lwd = 0, lwd.ticks = 1)
#' text(x=xtick-9,  par("usr")[3]-2.7, labels = seq_date[xtick], srt=45, cex=0.7, xpd=TRUE)
"MITphone"
