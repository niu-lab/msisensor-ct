args <- commandArgs(TRUE)
library(beeswarm)
library(boxplotdbl)
library(RColorBrewer,quiet=T)
m<-read.table(args[1])
nn<-ncol(m)
if(nn == 2){
    colnames(m)<-c( 'state', 'score')
    #colnames(m)<-c('sample','type','state', 'score')
    #colnames(m)<-c('state', 'score')
    n<-length(unique(m$state))
    cols <- colors<-colorRampPalette(brewer.pal(12, "Paired"))(n)
    xs<-1:n
    point_size<-0.4
    pdf(file=args[2])
    #jpeg(file=args[2])
    par.lend=par()$lend;par(tck=0.01, lend=1, xpd=T)
    beeswarm(score~state, data=m, method='swarm', pch=16, ylim=c(0, 1), xlab='Sample group', ylab='XGBoost probability value', main=args[1], cex=point_size, las=0, do.plot=T, add=F, at=xs, spacing=point_size, xaxs='i', yaxs='i', lwd=.75, col=c("blue","red"))
    legend("topright", legend=c("MSI","MSS"), col=c("blue","red"),pch=16)
    #legend("topright", legend=c("0:MSS","1:MSI"), col=c("blue","red"), pch=16)
    #boxplot(score~state, data = m, add = T)
    #abline( h = 3.5, lty = 2, col = "grey", lwd = 1)
    #abline( h = 8.5, lty = 2, col = "grey", lwd = 1, xpd=F)
    #beeswarm(score~state, data=m, method='swarm', pch=16, ylim=c(0, 100), xlab='Sample group', ylab='MSI score(%)', main='MSI score', cex=0.8, las=0, do.plot=T, add=F, at=xs, col=cols[xs], spacing=0.8, xaxs='i', yaxs='i', lwd=.75)
}else if(nn == 1){
    colnames(m)<-c('score')
    xs<-nrow(1)

    pdf(file=args[2])
    par.lend=par()$lend;par(tck=0.01, lend=1, xpd=T)
    #plot(c(1:xs), m$score, main='MSI score', xlab='Sample group', ylab="MSI score(%)", col="red", pch=19, spacing=1, xaxs='i', yaxs='i')
    beeswarm(m$score, data=m, method='swarm', pch=16, xlab='Sample group', ylab='MSI score(%)', main='MSI score distribution(tumor only)', cex=1, las=0, do.plot=T, add=F, at=xs, col=cols[1], spacing=1, xaxs='i', yaxs='i', lwd=.75)
}

xmin<-par('usr')[1];
xmax<-par('usr')[2];
ymin<-par('usr')[3];
ymax<-par('usr')[4];
chw<-par()$cxy[1]

axis(1,lwd=NA,tick=T,lwd.ticks=.5,labels=F)
axis(2,lwd=NA,tick=T,lwd.ticks=.5,labels=F)
axis(3,lwd=NA,tick=T,lwd.ticks=.5,labels=F)
axis(4,lwd=NA,tick=T,lwd.ticks=.5,labels=F)

