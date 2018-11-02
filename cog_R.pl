#!/usr/bin/perl
use strict;

my $usage = <<USAGE;
Usage:
    perl $0 cog_cog_class_annot.xls out_prefix species_name

USAGE
if (@ARGV == 0){die $usage}

my $in=shift;
my $out=shift;
my $seqname=shift;

my $cog_R=<<COGR;
cog <- scan("$in",what="character",sep="\\n")
cog <- as.character(sapply(cog,function(x) strsplit(x,"\\t")[[1]][1:3]))
cog <- as.data.frame.list(tapply(cog,rep(1:3,t=(length(cog)/3)),as.character))
cog <- cog[-1,]
cog.des <- cog[,2]
cog.des <- gsub(" \$","",cog.des)

cog.des <- paste(as.character(cog[,1]),cog.des,sep=": ")
cog.num <- as.numeric(as.matrix(cog[,3]))
names(cog.num) <- as.character(cog[,1])

pdf(file="$out.pdf",10,7)
#plot.new()
def.par <- par()
layout(matrix(c(1,2),nr=1),widths=c(20,15))
op1 <- par(mar=c(4.1,5.1,4.1,0))
a <- barplot(cog.num,col=rainbow(24),ylim=c(0,max(cog.num)*1.1),cex.names=0.8,las=1,
axis.lty=0,cex.axis=0.6,xaxt="n"
)
box()
mtext("Function Class",side=1,line=2,font=2,cex=0.8)
mtext("Number of Contigs",side=2,line=2.5,font=2,cex=0.8)
par(cex.axis=0.6)
axis(side=1,at=a,labels=cog[,1],padj=-2,tick=FALSE)
op1 <- par(mar=c(0,0,1.1,1.1))
plot(0,0,pch="",xlim=c(0,1),ylim=c(0,26),bty="n",axes=F,xlab="",ylab="")
for(i in length(cog.des):1){
text(0.01,i,cog.des[length(cog.des)+1-i],pos=4,cex=0.7)
}
par(op1)
par(def.par)
title("COG Function Classification of $seqname\'s expressed contig sequences")

dev.off()

COGR

open COGRSH,'>',$out.'.R' or die "can't open the R file $out.R";
print COGRSH $cog_R;
close COGRSH;



my $cog_R_out=$out.'.Rout';
system("R CMD BATCH $out.R $cog_R_out ");
system("convert -density 98 $out.pdf $out.png ");
