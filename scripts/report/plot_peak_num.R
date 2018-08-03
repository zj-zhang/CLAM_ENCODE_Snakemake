## Compute the overlapping status of peak files
## and plot a heatmap
## Zijun Zhang
## 10.23.2017
## revised 1.8.2018: better visualization

library(ggplot2)

argv = commandArgs(trailingOnly=T)
print(argv)

clam_mpeak = argv[1]
clam_upeak = argv[2]

peak_idx = c('CLAM_m', 'CLAM_u')

peak_df = data.frame(
	x=c(rep('CLAM_m',2), rep('CLAM_u',2)),
	y=rep(c('CLAM_m','CLAM_u'),2),
	value=NA,
	stringsAsFactors=F
	)

for(i in 1:nrow(peak_df))
{
	p1 = which(peak_idx == peak_df$x[i])
	p2 = which(peak_idx == peak_df$y[i])
	cmd = paste0('bedtools intersect -a ', argv[p1], ' -b ', argv[p2], ' -u | wc -l')
	resp = system(cmd, intern=T)
	peak_df$value[i] = as.numeric(resp)
}


p = ggplot(data=peak_df, aes(x=x, y=y)) + 
	geom_tile(aes(fill=value)) + 
	scale_fill_gradient(low = "white",  high = "steelblue", guide=FALSE) +
	geom_text(aes(label = value), size=3) +
	xlab('Reference') + ylab('Target') +
	theme_bw() +
	theme(axis.text.x=element_text(size=9), axis.title=element_text(size=10), axis.text.y=element_text(size=9, angle = 90, hjust=0.5)) 
	
ggsave(argv[3], plot=p, width=2, height=2)