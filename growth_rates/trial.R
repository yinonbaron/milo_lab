raw.tmp.data = read.csv('/home/yinonbaron/Documents/Experiments/pyruvate transporter/20151220_pyr_trans_clones_different_cs.csv',sep=',')
cropped.data = raw.tmp.data[c(3:dim(raw.tmp.data)[1]),c(2:dim(raw.tmp.data)[2])]
