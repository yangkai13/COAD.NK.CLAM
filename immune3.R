######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

install.packages("RColorBrewer")

#读取两个文件：KM.txt 包含风险数据，Subtype_Immune_Model.txt 包含患者的临床亚型。
riskFile="KM.txt"                          
library(RColorBrewer)        
subtypeFile="Subtype_Immune_Model.txt"      
setwd("/Users/q")     

#ChischisqTest 函数用于对不同组的分布进行卡方检验，并返回 p 值。
ChischisqTest <- function(tabledata){
	ct = chisq.test(tabledata)
	ct.pvalue = ct$p.value
	ct.pvalue = ifelse(ct.pvalue<0.001,0.001,round(ct.pvalue,3))
}

#rect2text：一个辅助函数，用于在图形中绘制带有文本标签的矩形块，设置颜色和文本大小。
rect2text <- function(x1,y1,x2,y2,text,rect.col,text.cex,text.col="white",...){
  rect(x1,y1,x2,y2,col=rect.col,border=NA,...)
  text((x1+x2)/2,(y1+y2)/2,text,cex=text.cex,col=text.col,...)
}

#risk=read.table(riskFile,...)：读取预后风险数据
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#subtype=read.table(subtypeFile,...)：读取临床亚型数据，处理行名以便与预后风险数据匹配。
subtype=read.table(subtypeFile, header=T, sep="\t", check.names=F, row.names=1)
#rownames(subtype)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(subtype))

#sameSample=intersect(...)：提取两者共同的样本，保证数据一致性。
sameSample=intersect(row.names(subtype), row.names(risk))
subtype=subtype[sameSample,]
subtype=gsub(".+Immune |\\)","",subtype)
risk=risk[sameSample,]
data=cbind(subtype, as.data.frame(risk))
typeTab=table(data$subtype)    #typeTab=table(data$subtype)：计算每种临床亚型的患者数量。
typeName=names(typeTab[typeTab>6])
#merge.data=data[which(data[,"subtype"] %in% typeName),]。#merge.data=data[which(data[,"subtype"] %in% typeName),]：合并数据，将那些患者数大于6的亚型提取出来。
merge.data = data[which(data[,"subtype"] %in% typeName),]

cliName=colnames(merge.data)[1]
colnames(merge.data)[1]="Clinical"

#pergrouptable = table(merge.data$Clinical, merge.data$Risk)：生成临床亚型与风险的列联表。
#ct.pvalue = ChischisqTest(pergrouptable)：对亚型与预后风险进行卡方检验，计算 p 值。
groups = sort(unique(merge.data$Clinical))
groups.table = table(merge.data$Clinical)
groups.num = length(groups)
samples.num = nrow(merge.data)
groups.lownum = sum(merge.data$Risk=="L")
groups.highnum = sum(merge.data$Risk=="H")
pergrouptable = table(merge.data$Clinical, merge.data$Risk)
ct.pvalue = ChischisqTest(pergrouptable)

#设置了多个变量用于定义图形的坐标轴和颜色，比如 groups.col 是不同临床亚型的颜色，high_low.col 是高风险和低风险的颜色。
#调用 pdf(file="immSubtype.pdf",...)：将输出保存为 PDF 文件。
#通过 rect2text 函数绘制多个矩形块，表示不同的患者分组及其风险分布，并添加文本标签。
xlim1=0; xlim2=samples.num; ylim1=1; ylim2=10
space = dist.up = 0.1
space.x = 1/100*xlim2
groups.col = col = colorRampPalette(brewer.pal(9, "Paired"))(groups.num)
high_low.col = c(rgb(235/255,116/255,106/255),rgb(30/255,181/255,184/255))
bottom.unitwidth = xlim2 / (groups.num+1+1.3)
bottom.left.width = bottom.unitwidth*1.3


pdf(file="immSubtype.pdf", width=10, height=7)
plot(1,type="n",xlim=c(xlim1,xlim2),ylim=c(ylim1,ylim2),axes=F,xlab="",ylab="")

header.text.cex = 1.5
header.text = gettextf("%s TCGA patients",samples.num)
header.text.width = strwidth(header.text,cex=header.text.cex)
arrows.len = (xlim2 - header.text.width - space.x*2)/2
text(xlim2/2,9.5,header.text,cex=header.text.cex,font=2)
arrows(0,9.5,arrows.len,9.5,code=1,lwd=3,length=0.2)
arrows(xlim2-arrows.len,9.5,xlim2,9.5,code=2,lwd=3,length=0.2)
rect2text(x1=0,y1=4.5+space,x2=bottom.left.width,y2=7-space,text.cex=1.2,text = 
	"Riskscore\n groups",rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)
#draw bottom header
bottom.header.text = paste0(cliName, " group (n=" , samples.num, ")")
rect2text(x1=bottom.left.width+space.x,y1=6+space/2,x2=xlim2,y2=7-space,text.cex=1.2,
	text= bottom.header.text,rect.col=rgb(64/255,109/255,181/255),text.col="white",font=2)


#draw bottom left 2
bottom.left2.text = gettextf("Riskscore-low\n(n=%s)",groups.lownum)
rect2text(x1=0,y1=3+space/2,x2=bottom.left.width,y2=4.5-space/2,text.cex=1.2,text= 
	bottom.left2.text,rect.col=high_low.col[2],text.col="white",font=2)
#draw bottom left 3
bottom.left3.text = gettextf("Riskscore-high\n(n=%s)",groups.highnum)
rect2text(x1=0,y1=1.5+space/2,x2=bottom.left.width,y2=3-space/2,text.cex=1.2,text= 
	bottom.left3.text,rect.col=high_low.col[1],text.col="white",font=2)

#
start = 0
for(i in 1:length(groups.table)){
	end = groups.table[i]+start
	namesi = names(groups.table)[i]
	# up rect
	rect2text(x1=start,y1=8+space,x2=end,y2=9-space,text.cex=1.1,text= namesi, 
	  	rect.col=groups.col[i],text.col="white")
	merge.datai = merge.data[merge.data$Clinical==namesi,,drop=F]
	high.num = sum(merge.datai$Risk=="H")
	low.num = sum(merge.datai$Risk=="L")
	# middle low rect
	rect(start,7+dist.up,start+low.num,8-dist.up,col=high_low.col[2],border=NA)
	# middle high rect
	rect(start+low.num,7+dist.up,end,8-dist.up,col=high_low.col[1],border=NA)
	# bottom 1
	bottom.starti = bottom.left.width+(i-1)* bottom.unitwidth
	bottom.endi = bottom.starti+bottom.unitwidth
	subheader.text = gettextf("%s\n(n=%s, %s)",namesi,nrow(merge.datai),paste0(round(nrow(merge.datai)/samples.num*100),"%"))
	rect2text(x1=bottom.starti+space.x,y1=4.5+space,x2=bottom.endi,y2=6-space,text.cex=1.1,text= subheader.text,rect.col=groups.col[i],text.col="white")
	# bottom 2
	low.texti = gettextf("%s(%s)",low.num,paste0(round(low.num/groups.lownum*100),"%"))
	rect2text(x1=bottom.starti+space.x,y1=3+space/2,x2=bottom.endi,y2=4.5-space/2,text.cex=1.1,text= low.texti,rect.col="grey90",text.col="black")
	# bottom 3
	high.texti = gettextf("%s(%s)",high.num,paste0(round(high.num/groups.highnum*100),"%"))
	rect2text(x1=bottom.starti+space.x,y1=1.5+space/2,x2=bottom.endi,y2=3-space/2,text.cex=1.1,text= high.texti,rect.col="grey70",text.col="black")
	start = end
}

#????pvalue?ĳ?????
rect2text(x1=bottom.endi+space.x,y1=4.5+space,x2=xlim2,y2=6-space,text.cex=1.1,text= "P-value",rect.col="grey70",text.col="black",font=3)
rect2text(x1=bottom.endi+space.x,y1=1.5+space/2,x2=xlim2,y2=4.5-space/2,text.cex=1.1,text= ct.pvalue,rect.col="grey70",text.col="black")

#?????μ??ϱ߿?
rect(0,8+space,xlim2,9-space,lwd=2,border="grey30")
rect(0,7+space,xlim2,8-space,lwd=2,border="grey30")
rect(0,1.5+space/2,xlim2,7-space,lwd=2,border="grey30")

#????ͼ??
legend("bottom",legend=c("Riskscore-low","Riskscore-high"),col=rev(high_low.col),bty="n",ncol=2,pch=15,pt.cex=1.3,cex=1.3)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056


