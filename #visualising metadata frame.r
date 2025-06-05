#visualising metadata frame 
head(pbmc.tox.azi@meta.data)
colnames(pbmc.tox.azi[[]])

#add data column
#timepoint - baseline, C2, C3

#add data column
#Cancer - RCC, melanoma, NSCLC, cSCC

#set Idents for Cell type
unique(pbmc.tox.azi$predicted.celltype.l2)
Idents(pbmc.tox.azi) <- 'predicted.celltype.l2'
levels(pbmc.tox.azi)RCC