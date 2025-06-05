#CD4_TEM_subset_analysis
#subset out a cell type
CD4TEM <- subset(x=pbmc.tox.azi, subset= predicted.celltype.l2 == c("CD4 TEM"))
CD4TEM

#set Idents
Idents(CD4TEM) <- 'Population'
levels(CD4TEM)

#DGE CD4 cells
CD4.de.markers <- FindMarkers(CD4TEM, ident.1="Control", ident.2="Cohort")
head(CD4.de.markers)

#heatmap
data("CD4.de.markers")
plot3 <- DoHeatmap(object=CD4TEM, features=NULL, group.by="ident")
ggsave(plot3, filename="HeatMapCD4TEM2.png",width = 7.5, height = 10)

getwd()