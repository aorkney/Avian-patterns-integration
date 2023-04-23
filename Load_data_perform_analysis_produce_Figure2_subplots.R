# This script will perform the essential Analyses in 
# Orkney et al., 2021 'patterns of integration...' : https://doi.org/10.1038/s41559-021-01509-w
# It will produce a facsimile of Figure 2 from this manuscript
# The code has been updated for R version 4.2.3
# This has resulted in subtle changes to some test statistics, but no major changes
# to the overall results. 
# The results demonstrate that the organisation of evolutionary integration between bone sizes is strongly modular.
# 'Head', 'Body', 'Wing' and 'Leg' modules are clearly defined. 
# However, the evolutionary integration between bone shapes does not exhibit such clear modular organisation, which
# defied expectation. 
# The shapes and sizes of bird bones exhibit different patterns of strength of ecological adaptation, 
# and it is possible that ecological adaptatation explains why shape and size exhibit different organisations.


# Our analysis will require the following packages:

require( geomorph )
require( stringdist )
require( ape )
require( abind )
require( nlme )
require( ggplot2 )
require( reshape2 )
require( ggnewscale )
require( cluster )
require( ggrepel )
require( vegan )

# If you have not yet installed these packages, run 'install.packages('insert package name here')'


load("D:/Documents/Alex_birds/NEW_analysis/2021_Bjarnason_125_SIdata/Supp 3 Datasets for Bjarnason & Benson 2021 revised submission 26 January 2021/Bird_landmarks_clean_26January2021_Bjarnason&Benson taxon set.RData")
# Load landmark constellation dataset 
# You will need to download this from 
# You will need to change the file pathway in the 'load' command appropriately. 
# This dataset can be sourced from:
# Bjarnason and Benson (2021) https://doi.org/10.18563/journal.m3.125
# The file names inside this dataset may be slightly different to those I have used here

landmark.info<-read.csv("D:\\Documents\\Alex_birds\\Bird_landmark_info_23Nov2019.csv" , row.names = 1)
# Load landmark constellation metadata
# This file should be inside the dataset files from Bjarnason and Benson (2021)
# https://doi.org/10.18563/journal.m3.125
# The file may have a slightly different name. 
# If you cannot find the file, I can upload it to my GitHub repository at user request. 

setwd("D:/Documents/Alex_birds/")
# I am changing the work directory. You may or may not need to change this line, depending on
# how you store your files on your computer. 

PrumTree <- read.nexus(file = "Avian-TimeTree.tre")
# Load the Avian Phylogeny of Prum et al., (2015)
# You can download this phylogeny here: 
# Avian Phylogeny: Prum et al., (2015) https://doi.org/10.1038/nature15697

setwd('D:/Documents/Alex_birds/NEW_analysis/')
# I am changing the working directory again, because I need to import another spreadsheet, 
# which is stored in a different location on my computer. 
metadata <- read.csv( "Eco_meta_data.csv" , row.names = 1 )
# Load metadata (this contains information about bird masses and ecological properties)
# You can find this spreadsheet here: https://github.com/aorkney/Bird-bodymass-2023/blob/main/Eco_meta_data.csv

SL.counts <- "min"
# The minimum number of slinding landmarks will be used 
# 'Sliding landmarks' are a type of landmark used to approximate continuous curves in 3-d Geometries. 
# Bjarnason and Benson published several versions of their data, with different numbers of sliding landmarks. 

elements <- c( "skull" , "mandible" , "scapula", "coracoid", "sternum" , "humerus", "ulna", "radius", "carpometacarpus", "synsacrum" , "femur", "tibiotarsus", "tarsometatarsus" )
# This is a vector of all skeletal elements / bones that were studied


# The following lines will produce a table of all bird species in the dataset
	name_list<-list()
	for(i in 1:length(elements)){
		name_list[[i]]<-dimnames(compiled.bird.landmarks.list[[ elements[ i ] ]][[ SL.counts ]][[ "output.landmarks.array" ]] )[[ 3 ]]
	}
		AA<-table(unlist(name_list))
			birds <- names( AA )[ AA == length( elements ) ]
				birds <- birds[ birds != "Falco_sparverius_UMMZ_154452" ]
				# !Idiosyncracy: I have removed a duplicate entry
				# You may not need to do this; if your code stops working at this point, comment out this command.

# The following lines will perform a generalised procrustes alignment on the landmark constellations. 
	GPA.results <- list()
		for( i in 1:length( elements ) )	{
			element <- elements [ i ]
				GPA.results[[ i ]] <- gpagen(approxBE=T, compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "output.landmarks.array" ]][ , , birds ] ,  curves = compiled.bird.landmarks.list[[ element ]][[ SL.counts ]][[ "sliders" ]]  )
			}
	
names( GPA.results ) <- elements
# Set the names of the aligned landmark constellations to the skeletal unit names


# !Idiosyncracy: The following lines will correct spelling mistakes in the original data
# The data from Bjarnason and Benson 2021 may have been corrected already: 
# You should check that the original names correspond to the same genus as the new names
birds[7]<-'Anseranas semipalmata NHMUK 1852.7.22.1'
birds[10]<-'Archilochus_colubris_FMNH_484738'
birds[12]<-'Choriotis_kori_UMMZ_210444'
birds[13]<-'Arenaria_interpres_FMNH_313992'
birds[30]<-'Pipra_erythrocephala_UMMZ_157210'
birds[31]<-'Chaetura_brachyrua_UMMZ_157689'
birds[32]<-'Charadris_vociferus_FMNH_470173'
birds[37]<-'Larus_novae-hollandiae_UMZC_274.C_'
birds[44]<-'Columbina_minuta_FMNH_289160'
birds[48]<-'Chloropipio_holochlora_FMNH_288169'
birds[66]<-'Podica_senegalensis_UMZC_209A_'
birds[75]<-'Leptotila_rufazilla_FMNH_318659'
birds[84]<-'Myiobius_varbatus_FMNH_386812'
birds[105]<-'Poecile_atricapillis_FMNH_504323'
birds[107]<-'Probosciger_aterriumus_NHMUK_S2006.15.5'
birds[117]<-'Rhamphastos_ambiguus_NHMUK_S2002.21'
birds[124]<-'Rupicola_peruviana_UMMZ_119210'


# We have performed the rotation to align the skeletons. We have rejected birds that do not appear in the Meta data 
# All of the rotated constellations host 149 birds

# The following lines organise the landmark and centroid size data
get.item <- function( X , item ) { X[[ item ]] }	
GPA.coords <- lapply( GPA.results , get.item , item = "coords" )
GPA.Csize <- lapply( GPA.results , get.item , item = "Csize" )


# The following lines prune the original phylogeny to the taxa of interest
tree_genera <- as.character( metadata[ birds , "Prum_genus" ] )
tree_names <- PrumTree$tip.label[ match( tree_genera , sub( "_.*" , "\\1" , PrumTree$tip.label ) ) ]
pruned.tree <- drop.tip( PrumTree , PrumTree$tip.label[ !PrumTree$tip.label %in% tree_names ] )

# The following lines produce a vector of body masses for the taxa of interest
masses <- apply( metadata[ birds , c( "Mass_F_.hbw_alive." , "Mass_M_.hbw_alive." ) ] , 1 , mean )
names( masses ) <- tree_names

# Now we can begin analysis
# This loop will remove shape variation associated with allometric scaling
# Size is quantified as the logarithm of the mass of the bird species
# We will use a Procrustes GLS regression to attain residuals. It is phylogenetically gnostic.
# That means that statistical non-independence caused by relatedness is controlled for.
# A Brownian model of evolution is assumed (this is the only currently available method). 

allometry_list <- list()
allometry.residuals <- list()

for(i in 1:length( GPA.coords ) ) {
	shape.temp <- GPA.coords[[ i ]]
		dimnames( shape.temp )[[ 3 ]] <- tree_names
			gdf.temp <- geomorph.data.frame( shape = shape.temp , size = log10( masses ) , phy = pruned.tree )
	allometry_list[[ i ]] <- procD.pgls( shape ~ size , phy = pruned.tree , data = gdf.temp ) # This is the line that performs the GLS regression
	names(allometry_list)[[i]]<-elements[i]
		for( j in 1:nrow( allometry_list[[ i ]]$pgls.residuals ) ) {
			if( j == 1 ) {
				allometry.residuals[[ i ]] <- matrix( allometry_list[[ elements[i] ]]$pgls.residuals[ j , ] , ncol = 3 , byrow = T ) + matrix (t(GPA.results[[ elements [ i ] ]]$consensus), ncol = 3, byrow = T)
					rownames( allometry.residuals[[ i ]] ) <- dimnames( shape.temp )[[ 1 ]]
					colnames( allometry.residuals[[ i ]] ) <- c( "x" , "y" , "z" )
					} else {
			 	allometry.residuals[[ i ]] <- abind( allometry.residuals[[ i ]] , matrix( allometry_list[[ elements[i] ]]$pgls.residuals[ j , ] , ncol = 3 , byrow = T ) + matrix( t(GPA.results[[ elements [ i ] ]]$consensus), ncol = 3, byrow = T), along = 3 )
		}	}	
					dimnames( allometry.residuals[[ i ]] )[[ 3 ]] <- dimnames( shape.temp )[[ 3 ]]
	}
						names( allometry.residuals ) <- names( allometry_list ) <- names( GPA.coords )



# If you want to view the residual skull shapes, then remove the '#' marks for the following 
# lines and run them. 
# (red ones are the originals and pink are the residuals. Differences are subtle. 
#for(i in 1:149){
#	plot3d(GPA.coords[['skull']][,,i],size=15,col='red'); aspect3d("iso")
#	Sys.sleep(.5)
#	plot3d(allometry.residuals[['skull']][,,i],size=15,col='hotpink'); aspect3d("iso");i<-i+1
#	Sys.sleep(.5)
#}

# The following lines turn the 3 dimensional array into a 2 dimensional array
# with dimensions (n , p*3)

allometry.residuals.matrices<-list()
	for( i in 1:length( allometry.residuals ) ) {
		element <- elements[ i ]
			n <- dim( allometry.residuals[[ i ]] )[3] # number of birds 
			p <- dim( allometry.residuals[[ i ]] )[1] # number of landmarks
		temporary_matrix <- matrix( NaN , n , p*3 )
			for(j in 1:n) {
				temporary_matrix[j,] <- c(t(allometry.residuals[[ elements[i] ]][,,j]))
			}
				colnames(temporary_matrix) <- c( t( matrix( c(paste("x",dimnames(allometry.residuals[[ elements[i] ]])[[1]],sep="_"),
					paste("y",dimnames(allometry.residuals[[ elements [i ] ]])[[1]],sep="_"),
					paste("z",dimnames(allometry.residuals[[ elements [i ] ]])[[1]],sep="_")) , ncol = 3 ) ) )
						rownames(temporary_matrix)<-tree_names
						
				allometry.residuals.matrices[[i]]<-temporary_matrix
	}

		names( allometry.residuals.matrices ) <- elements

# The residual shape data is all in a tidy list of (n, p*3) matrices now
# 1 matrix represents 1 bone


# We now want to look at size variation of bones; variation in skeletal proportions may tell us a different story
# to variation in skeletal shape
# We will quantify size as 'centroid size' of the landmark constellations
# We will need to remove the effects of variance in bird mass
# Use gls() from the package nlme and the "correlaiton = corPagel" argument
# from the package ape
# This is a phylogenetically gnostic way of getting at 'allometry-adjusted' centroid size

	allometry.Csize <- list()
		for( i in 1:length( GPA.Csize ) ) {
			df <- data.frame( mass = log10( masses ) , Csize = log10( GPA.Csize[[ elements[ i ] ]] ) )
			species<- rownames(df)
				allometry.Csize[[ i ]] <- gls( Csize ~ mass , correlation = corPagel( 0 , phy = pruned.tree, form= ~species   ) , data = df )
		}

		names( allometry.Csize ) <- names( GPA.Csize )
			lapply( allometry.Csize , summary )

		residual_Csize <- lapply( lapply( allometry.Csize , get.item , item = "residuals" ) , c )


# ecological data will now be prepared for analysis
# Please check that the objects 'food', 'flight' and 'hindlimb' possess column names that make sense;
# you may need to change the indices (the numbers in the square brackets) that occur after 'metadata'

food<-metadata[birds,38:57]
rownames(food)<-tree_names

food[,"Beak_use_Andrew"]<-as.numeric(as.factor(food[,"Beak_use_Andrew"]))
food$Diet.5Cat<-as.numeric(as.factor(food$Diet.5Cat))

flight<-metadata[birds,58:78]
rownames(flight)<-tree_names
flight<-flight[,-which(complete.cases(colSums(flight))==F)]

hindlimb<-metadata[birds,c(14,16,17,18,19,20,21,22,23,25,26,27,28,29,30,31,32,33,34,36)]
rownames(hindlimb)<-tree_names


# These lines will transform the ecological data, which is of mixed data types, 
# into a continuous multivariate representation.
# The method used is 'Principal Coordinate Analysis', and this is a non-linear 
# transformation. 
# Only those axes explaining more than 5% of variance will be retained. 
# It is necessary to perform this transformation because the subsequent statistical analysis we will apply assumes data is continuous.

fooddist <- daisy(food[,-12],metric="gower",stand=T,type=list(asymm = 19))
	foodpcoa <- cmdscale(fooddist,eig=T)
k<-max(which(foodpcoa$eig/sum(foodpcoa$eig)>=0.05))
foodpcoa <- cmdscale( fooddist, k=k  ,eig = T)

flightdist <- daisy(flight,metric="gower",stand=F,type=list(asymm=c(1:8)))
	flightdist[which(is.na(flightdist)==T)]<-1
	flightpcoa <- cmdscale(flightdist,eig=T)
k<-max(which(flightpcoa$eig/sum(flightpcoa$eig)>=0.05))
flightpcoa <- cmdscale( flightdist, k=k  ,eig = T)

hindlimbdist <- daisy(hindlimb,metric="gower",stand=F,type=list(asymm=c(1,2,6,7,8,9,10,11,12,13,14,15,16,18,20)))
	hindlimbdist[which(is.na(hindlimbdist)==T)]<-1
	hindlimbpcoa <- cmdscale(hindlimbdist,eig=T)
k<-max(which(hindlimbpcoa$eig/sum(hindlimbpcoa$eig)>=0.05))
hindlimbpcoa <- cmdscale( hindlimbdist, k=k  ,eig = T)

# The ecological metadata will now be added to the list of bone shape data, for subsequent analysis.

allometry.residuals.matrices[[14]]<-foodpcoa$points
allometry.residuals.matrices[[15]]<-flightpcoa$points
allometry.residuals.matrices[[16]]<-hindlimbpcoa$points
names(allometry.residuals.matrices)[14:16]<-c("food","flight","foot-use")
bones<-c(elements,"food","flight","foot-use")

# Now we will look at the patterns of evolutionary covariance between different sets of the bone-shape objects we have prepared
# We will use phylo.integration, so that we can do this in a phylogenetically gnostic way

phylo_store_allometry_removed<-list()
bone_combinations<-combn(c("skull","mandible","scapula","coracoid","sternum","humerus",
"ulna","radius","carpometacarpus","synsacrum","femur","tibiotarsus","tarsometatarsus","food","flight","foot-use"),2,simplify=T)
for(i in 1:max(dim(bone_combinations)) ){
	x<-paste(bone_combinations[,i][1])
	y<-paste(bone_combinations[,i][2])
	possibleError<-tryCatch(
	object<-phylo.integration(as.matrix(allometry.residuals.matrices[[paste(x)]]),as.matrix(allometry.residuals.matrices[[paste(y)]]),phy=pruned.tree,iter=999) # This is the line performing the evolutionary integration test
	,error=function(e)e)
	if(inherits(possibleError, "error")) next
	else {
		phylo_store_allometry_removed[[i]]<-object
	} 
	rm(object)
	print(i)
}
errors<-list()
for(i in 1:length(phylo_store_allometry_removed)){
	if(class(phylo_store_allometry_removed[[i]])!="pls"){
		errors[[i]]<-1
	} else {
		errors[[i]]<-0
	}
}
errors<-unlist(errors)

# We are now going to extract the effect sizes and p-values of the covariance statistics 
compar<-compare.pls(phylo_store_allometry_removed[which(errors==0)])

P_values_shape<-list()
for(i in 1:length(errors)){
P_values_shape[[i]]<-summary(phylo_store_allometry_removed[[i]])$P.value
}
P_values_shape<-unlist(P_values_shape)

effect_sizes<-list()
number_of_errors<-0
for(i in 1:length(errors)){
	if(errors[i]==0){
		effect_sizes[[i]]<-compar$sample.z[i-number_of_errors]
	} else {
		effect_sizes[[i]]<-NaN
		number_of_errors<-number_of_errors+1
	}
}
effect_sizes<-unlist(effect_sizes)

# the resultant test statistics will now be formatted into matrices

image<-matrix(NaN,16,16)
rownames(image)<-(c("skull","mandible","carpometacarpus","radius","ulna","humerus","sternum","coracoid","scapula","synsacrum","femur","tibiotarsus","tarsometatarsus","food","flight","foot-use"))
colnames(image)<-rownames(image)
for(i in 1:length(bones)){
	for(j in 1:length(bones)){
		column<-which(bone_combinations[1,]==colnames(image)[i])
		row<-which(bone_combinations[2,column]==rownames(image)[j])
		if(length(row)>0){
			image[i,j]<-as.numeric(effect_sizes[column[row]])
		}
	}
}

Pimage<-matrix(NaN,16,16)
rownames(Pimage)<-(c("skull","mandible","carpometacarpus","radius","ulna","humerus","sternum","coracoid","scapula","synsacrum","femur","tibiotarsus","tarsometatarsus","food","flight","foot-use"))
colnames(Pimage)<-rownames(Pimage)
for(i in 1:length(bones)){
	for(j in 1:length(bones)){
		column<-which(bone_combinations[1,]==colnames(Pimage)[i])
		row<-which(bone_combinations[2,column]==rownames(Pimage)[j])
		if(length(row)>0){
			Pimage[i,j]<-as.numeric(P_values_shape[column[row]])
		}
	}
}


# Now we will repeat the investigate the pattern of covariance, in a phylogenetic framework, for the skeletal proportion objects
# These are based on centroid ssize, rather than 3-D landmark constellation shape.

residual_Csize[[14]]<-foodpcoa$points
residual_Csize[[15]]<-flightpcoa$points
residual_Csize[[16]]<-hindlimbpcoa$points
names(residual_Csize)[14:16]<-c("food","flight","hindlimb")
phylo_size<-list()
bone_combinations<-combn(c("skull","mandible","scapula","coracoid","sternum","humerus",
"ulna","radius","carpometacarpus","synsacrum","femur","tibiotarsus","tarsometatarsus","food","flight","foot-use"),2,simplify=T)
for(i in 1:max(dim(bone_combinations)) ){
	x<-paste(bone_combinations[,i][1])
	y<-paste(bone_combinations[,i][2])
	possibleError<-tryCatch(
	object<-phylo.integration(as.matrix(residual_Csize[[which(bones==x)]]),as.matrix(residual_Csize[[which(bones==y)]]),phy=pruned.tree,iter=999) # This is the line performing the evolutionary integration test
	,error=function(e)e)
	if(inherits(possibleError, "error")) next
	else {
		phylo_size[[i]]<-object
	} 
	rm(object)
	print(i)
}

errors<-list()
for(i in 1:length(phylo_size)){
	if(class(phylo_size[[i]])!="pls"){
		errors[[i]]<-1
	} else {
		errors[[i]]<-0
	}
}
errors<-unlist(errors)

# Extract the effect sizes of the covariances
compar<-compare.pls(phylo_size[which(errors==0)])

P_values_size<-list()
for(i in 1:length(errors)){
P_values_size[[i]]<-summary(phylo_size[[i]])$P.value
}
P_values_size<-unlist(P_values_size)


effect_sizes<-list()
number_of_errors<-0
for(i in 1:length(errors)){
	if(errors[i]==0){
		effect_sizes[[i]]<-compar$sample.z[i-number_of_errors]
	} else {
		effect_sizes[[i]]<-NaN
		number_of_errors<-number_of_errors+1
	}
}
effect_sizes<-unlist(effect_sizes)

# Store the results in matrices

imageb<-matrix(NaN,16,16)
rownames(imageb)<-(c("skull","mandible","carpometacarpus","radius","ulna","humerus","sternum","coracoid","scapula","synsacrum","femur","tibiotarsus","tarsometatarsus","food","flight","foot-use"))
colnames(imageb)<-rownames(imageb)
for(i in 1:length(bones)){
	for(j in 1:length(bones)){
		column<-which(bone_combinations[1,]==colnames(imageb)[i])
		row<-which(bone_combinations[2,column]==rownames(imageb)[j])
		if(length(row)>0){
			imageb[i,j]<-as.numeric(effect_sizes[column[row]])
		}
	}
}


Pimageb<-matrix(NaN,16,16)
rownames(Pimageb)<-(c("skull","mandible","carpometacarpus","radius","ulna","humerus","sternum","coracoid","scapula","synsacrum","femur","tibiotarsus","tarsometatarsus","food","flight","foot-use"))
colnames(Pimageb)<-rownames(Pimageb)
for(i in 1:length(bones)){
	for(j in 1:length(bones)){
		column<-which(bone_combinations[1,]==colnames(Pimageb)[i])
		row<-which(bone_combinations[2,column]==rownames(Pimageb)[j])
		if(length(row)>0){
			Pimageb[i,j]<-as.numeric(P_values_size[column[row]])
		}
	}
}

# The following lines will manipulate the test statistic matrices, to make them ready for plotting.
# I was quite inexperienced at writing code when I originally wrote this script, and have not updated this section.
# It is probably much longer than it needs to be. 
# I am not going to commend this section. 

temporary<-image
temporary[is.na(image)]<- -0
temporary2<-t(image)
temporary2[is.na(t(image))]<- -0

temporary3<-temporary+temporary2
diag(temporary3)<-max(temporary3)
temporary3[which(temporary3==0)]<-mean(temporary3[which(temporary3!=0)])
diag(temporary3)<-max(effect_sizes[complete.cases(effect_sizes)])

Ptemporary<-Pimage
Ptemporary[is.na(Pimage)]<- -0
Ptemporary2<-t(Pimage)
Ptemporary2[is.na(t(Pimage))] <- -0
Ptemporary3<-Ptemporary+Ptemporary2
diag(Ptemporary3)<-max(Ptemporary3)
Ptemporary3[which(Ptemporary3==0)]<-mean(Ptemporary3[which(Ptemporary3!=0)])
diag(Ptemporary3)<-min(P_values_shape[complete.cases(P_values_shape)])


temporary4<-apply(t(temporary3), 2, rev)
rownames(temporary4)<-rev(rownames(temporary3))

Ptemporary4<-apply(t(Ptemporary3), 2, rev)
rownames(Ptemporary4)<-rev(rownames(Ptemporary3))
Ptemporary5<-Ptemporary4


temporary5<-temporary4
temporary5[which(Ptemporary5>0.05)]<-NA


temporaryb<-imageb
temporaryb[is.na(imageb)]<- -0
temporary2b<-t(imageb)
temporary2b[is.na(t(imageb))]<- -0

temporary3b<-temporaryb+temporary2b
diag(temporary3b)<-max(temporary3b)
temporary3b[which(temporary3b==0)]<-mean(temporary3b[which(temporary3b!=0)])
diag(temporary3b)<-max(effect_sizes[complete.cases(effect_sizes)])


Ptemporaryb<-Pimageb
Ptemporaryb[is.na(Pimageb)]<- -0
Ptemporary2b<-t(Pimageb)
Ptemporary2b[is.na(t(Pimageb))] <- -0
Ptemporary3b<-Ptemporaryb+Ptemporary2b
diag(Ptemporary3b)<-max(Ptemporary3b)
Ptemporary3b[which(Ptemporary3b==0)]<-mean(Ptemporary3b[which(Ptemporary3b!=0)])
diag(Ptemporary3b)<-min(P_values_size[complete.cases(P_values_size)])

Ptemporary4b<-apply(t(Ptemporary3b), 2, rev)
rownames(Ptemporary4b)<-rev(rownames(Ptemporary3b))
Ptemporary5b<-Ptemporary4b

temporary4b<-apply(t(temporary3b), 2, rev)
rownames(temporary4b)<-rev(rownames(temporary3b))

temporary5b<-temporary4b
temporary5b[which(Ptemporary5b>0.05)]<-NA

# Now the test statistics are almost read to plot

# create colour palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col.o <- colorRampPalette(c('white','orange','darkorange2','darkorange4','darkred','black'))(20) #'darkorange2',,'darkorange4'
col.b <- colorRampPalette(c('white','turquoise1','turquoise4','blue','darkblue','black'))(20)
bone_colours<-c(rep(cbbPalette[7],2),rep(cbbPalette[4],4),rep(cbbPalette[6],4),rep(cbbPalette[8],3),rep(cbbPalette[1],3))


# Plot the results: 

Image_data<-as.matrix(t(temporary5[1:16,1:16]))
Image_data<-Image_data[,-c(1:3)]

longData<-melt(t(Image_data))
longData$Diag<-as.character(longData$Var1) == as.character(longData$Var2)
longData$Diag[!longData$Diag] <- NA

# 'mp1' is a plot of pairwise integration for skeletal 3-d shapes

mp1<-
ggplot(longData[which(is.na(longData$Diag==T)),], aes(x = Var2, y = Var1)) + 
	coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  geom_tile(aes(fill=value)) + 
   scale_fill_gradientn(colours = col.o[3:20],na.value="white") +
	labs(fill = expression(italic('Z')))+
 new_scale("fill") +
  geom_tile(data= longData, aes(fill=Diag)) + 
    scale_fill_manual(guide = FALSE,values = c('TRUE' = "lightgrey"), na.value="transparent")+
  labs(x="", y="") +
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, vjust=0.3,colour=bone_colours,face = "bold"),
plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
legend.title=element_text(size=20),
                     axis.text.y=element_text(size=10,colour=rev(bone_colours)[4:16],face = "bold"),
                     plot.title=element_text(size=11,hjust=.5), legend.position="right",panel.background = element_blank())+
	labs(fill = "Effect Size")+
geom_hline(yintercept = 3.5,size=1,colour="black")+
geom_hline(yintercept = 11.5,size=1,colour="black")+
geom_vline(xintercept = 2.5,size=1,colour="black")+
geom_vline(xintercept = 6.5,size=1,colour="black")+
geom_hline(yintercept = 7.5,size=1,colour="black")+
geom_vline(xintercept = 13.5,size=1,colour="black")+
geom_vline(xintercept = 10.5,size=1,colour="black")+
scale_x_discrete(position = "top") 


Image_data<-as.matrix(t(temporary5b[1:16,1:16]))
Image_data<-Image_data[,-c(1:3)]
longData<-melt(t(Image_data))
longData$Diag<-as.character(longData$Var1) == as.character(longData$Var2)
longData$Diag[!longData$Diag] <- NA

# 'mp2' is a plot of pairwise integration for skeletal centroid sizes

mp2<-
ggplot(longData[which(is.na(longData$Diag==T)),], aes(x = Var2, y = Var1)) + 
	coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  geom_tile(aes(fill=value)) + 
   scale_fill_gradientn(colours = col.b[3:20],na.value="white") +
	labs(fill = expression(italic('Z')))+
 new_scale("fill") +
  geom_tile(data= longData, aes(fill=Diag)) + 
    scale_fill_manual(guide = FALSE,values = c('TRUE' = "lightgrey"), na.value="transparent")+
  labs(x="", y="") +
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, vjust=0.3,colour=bone_colours,face = "bold"),
plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
legend.title=element_text(size=20),
                     axis.text.y=element_text(size=10,colour=rev(bone_colours)[4:16],face = "bold"),
                     plot.title=element_text(size=11,hjust=.5), legend.position="right",panel.background = element_blank())+
	labs(fill = "Effect Size")+
geom_hline(yintercept = 3.5,size=1,colour="black")+
geom_hline(yintercept = 11.5,size=1,colour="black")+
geom_vline(xintercept = 2.5,size=1,colour="black")+
geom_vline(xintercept = 6.5,size=1,colour="black")+
geom_hline(yintercept = 7.5,size=1,colour="black")+
geom_vline(xintercept = 13.5,size=1,colour="black")+
geom_vline(xintercept = 10.5,size=1,colour="black")+
scale_x_discrete(position = "top") 

sd1<-(temporary5-mean(temporary5,na.rm=T))/sd(temporary5,na.rm=T)
sd1[is.na(temporary5)]<-0
sd2<-(temporary5b-mean(temporary5b,na.rm=T))/sd(temporary5b,na.rm=T)
sd2[is.na(temporary5b)]<-0
object<-sd2-sd1

Image_data<-object
Image_data<-Image_data[-c(1:3),]

longData<-melt((Image_data))
longData$Diag<-as.character(longData$Var1) == as.character(longData$Var2)
longData$Diag[!longData$Diag] <- NA

first_values<-seq(0,0.5,length.out=11)
adjustment1<-seq(0,0.3969191-0.5,length.out=11)
first_values<-first_values+adjustment1

second_values<-seq(0.5,1,length.out=11)
adjustment<-seq(0.3969191-0.5,0,length.out=11)
second_values<-second_values+adjustment

values=rev(c(first_values,second_values[2:11]))

values=seq(10,-5,length.out=10)

col.com<-c(rev(col.o),col.b)

z<-max(abs(object))

# 'mp3' is a plot that visualises the difference between the pairwise integraiton of skeletal shape and size

mp3<-
ggplot(longData[which(is.na(longData$Diag==T)),], aes(x = Var2, y = Var1)) + 
	coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")+
  geom_tile(aes(fill=value)) + 
   scale_fill_gradientn( colours=(col.com[15:27]), na.value = "white",limits=c(-2.0,2.0))+ 
	#scale_fill_gradient2(low=col.t[15],high=col.t[7],na.value='white')+
	labs(fill = expression(Delta))+
 new_scale("fill") +
  geom_tile(data= longData, aes(fill=Diag)) + 
    scale_fill_manual(guide = FALSE,values = c('TRUE' = "lightgrey"), na.value="transparent")+
  labs(x="", y="") +
  theme_bw() + theme(axis.text.x=element_text(size=10, angle=90, vjust=0.3,colour=bone_colours,face = "bold"),
plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
legend.title=element_text(size=20),
                     axis.text.y=element_text(size=10,colour=rev(bone_colours)[4:16],face = "bold"),
                     plot.title=element_text(size=11,hjust=.5), legend.position="right",panel.background = element_blank())+
	labs(fill = "Effect Size")+
geom_hline(yintercept = 3.5,size=1,colour="black")+
geom_hline(yintercept = 11.5,size=1,colour="black")+
geom_vline(xintercept = 2.5,size=1,colour="black")+
geom_vline(xintercept = 6.5,size=1,colour="black")+
geom_hline(yintercept = 7.5,size=1,colour="black")+
geom_vline(xintercept = 13.5,size=1,colour="black")+
geom_vline(xintercept = 10.5,size=1,colour="black")+
scale_x_discrete(position = "top") 




asp<-13/16
# define aspect ratio of PCoA plots


margins<-c(2.5,0,0,0)
# define margins of PCoA plots

# The following lines will perform Principal Coordinate Analysis to help visualise the clustering of 
# the relationships of bones to one another 

shape_dist<-vegdist(t(temporary4[4:16,1:13]),method='gower')
shapepcoa<-pcoa(as.matrix(shape_dist))

size_dist<-vegdist(t(temporary4b[4:16,1:13]),method='gower')
sizepcoa<-pcoa(size_dist)


plot_data<-as.data.frame(cbind(sizepcoa$vectors[,1],sizepcoa$vectors[,2],sizepcoa$vectors[,3]))
bone_colours<-c(rep(cbbPalette[7],2),rep(cbbPalette[4],4),rep(cbbPalette[6],4),rep(cbbPalette[8],3),rep(cbbPalette[1],3))


# This is a PCoA plot of centroid size clustering 
Size_1_2<-
ggplot()+
geom_point(data=plot_data,aes(x=V1/.25,y=V2),size=3,stroke=1/2,fill=bone_colours[1:13],shape=c(rep(21,13)))+
geom_text_repel(data=plot_data,aes(x=V1/.25,y=V2,label = rownames(plot_data)),vjust=-1.2,
colour=bone_colours[1:13], fontface='bold')+
	xlim(c(min(plot_data[,1])-.05,max(plot_data[,1])+.05))+
	ylim(c(min(plot_data[,2])-.05,max(plot_data[,2])+.05))+
  labs(x=paste("PCO-1",round(sizepcoa$values[1,3]*100,2),"% relative eigenvalue"), 
y=paste("PCO-2",round(sizepcoa$values[2,3]*100,2),"% relative eigenvalue")) +
  theme_bw(base_size=6) +theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	plot.margin=unit(margins,'cm'),
	#axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_text(size=7,colour="black"),
	axis.title.x.top=element_text(size=16,colour="black"),
	axis.title.y.left=element_text(size=16,colour="black"),
      axis.ticks=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      #panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	panel.border = element_rect(linetype = "solid", fill = NA,size=1),
aspect.ratio=asp)+
scale_x_discrete(position = "top") 

plot_data<-as.data.frame(cbind(shapepcoa$vectors[,1],shapepcoa$vectors[,2],shapepcoa$vectors[,3]))
bone_colours<-c(rep(cbbPalette[7],2),rep(cbbPalette[4],4),rep(cbbPalette[6],4),rep(cbbPalette[8],3),rep(cbbPalette[1],3))

# This is a PCoA plot of 3-D skeletal unit shape clustering 
Shape_1_2<-
ggplot()+
geom_point(data=plot_data,aes(x=V1/.25,y=V2),size=3,stroke=1/2,fill=bone_colours[1:13],shape=c(rep(21,13)))+
geom_text_repel(data=plot_data,aes(x=V1/.25,y=V2,label = rownames(plot_data)),vjust=-1.2,
colour=bone_colours[1:13], fontface='bold')+
	xlim(c(min(plot_data[,1])-.05,max(plot_data[,1])+.05))+
	ylim(c(min(plot_data[,2])-.05,max(plot_data[,2])+.05))+
  labs(x=paste("PCO-1",round(shapepcoa$values[1,3]*100,2),"% relative eigenvalue"), 
y=paste("PCO-2",round(shapepcoa$values[2,3]*100,2),"% relative eigenvalue")) +
   theme_bw(base_size=6) + theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	plot.margin=unit(margins,'cm'),
	#axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_text(size=7,colour="black"),
	axis.title.x.top=element_text(size=16,colour="black"),
	axis.title.y.left=element_text(size=16,colour="black"),
      axis.ticks=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      #panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	panel.border = element_rect(linetype = "solid", fill = NA,size=1),
aspect.ratio=asp)+
scale_x_discrete(position = "top") 


sd1<-(t(temporary4[4:16,1:13])-mean(t(temporary4[4:16,1:13]),na.rm=T))/sd(t(temporary4[4:16,1:13]),na.rm=T)
sd2<-(t(temporary4b[4:16,1:13])-mean(t(temporary4b[4:16,1:13]),na.rm=T))/sd(t(temporary4b[4:16,1:13]),na.rm=T)


object<-sd1-sd2

diff_dist<-vegdist(object,method='gower')
diffpcoa<-pcoa(as.matrix(diff_dist))

plot_data<-as.data.frame(cbind(diffpcoa$vectors[,1],diffpcoa$vectors[,2],diffpcoa$vectors[,3]))

# This plot visualises the difference in the clustering of skeletal shape and size organisation
Diff_1_2<-
ggplot()+
geom_point(data=plot_data,aes(x=V1/.25,y=V2),size=3,stroke=1/2,fill=bone_colours[1:13],shape=c(rep(21,13)))+
geom_text_repel(data=plot_data,aes(x=V1/.25,y=V2,label = rownames(plot_data)),vjust=-1.2,
colour=bone_colours[1:13], fontface='bold')+
	xlim(c(min(plot_data[,1])-.05,max(plot_data[,1])+.05))+
	ylim(c(min(plot_data[,2])-.05,max(plot_data[,2])+.05))+
  labs(x=paste("PCO-1",round(diffpcoa$values[1,3]*100,2),"% relative eigenvalue"), 
y=paste("PCO-2",round(diffpcoa$values[2,3]*100,2),"% relative eigenvalue")) +
    theme_bw(base_size=6) + theme(axis.line.x.bottom=element_line(size = 1, colour = "black"),
	plot.margin=unit(margins,'cm'),
	#axis.line.y.left=element_line(size = 1, colour = "black"),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
	axis.title.x.bottom=element_text(size=7,colour="black"),
	axis.title.x.top=element_text(size=16,colour="black"),
	axis.title.y.left=element_text(size=16,colour="black"),
      axis.ticks=element_blank(),
      legend.position="none",
      panel.background=element_blank(),
      #panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
	panel.border = element_rect(linetype = "solid", fill = NA,size=1),
aspect.ratio=asp)+
scale_x_discrete(position = "top") 


## All analysis and preparation of plots has concluded 
##
# Type the following names into the command line to view plots:
# mp1
# 'mp1' provides a plot of pairwise evolutionary integration between different bone 3-d shapes
# mp2
# 'mp2' provides a plot of pairwise evolutionary integration between different bone centroid sizes 
# mp3
# 'mp3' provides a qualitative impression of the difference between mp1 and mp2

# Size_1_2
# 'Size_1_2' provides a PCoA plot of evolutionary integration of bone centroid sizes
# Shape_1_2
# 'Shape_1_2' provides a PCoA plot of evolutionary integration of bone shapes
# Diff_1_2
# 'Diff_1_2' visualises the difference in size and shape organisation

