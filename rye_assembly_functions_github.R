
#wrappers and modifications to extant functions used in the rye genome assembly

#required libraries
library(igraph)
library(zoo)
library(stringi)

#integrates with TRITEX assembly objects
#maps coordinates on an agp object to scaffold-based coordinates
agp_coords_to_scaffs <- function(agp,sscaff_info,sd){ #sd with cols "chr" (int) and agp_coords (number)
  a <- copy(agp)[ , .(superscaffold=scaffold,scaffold_length,agp_orientation=orientation,agp_start,chr,agp_chr)] #agp
  a[ agp_chr=="chrUn" , chr := 0L ]
  a[ is.na(chr) , chr := sub("chr(\\d)R","\\1",agp_chr) %>% as.integer ]
  
  a[ is.na(agp_orientation) , agp_orientation := 1 ]
  s <- copy(sscaff_info) #all scaffs
  s[ , scaffold_length := NULL ]
  s[is.na(superscaffold) , superscaffold := scaffold]
  t <- s[ a , on=.(superscaffold) ] #transpose
  t[is.na(scaffold) , scaffold := "gap" ]
  t[is.na(scaffold)]
  t[ is.na(ss_orientation) , ss_orientation := 1 ]
  t[ superscaffold=="gap" , ss_length := scaffold_length ] #these are agp gaps (as opposed to intra-ss gaps)
  t[ superscaffold=="gap" , ss_start := 1 ]
  t[ , scaffold_agp_start := ifelse( agp_orientation == 1 , agp_start + ss_start - 1 , (agp_start + ss_length) - ss_end ) ] #add for each row a coordinate of the start on the agp
  t[ , agp_coords := scaffold_agp_start]
  o <- t[ sd , on=.( chr , agp_coords ) , roll=TRUE ] #output
  o[ , scaffold_coord := ifelse( ss_orientation == 1 , (agp_coords-scaffold_agp_start) + 1  , scaffold_length - (agp_coords-scaffold_agp_start) ) ]
  o[ , scaffold_orientation := ifelse(ss_orientation==agp_orientation , 1 , -1 ) ]
  
  o[ , .(scaffold,scaffold_coord,scaffold_orientation)]
  o
}

#integrates with TRITEX assembly objects
#maps coordinates on superscaffolds to AGP coordinates
ss_orig_coord_to_psmol <- function(agp,ss_info,orig_scaffs_and_posns,bin_length=0){
  #connect orig scaffs to positions (and to sscaffs and sscaff positions/orientations)
  osp <- copy(orig_scaffs_and_posns)[,idx := 1:.N] #necessary for some bullshit where asterisks end up in there.
  s <- ss_info[scaffold != "gap" , .(orig_scaffold,orig_pos = orig_start,orig_scaff_start=orig_start,scaffold,superscaffold,ss_orientation,ss_start,ss_end) ][osp,on=.(orig_scaffold,orig_pos),roll=TRUE]
  #s[orig_scaffold=="scaffold38"] %>% ggplot(aes(x=orig_pos,y=1)) + geom_point()
  
  #get agp details (start pos, orient) of sscaffs from agp
  s[is.na(superscaffold) , superscaffold := scaffold ]
  ss <- agp[scaffold != "gap" , .(agp_chr,superscaffold = scaffold,agp_orientation = orientation , agp_start , agp_end )][ s , on = .(superscaffold) ]
  ss[is.na(agp_orientation),agp_orientation := 1]
  ss[is.na(ss_orientation),ss_orientation := 1]
  
  
  #make scaff_coords col
  ss[ , scaff_coords := orig_pos - orig_scaff_start + 1 ]
  #ss[orig_scaffold=="scaffold38"] %>% ggplot(aes(x=scaff_coords,y=scaffold)) + geom_point() + xlim(c(1,5e6))
  #make sscaff_coords col
  ss[ , sscaff_coords := ifelse(ss_orientation>0 , ss_start + scaff_coords - 1 , ss_end - scaff_coords + 1 - bin_length + 1 ) ]
  #ss[orig_scaffold=="scaffold38"] %>% ggplot(aes(x=sscaff_coords,y=superscaffold,colour=orig_pos)) + geom_point()
  #make agp_coords col
  ss[ , agp_coords := ifelse(agp_orientation>0 , agp_start + sscaff_coords - 1 , agp_end - sscaff_coords + 1 - bin_length + 1 ) ]
  #ss[orig_scaffold=="scaffold38"] %>% ggplot(aes(x=agp_coords,y=agp_chr,colour=superscaffold,alpha=orig_pos)) + geom_jitter(height = .2)
  ss <- ss[order(idx)]
  ss
}


#integrates with TRITEX assembly objects
#maps coordinates of a genetic map ("popseq") object to those of the AGP
map_to_agp_coords <- function(assembly,hicagp,transpose=FALSE,map=NULL){
  if(is.null(map)){
    map <- copy(assembly$popseq)
  }
  if(transpose==TRUE){
    copy(map) -> ps
    if(is.null(ps$orig_scaffold)) { ps[, orig_scaffold := scaffold] }
    if(is.null(ps$orig_scaffold_length)) { ps[, orig_scaffold_length := scaffold_length] }
    if(is.null(ps$orig_pos)) { ps[, orig_pos := pos] }
    ps[,scaffold := NULL]
    ps[,scaffold_length := NULL]
    ps[,pos := NULL]
    assembly$info[, .(scaffold, scaffold_length=length, orig_scaffold, orig_start, orig_pos=orig_start) ][ps, on=c("orig_scaffold", "orig_pos"), roll=T]->ps
    ps[, pos := orig_pos - orig_start + 1]
    ps[, orig_start := NULL]
    map <- ps
  }
  if(is.null(map$scaffold) & !is.null(map$css_contig)){
    #get pos and scaffold from cssaln
    p <- assembly$cssaln[,.(marker_name = css_contig, scaffold, marker_pos_on_scaffold = pos , marker_css_chr = sorted_chr )][map[, .( marker_name = css_contig , marker_popseq_lg_chr = popseq_chr , marker_popseq_cM = popseq_cM ) ] , on="marker_name" ]
  } else {
    p <- map[,.(scaffold, marker_pos_on_scaffold = pos, marker_popseq_lg_chr = popseq_chr , marker_popseq_cM = popseq_cM)]
  }
  hicagp[,.(scaffold,agp_chr,agp_popseq_cM = popseq_cM , scaffold_length , agp_orientation = orientation ,  agp_start , agp_end )] -> a
  agp_popseq <- a[p,on="scaffold"]; rm(p,a)
  agp_popseq[, marker_agp_pos := ifelse( is.na(agp_orientation) , agp_start + marker_pos_on_scaffold - 1 , ifelse ( agp_orientation > 0 , agp_start + marker_pos_on_scaffold - 1 , agp_end - marker_pos_on_scaffold + 1 )  )  ][]
}
#plot_optigs_agp_detail(data = rye_v4p3$optig_mapping , hicmap=rye_v4p3$hic_maps$rye_v4_hic_map_v1 , selchr = "chr1R" , range = c(160e6,165e6) , min_nonagp_scafflen=1e4 , obreaks = NULL , tbreaks = NULL , agp = NULL )

#plot_genetic_map_vs_agp( assembly = r , hicagp = r$hic_maps$rye_v4p5_hic_map_v2$agp , transpose = FALSE );
plot_genetic_map_vs_agp <- function(assembly,hicagp,transpose=FALSE,map=NULL){
  m <- map_to_agp_coords(assembly,hicagp,transpose,map)
  #write.csv( m[,.(agp_chr=agp_chr,agp_position=marker_agp_pos,cM=marker_popseq_cM,linkage_group=marker_popseq_lg_chr)]   ,file="/filer/transfer/bauer17_gmap_to_Sc_v4p5_agp_coords.csv")
  
  ggplot( data=m , aes(x = marker_agp_pos , y = marker_popseq_cM , colour = as.factor(marker_popseq_lg_chr))) + 
    geom_point(data=m , size=.2)  + 
    geom_vline(data=hicagp , mapping=aes(xintercept = agp_start),alpha = .2 , colour="black") + 
    facet_grid(agp_chr~.) +
    geom_text(data=hicagp , aes(x=agp_start,y=100,label=ifelse(scaffold_length > 500 , sub(".*_","",scaffold) , "")),angle=90,size=3,colour="black")
}
#write.csv( m[,.(agp_chr=agp_chr,agp_position=marker_agp_pos,cM=marker_popseq_cM,linkage_group=marker_popseq_lg_chr)]   ,file="bauer17_gmap_to_Sc_v4p5_agp_coords.csv")




#integrates with TRITEX assembly objects
#reads hic objects as per the TRITEX 'read_fragdata' but accounts for superscaffold level assembly
read_fragdata_ss<-function(info, file, ss_info){
  fragbed<-fread(file, head=F, col.names=c("orig_scaffold", "start", "end"))
  fragbed[, length := end - start]
  fragbed[, start := start + 1]
  
  #first transpose to broken
  fragbed_tobroken <-  ss_info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]
  fragbed_tobroken[, start := start - orig_start + 1]
  fragbed_tobroken[, end := end - orig_start + 1]
  fragbed_tobroken[, orig_start := NULL]
  fragbed_tobroken[, orig_scaffold := NULL]
  
  #then to ss
  fragbed_tosuperscaffs <- ss_info[ fragbed_tobroken , on="scaffold" ]
  fragbed_tosuperscaffs_ss <- fragbed_tosuperscaffs[!is.na(superscaffold)]
  fragbed_tosuperscaffs_ss[,scaffold := superscaffold]
  fragbed_tosuperscaffs_ss[,prev_start := start ]
  fragbed_tosuperscaffs_ss[,prev_end := end ]
  fragbed_tosuperscaffs_ss[,start := ifelse(ss_orientation>0 , (ss_start-1) + prev_start , (ss_end+1) - prev_end ) ]
  fragbed_tosuperscaffs_ss[,end := ifelse(ss_orientation>0 , (ss_start-1) + prev_end , (ss_end+1) - prev_start ) ]
  
  fragbed_tosuperscaffs_no_ss <- fragbed_tosuperscaffs[is.na(superscaffold)]
  
  fragbed <- bind_rows( fragbed_tosuperscaffs_no_ss[,.(scaffold,start,end)] , fragbed_tosuperscaffs_ss[,.(scaffold,start,end)] )
  
  
  info <- fragbed[, .(nfrag = .N), keyby=scaffold][info, on="scaffold"][is.na(nfrag), nfrag := 0]
  
  list(bed=fragbed[], info=info[])
}
#read_fragdata_ss(info=rye_v4$info,file=hic_fragfile,ss_info=rye_v4$ss_info)



#integrates with TRITEX assembly objects
#breaks scaffolds per the TRITEX 'break_scaffolds' but accounts for superscaffold level assembly
break_scaffolds_ss<-function(breaks, assembly, prefix, slop, cores=1, species="wheat"){
  regex1 <- "(^.*[^-0-9])(([0-9]+)(-[0-9]+)?$)"
  regex2 <- "(^.*[^-0-9])([0-9]+(-[0-9]+)?$)"
  
  info <- assembly$info
  cov <- assembly$cov
  fpairs <- assembly$fpairs
  cssaln <- assembly$cssaln
  molecules <- assembly$molecules
  
  br <- copy(breaks)[,.(scaffold,br)]
  fai <- info[, .(scaffold, orig_scaffold, old_scaffold=scaffold, orig_start, orig_end, length)]
  
  cat("Split scaffolds\n")
  j <- 0
  while(nrow(br) > 0){
    j <- j + 1
    o <- nrow(fai)
    fai[br, on="scaffold"] -> br
    br[, orig_br := orig_start + br - 1]
    br[order(scaffold, br)] -> br
    br[duplicated(scaffold)] -> nbr
    br[!duplicated(scaffold)] -> br
    
    max(as.integer(sub(regex1, "\\3", fai$scaffold)))->maxidx
    br[, idx := 3*1:.N-2]
    br[, scaffold1 := paste0(prefix, maxidx+idx)]
    br[, start1 := 1]
    br[, end1 := pmax(0, br - slop - 1)]
    br[, scaffold2 := paste0(prefix, maxidx+idx+1)]
    br[, start2 := pmax(1, br - slop)]
    br[, end2 := pmin(br + slop - 1, length)]
    br[, scaffold3 := paste0(prefix, maxidx+idx+2)]
    br[, start3 := pmin(length + 1, br + slop)]
    br[, end3 := length]
    br[, length1 := 1 + end1 - start1]
    br[, length2 := 1 + end2 - start2]
    br[, length3 := 1 + end3 - start3]
    rbind(
      br[, .(scaffold=scaffold1, length=length1, orig_scaffold, orig_start=orig_start+start1-1, orig_end=orig_start+end1-1, old_scaffold)],
      br[, .(scaffold=scaffold2, length=length2, orig_scaffold, orig_start=orig_start+start2-1, orig_end=orig_start+end2-1, old_scaffold)],
      br[, .(scaffold=scaffold3, length=length3, orig_scaffold, orig_start=orig_start+start3-1, orig_end=orig_start+end3-1, old_scaffold)],
      fai[!scaffold %in% br$scaffold, .(orig_scaffold, length,  orig_start, orig_end, old_scaffold,
                                        scaffold=paste0(prefix, sub(regex2, "\\2", scaffold)))]
      
    ) -> fai
    
    cat(paste0("Iteration ", j, " finished. "))
    fai[length > 0]->fai
    cat(paste0("The number of scaffolds increased from ", o, " to ", nrow(fai), ".\n"))
    fai[, .(scaffold, orig_scaffold, orig_start, orig_br=orig_start)][nbr[, .(orig_scaffold, orig_br)], on=c("orig_scaffold", "orig_br"), roll=T]->nbr
    nbr[, br := orig_br - orig_start + 1]
    nbr[, .(scaffold, br)] -> br
    cat("Br remaining:\n")
    print(br)
  }
  
  fai[, split := F]
  fai[old_scaffold %in% breaks$scaffold, split := T]
  fai[, old_scaffold := NULL]
  
  assembly_new<-list(info=fai)
  
  cat("Transpose optig mapping\n")
  
  #browser() ############################################################################################################################################################
  ########################################################################################################################################################################################
  
  
  if(!is.null(assembly$optig_mapping)){
    oplinks <- copy(assembly$optig_mapping) #update scaffold , scaffold pos
    
    if (is.null(oplinks$orig_scaffold)) {oplinks[,orig_scaffold := scaffold ]}
    if (is.null(oplinks$orig_scaffold_pos)) {oplinks[,orig_scaffold_pos := scaffold_pos ]}
    if (!is.null(oplinks$orig_scaffold_start)) {oplinks[,orig_scaffold_start := NULL ]}
    
    oplinks[,scaffold := NULL]
    oplinks[,scaffold_pos := NULL]
    
    #dev
    # oplinks[,old_cmap_pos := cmap_pos]
    # oplinks[,old_cmap_length := cmap_length]
    
    #work on all
    #create a handy converter and attach new scaffold, start, length
    s1 <- fai[, .( orig_scaffold , orig_scaffold_pos = orig_start ,              scaffold, orig_scaffold_start = orig_start , ref_new_cmap_length = length )]
    t <- s1[oplinks, on=c("orig_scaffold","orig_scaffold_pos") , roll=T]
    
    #dev. in the second it both s1 and oplinks have orig_scaffold_start cols. The important one is from fai. So always delete the old and import from oplinks and import the new from fai.
    
    
    #then pick only refs and update cmap pos (=scaffold_pos) and cmap_length (=)
    t[ ref_query=="ref", cmap_pos := orig_scaffold_pos - orig_scaffold_start + 1 ]
    t[ ref_query=="ref", cmap_length := ref_new_cmap_length ]
    t[,ref_new_cmap_length := NULL ]
    
    #update xmap_match_id
    #list broken
    broken_scaffolds <- fai[split == TRUE ,orig_scaffold] %>% unique
    #new numbers
    global_max_xmap_id <- max(t$xmap_match_id) %>% as.integer
    assign_xmap_newid <- function() {(global_max_xmap_id <<- global_max_xmap_id + 1L)}
    t[orig_scaffold %in% broken_scaffolds, xmap_match_id := assign_xmap_newid()  ,.(scaffold)]
    
    assembly_new$optig_mapping <- t
  }
  
  
  
  
  cat("Transpose cssaln\n")
  
  copy(cssaln) -> z
  z[, scaffold_length := NULL]
  z[, scaffold := NULL]
  fai[, .(scaffold, scaffold_length=length, orig_scaffold, orig_start, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_pos"), roll=T]->z
  z[, pos := orig_pos - orig_start + 1]
  z[, orig_start := NULL]
  assembly_new$cssaln <- z
  
  
  #detect scaff col in popseq. if there, update table to new scaffs
  if( !is.null(assembly$popseq$scaffold) ){
    cat("Transposing popseq\n")
    copy(assembly$popseq) -> ps
    if(is.null(ps$orig_scaffold)) { ps[, orig_scaffold := scaffold] }
    if(is.null(ps$orig_scaffold_length)) { ps[, orig_scaffold_length := scaffold_length] }
    if(is.null(ps$orig_pos)) { ps[, orig_pos := pos] }
    ps[,scaffold := NULL]
    ps[,scaffold_length := NULL]
    ps[,pos := NULL]
    fai[, .(scaffold, scaffold_length=length, orig_scaffold, orig_start, orig_pos=orig_start) ][ps, on=c("orig_scaffold", "orig_pos"), roll=T]->ps
    ps[, pos := orig_pos - orig_start + 1]
    ps[, orig_start := NULL]
    assembly_new$popseq <- ps
    
    #check it worked
    
  } else {
    stop("Popseq info tranpose error: No scaffold column detected")
  }
  
  print("New popseq  :")
  print(assembly_new$popseq)
  
  if("fpairs" %in% names(assembly) && nrow(fpairs) > 0){
    cat("Transpose fpairs\n")
    copy(assembly$fpairs)[, .(orig_scaffold1, orig_scaffold2, orig_pos1, orig_pos2)]->z
    fai[, .(scaffold1=scaffold, orig_scaffold1=orig_scaffold, orig_start1=orig_start, orig_pos1=orig_start)][z, on=c("orig_scaffold1", "orig_pos1"), roll=T]->z
    fai[, .(scaffold2=scaffold, orig_scaffold2=orig_scaffold, orig_start2=orig_start, orig_pos2=orig_start)][z, on=c("orig_scaffold2", "orig_pos2"), roll=T]->z
    z[, pos1 := orig_pos1 - orig_start1 + 1]
    z[, pos2 := orig_pos2 - orig_start2 + 1]
    z[, orig_start1 := NULL]
    z[, orig_start2 := NULL]
    assembly_new$fpairs <- z
  } else {
    assembly_new$fpairs <- data.table()
  }
  
  # print("New fpairs  :")
  # print(assembly_new$fpairs)
  
  if("molecules" %in% names(assembly) && nrow(molecules) > 0){
    cat("Transpose molecules\n")
    copy(molecules) -> z
    z[, scaffold := NULL]
    fai[, .(scaffold, orig_scaffold, orig_start, s_length=orig_end - orig_start + 1, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_start"), roll=T]->z
    z[, start := orig_start - orig_pos + 1]
    z[, end := orig_end - orig_pos + 1]
    z[end <= s_length]->z
    z[, orig_pos := NULL]
    z[, s_length  := NULL]
    assembly_new$molecules <- z
  } else {
    assembly_new$molecules <- data.table()
  }
  
  assembly_new$breaks <- breaks
  
  cat("Anchor scaffolds\n")
  
  # print("With fpairs as :")
  # print(assembly_new$fpairs)
  
  
  anchor_scaffolds(assembly_new, popseq=assembly_new$popseq, species=species )->assembly_new
  
  #print("Check [scaffold == s_112246]:")
  #print(assembly_new$info[scaffold == "s_112246"])
  
  if("mr_10x" %in% names(assembly_new$info)){
    assembly_new$info[, mr_10x := NULL]
  }
  
  if("mr" %in% names(assembly_new$info)){
    assembly_new$info[, mr := NULL]
    assembly_new$info[, mri := NULL]
  }
  
  if("cov" %in% names(assembly) & nrow(fpairs) > 0){
    cat("Hi-C coverage\n")
    add_hic_cov(assembly_new, scaffolds=fai[split == T]$scaffold, binsize=assembly$binsize, minNbin=assembly$minNbin, innerDist=assembly$innerDist, cores=cores)->cov
    assembly$cov[!scaffold %in% breaks$scaffold]->x
    x[, scaffold:=paste0(prefix, sub(regex2, "\\2", scaffold))]
    if(nrow(cov$cov) > 0){
      rbind(x, cov$cov)->assembly_new$cov
    } else {
      x -> assembly_new$cov
    }
    info[!scaffold %in% breaks$scaffold]->x
    x[, scaffold:=paste0(prefix, sub(regex2, "\\2", scaffold))]
    x[, split := F]
    rbind(x[, names(cov$info), with=F], cov$info)->assembly_new$info
  } else {
    assembly_new$cov <- data.table()
  }
  
  
  
  if("molecule_cov" %in% names(assembly) & nrow(molecules) > 0){
    cat("10X molecule coverage\n")
    add_molecule_cov(assembly_new, scaffolds=fai[split == T]$scaffold, binsize=assembly$mol_binsize, cores=cores)->cov
    info[!breaks$scaffold, on="scaffold"]->x
    x[, scaffold := paste0(prefix, sub(regex2, "\\2", scaffold))]
    x[, split := F]
    rbind(x[, names(cov$info), with=F], cov$info)->assembly_new$info
    assembly_new$mol_binsize <- assembly$mol_binsize
    
    assembly$molecule_cov[!breaks$scaffold, on="scaffold"]->x
    x[, scaffold := paste0(prefix, sub(regex2, "\\2", scaffold))]
    if(nrow(cov$molecule_cov) > 0){
      rbind(x, cov$molecule_cov)->assembly_new$molecule_cov
    } else {
      x -> assembly_new$molecule_cov
    }
  } else {
    assembly_new$molecule_cov <- data.table()
  }
  
  assembly_new$binsize <- assembly$binsize
  assembly_new$innerDist <- assembly$innerDist
  assembly_new$minNbin <- assembly$minNbin
  
  
  return(assembly_new)
}
#break_scaffolds( breaks =  rye_v4$optig_scaffbreaks, assembly = rye_v4, prefix="scaffolds_optigbreak_", slop=1e4, cores=32, species = "rye")



#integrates with TRITEX assembly objects
#plots Hi-C asymetry ratios along the pseudomolecules
plot_asymmetry_agp <- function(agp,asymmetry){
  ggplot(asymmetry$ratio,aes(x=bin+5e5,y=r)) + 
    geom_point(size=.1) +
    geom_vline(data=agp,aes(xintercept=agp_start),size=.1,alpha=.2) + 
    geom_text(data=agp , aes(x=agp_start,y=-100,label=scaffold),angle=90,size=1,colour="black" ) +
    facet_grid(chr~.)
}
#plot_asymmetry_agp(rye_v2p2$hic_maps$rye_v2p2_hic_map_v1$agp,rye_v2p2$hic_maps$rye_v2p2_hic_map_v1$hic_asymmetry)


#integrates with TRITEX assembly objects
#plots genetic map cM assignments and marker positions along the pseudomolecules
plot_genetic_map_vs_agp <- function(assembly,hicagp,transpose=FALSE,map=NULL){
  m <- map_to_agp_coords(assembly,hicagp,transpose,map)
  ggplot( data=m , aes(x = marker_agp_pos , y = marker_popseq_cM , colour = as.factor(marker_popseq_lg_chr))) + 
    geom_point(data=m , size=.05)  + 
    geom_vline(data=hicagp , mapping=aes(xintercept = agp_start),alpha = .2 , size=.2 , colour="black") + 
    facet_grid(agp_chr~.) +
    geom_text(data=hicagp[scaffold != "gap"] , aes(x=agp_start,y=100,label= scaffold),angle=90,size=.8,colour="black")
}
#plot_genetic_map_vs_agp( assembly = rye_v2p2 , hicagp = rye_v2p2$hic_maps$rye_v2p2_hic_map_v1$agp  )



#integrates with TRITEX assembly objects
#reads a optical map .cmap file for integration of optical map data
read_cmap <- function(file,keyfile=NULL){
  # cmap <- fread(paste0("type ", file, " | find /V \"#\"") ,
  cmap <- fread(paste0("cat ", file, " | grep -v '^#'") ,
                col.names = c(
                  "cmap_contig_id",
                  "cmap_contig_length",
                  "cmap_contig_nsites",
                  "cmap_label_id",
                  "cmap_label_channel",
                  "cmap_label_position",
                  "cmap_label_interval_sd",
                  "cmap_label_interval_coverage",
                  "cmap_label_occurence",
                  "cmap_label_ChimQuality",
                  "cmap_label_SegDupL",
                  "cmap_label_SegDupR",
                  "cmap_label_FragileL",
                  "cmap_label_FragileR",
                  "cmap_label_OutlierFrac",
                  "cmap_label_ChimNorm"
                )
  )
  if(!is.null(keyfile)){
    key <- fread(paste0("cat ", keyfile, " | grep -v '^#'") , col.names = c("cmap_contig_id","cmap_insilico_ref","keyfile_length"))
    #key <- fread(paste0("type ", keyfile, " | find /V \"#\"") , col.names = c("cmap_contig_id","cmap_insilico_ref","keyfile_length"))
    return(key[cmap,,on="cmap_contig_id"])
  } else {
    return(cmap)
  }
}

#integrates with TRITEX assembly objects
#reads a optical map cut file
read_opcut <- function(file,keyfile){
  cutfile_colnames <- c("xmap_match_id" , "ref_query" , "cmap_id" , "break_left" , "break_right" , "orientation" , "cut_left" , "cut_right" , "discard")
  key <- fread(paste0("cat ", keyfile, " | grep -v '^#'") , col.names = c("xmap_ref_id","cmap_insilico_ref","keyfile_length"))
  r <- fread(paste0('cat ',file,' | grep -v "^#"'),header=FALSE)
  r <- bind_rows(r[,1:9] %>% set_colnames(cutfile_colnames),r[,c(1,10:17)]  %>% set_colnames(cutfile_colnames))
  r[,cut_left := cut_left=="cut"]
  r[,cut_right := cut_right=="cut"]
  r[,discard := discard=="cut"]
  r[break_left==-1,break_left:=NA]
  r[break_right==-1,break_right:=NA]
  #r[discard==TRUE] #it's none anyway but this would tell us whether to throw away the bit to the left/right of the cut ... could be used to influence colur/shape in future ... ?
  r <- bind_rows(r[cut_left==TRUE,][,cut := break_left],r[cut_right==TRUE,][,cut := break_right])
  r <- bind_rows( key[,.(scaffold=cmap_insilico_ref,cmap_id = xmap_ref_id )][r[ref_query=="ref"],on="cmap_id"] , r[ref_query=="qry"][, scaffold := as.character(cmap_id) ] )
  r[,ref_query := sub("qry","query",ref_query)]
  r
}
#r <- read_opcut(file=optical_cutfile,keyfile=cmap_Sc_v1_keyfile)



#integrates with TRITEX assembly objects
#an adaptation of 'break_scaffolds' to apply breaks recorded in an optical map cut file
break_optigs <- function(optig_links,breaks){
  
  output_cols <- colnames(optig_links)
  
  optigs <- copy(optig_links)
  breaks <- copy(breaks)
  #undo all and re-establish later: this only works if we do one round of breaking (and why would we ever do more?)
  optigs[,orig_cmap_id := cmap_id][,cmap_id := NULL]
  optigs[,orig_cmap_pos := cmap_pos][,cmap_pos := NULL]
  optigs[,orig_cmap_length := cmap_length][,cmap_length := NULL]
  optigs[,orig_xmap_match_id := xmap_match_id][,xmap_match_id := NULL]
  optigs[,orig_query := query][,query := NULL]
  
  #lookup table for optig lengths
  optig_info <- optigs[ref_query=="query", .SD[,.(orig_cmap_length)][1] ,.(orig_cmap_id)]
  
  #adjust name to allow matching
  setnames(breaks,"scaffold","orig_cmap_id")
  
  #separate out non-changers from changers
  optigs[ref_query=="query" & orig_cmap_id %in% breaks$orig_cmap_id]$link_uniq_tag -> changers_linktags
  optigs[link_uniq_tag %in% changers_linktags] -> oplinks_tochange
  optigs[!link_uniq_tag %in% changers_linktags] -> oplinks_nochange
  
  #nonchangers don't need to record original info
  c <- colnames(oplinks_nochange)
  setnames(oplinks_nochange,c,sub("orig_","",c))
  
  #set up output object
  new_optig_cmaps <- data.table( 
    cmap_id = character(), #to add
    cmap_length = numeric(), #to add
    orig_cmap_id = character(), #to join on
    orig_cmap_pos = integer(), # to join on (roll)
    orig_cmap_start = integer() #to add (saves start position in orig coords during roll)
  )
  
  #FUNCTIONS #define here to assure parent.frame and enclosing env are appropriate
  
  #new naming function. will create n new unique cmap IDs when poked
  max_cmap_id <- as.integer(optigs$orig_cmap_id) %>% max
  make_new_cmap_ids <- function(n){
    sapply(1:n, function(x) {(max_cmap_id <<- max_cmap_id + 1)} ) %>% as.character
  }
  
  max_xmap_match_id <- as.integer(optigs$orig_xmap_match_id) %>% max
  make_new_xmap_match_ids <- function(n){
    sapply(1:n, function(x) {(max_xmap_match_id <<- max_xmap_match_id + 1)} ) %>% as.character
  }
  
  #browser()
  
  #breaking function, taking a scaffold and all its breakpoints
  break_cmap <- function(id,br){
    
    cat("Breaking",id,"at positions:",br,"\n")
    
    
    q <- oplinks_tochange[ref_query=="query" & orig_cmap_id==id]
    r <- oplinks_tochange[ref_query=="ref" & orig_query==id]
    
    
    
    q_len <- q$orig_cmap_length[1]
    r_len <- r$orig_cmap_length[1]
    
    q_joiner <- data.table( cmap_id = make_new_cmap_ids(length(br) + 1) , xmap_match_id = make_new_xmap_match_ids(length(br) + 1) ,  new_cmap_start = c(0,br) , orig_cmap_pos = c(0,br) , cmap_length = c(br,q_len) - c(0,br) )
    
    q_joined <- q_joiner[q,on="orig_cmap_pos",roll=T]
    q_joined[, query := cmap_id ]
    q_joined[, cmap_pos := orig_cmap_pos - new_cmap_start ]
    q_joined[, xmap_match_id := as.integer(xmap_match_id)]
    
    if (q[,.N,.(orig_cmap_id,orig_cmap_length)][,sum(orig_cmap_length)] != q_joined[,..output_cols][,.N,.(cmap_id,cmap_length)][,sum(cmap_length)]) {
      cat("Unmatched end removed from",id,"resulting in length reduction from",q[,.N,.(orig_cmap_id,orig_cmap_length)][,sum(orig_cmap_length)],"to",q_joined[,..output_cols][,.N,.(cmap_id,cmap_length)][,sum(cmap_length)],"... \n")
    }
    
    r_joiner <- q_joined[,.(query,xmap_match_id),by="link_uniq_tag"]
    
    r_joined <- r_joiner[r,on="link_uniq_tag"]
    r_joined[, cmap_id := orig_cmap_id]
    r_joined[, cmap_pos := orig_cmap_pos]
    r_joined[, cmap_length := orig_cmap_length]
    
    bind_rows( r_joined[,..output_cols] , q_joined[,..output_cols] )
  }
  
  
  breaks[,break_cmap(id=id,br=br),.(id=orig_cmap_id)][,id := NULL] %>% bind_rows ( oplinks_nochange[,..output_cols] )
}
#break_optigs(xmap_links_melt,optig_opbreaks)


#integrates with TRITEX assembly objects
#reads an optical map alignment file (xmap)
read_xmap <- function(file,keyfile){
  #xmap <- fread(paste0("type ", file) ,
  xmap <- fread(paste0("cat ", file, " | grep -v '^#'") ,
                col.names = c(
                  "xmap_match_id",
                  "xmap_query_id",
                  "xmap_ref_id",
                  "xmap_qstart",
                  "xmap_qend",
                  "xmap_rstart",
                  "xmap_rend",
                  "xmap_orientation",
                  "xmap_confidence",
                  "xmap_cigar",
                  "xmap_qlength",
                  "xmap_rlength",
                  "xmap_label_channel",
                  "xmap_alignment"
                ),
                skip = 10
  )
  key <- fread(paste0("cat ", keyfile, " | grep -v '^#'") , col.names = c("xmap_ref_id","cmap_insilico_ref","keyfile_length"))
  #key <- fread(paste0("type ", keyfile, " | find /V \"#\"") , col.names = c("xmap_ref_id","cmap_insilico_ref","keyfile_length"))
  key[xmap,,on="xmap_ref_id"][,xmap_query_id := as.character(xmap_query_id)][,xmap_ref_id := as.character(xmap_ref_id)]
}

#integrates with TRITEX assembly objects
#plots optig alignments in against the AGP in a requested range
plot_optigs_agp <- function(data , hicmap , selchr = 1 , range = c(0,30e6) , obreaks = NULL , tbreaks = NULL , agp = NULL ){
  #filter here so only matches to plotted refs are used to position query
  par(mar = c(0,0,0,0)) #for regular plot operations, here used for text only
  p <- copy(data)
  #put in a melter for the data here, so we can work with it in a more sensible form otherwise
  setkey(hicmap$agp,scaffold)
  setkey(p,scaffold)
  p <- merge(hicmap$agp[scaffold != "gap",.(scaffold , agp_start , agp_end , agp_chr , agp_orientation = orientation )],p,all=T)
  p <- p[grepl(selchr,agp_chr) & (agp_start %between% range | agp_end %between% range | ( agp_start < range[1] & agp_end > range[2] ) ) ]
  
  if(nrow(p) < 1) {return(  ggplot(data.table(x=1,y=1,text=paste0("No AGP scaffolds overlapping request range: ",range[1]," -- ",range[2],"\n")),aes(x=x,y=y))  + geom_text(aes(label = text)) +
                              theme_bw() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                 axis.title.x=element_blank(),
                                                 axis.title.y=element_blank(),legend.position="none",
                                                 panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                 panel.grid.minor=element_blank(),plot.background=element_blank()
                              )
  )
  }
  
  #work out hpos for labels on ref
  pr <- p[ref_query=="ref"][, hpos := as.numeric ( agp_start + cmap_pos ) ][]
  pr[ agp_orientation==-1 , hpos := agp_end - cmap_pos ]
  
  pq <- pr[ , .(mean_refpos = mean(hpos)) , query ][p[ref_query=="query",.(mean_querypos = mean(cmap_pos)),query] , on = "query"][p[ref_query=="query"],on="query"][,hpos := mean_refpos - mean_querypos + cmap_pos ][,mean_refpos := NULL][,mean_querypos := NULL][]
  #browser()
  #get display row onto pq
  
  if(nrow(pq) < 1){return(  ggplot(data.table(x=1,y=1,text=paste0("No optig data for scaffolds overlapping request range: ",range[1]," -- ",range[2],"\n")),aes(x=x,y=y))  + geom_text(aes(label = text)) +
                              theme_bw() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                 axis.title.x=element_blank(),
                                                 axis.title.y=element_blank(),legend.position="none",
                                                 panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                 panel.grid.minor=element_blank(),plot.background=element_blank()
                              )
  )
  }
  
  p_optigs <- pq[, .( hpos = c( min(hpos) - min(cmap_pos)  , min(hpos) - min(cmap_pos) + (cmap_length[1]) ) , cmap_length = rep(cmap_length[1],2)  )  , .(grouper = query) ]
  optig_info <- p_optigs[,.SD[1,.(cmap_length,hpos)],.(cmap_id = grouper)]
  
  #work out vpos
  setkey(optig_info,hpos)
  disp_rowinfo_aux <- rep(min(optig_info$hpos)-1,100)
  optig_info[, vpos := ldply(1:nrow(optig_info),function(i){
    bit <- copy(optig_info[i,])
    first_free_disprow <- which(disp_rowinfo_aux < bit$hpos)[1] #first row where the end of the last thing is far enough back to let this one in.
    #putting it in
    if(first_free_disprow == length(disp_rowinfo_aux)){
      disp_rowinfo_aux <<- c(disp_rowinfo_aux,rep(0,100))
    }
    disp_rowinfo_aux[first_free_disprow] <<- bit$hpos + bit$cmap_length + 10e6
    first_free_disprow
  })]; disp_rowinfo_aux <- NULL
  
  #get scaff info
  pr[,vpos := 0]
  pq <- optig_info[,.(cmap_id , vpos)][pq,on="cmap_id"]
  
  #final curation of plot data
  if (!is.null(obreaks)){
    p_obreaks <- hicmap$agp[,.(scaffold,chromosome = agp_chr,agp_start,agp_end)][obreaks,on="scaffold"][!is.na(agp_start) & grepl(selchr,chromosome) & (agp_start %between% range | agp_end %between% range | ( agp_start < range[1] & agp_end > range[2] ) ) ][,br := br+agp_start][]
  }
  if (!is.null(obreaks)){
    p_tbreaks <- hicmap$agp[,.(scaffold,chromosome = agp_chr,agp_start,agp_end)][tbreaks,on="scaffold"][!is.na(agp_start) & grepl(selchr,chromosome) & (agp_start %between% range | agp_end %between% range | ( agp_start < range[1] & agp_end > range[2] ) ) ][,br := br+agp_start][]
  }
  p_optigs <- optig_info[,.(grouper = cmap_id,vpos)][p_optigs,on="grouper"][,cmap_length := NULL]
  p_links <- rbind(pq,pr)[ , colourer := as.factor(xmap_match_id) ][ , grouper := link_uniq_tag ][]
  p_scaffolds <- pr[, .( hpos = c(agp_start[1],agp_end[1]) , vpos = 0 ) , .(grouper = scaffold)]
  p_requestrange <- data.table(vpos=c(0,0),hpos=range,grouper=c(1,1))
  
  ggplot(data=p_links,aes(x = hpos , y = vpos , group = grouper ) ) +
    geom_line(aes(colour = colourer),size=.05) +
    geom_line(data=p_optigs,size=1) +
    geom_line(data=p_scaffolds,size=2) +
    geom_line(data=p_requestrange,size=2,colour="red") +
    geom_point(data=p_scaffolds,size=5,shape=124) +
    geom_text(data=p_scaffolds[,.SD[which(hpos==min(hpos))],by="grouper"],aes(label=grouper),hjust=-.3,vjust=.8,size=3,angle=90) +
    geom_text(data=p_optigs[,.SD[which(hpos==min(hpos))],by="grouper"],aes(label=grouper),hjust=-.3,vjust=.8,size=3,angle=90) +
    theme(legend.position="none") +
    ggtitle(paste0("Chromosome: ",selchr,"; Range: ",(min(p_scaffolds$hpos)/1e6) %>% round," -- ",(max(p_scaffolds$hpos)/1e6) %>% round," Mb")) +
    switch(is.null(obreaks)+1,geom_vline(data=p_obreaks,aes(xintercept=br),colour="orange",linetype="dashed",size=1),NULL) +
    switch(is.null(tbreaks)+1,geom_vline(data=p_tbreaks,aes(xintercept=br),colour="turquoise4"),NULL)
}
#plot_optigs_agp( data = rye_v1$optig_mapping , hicmap = rye_v1_hic_map_v1 , selchr = chr , range = rge , obreaks = rye_v1$optig_opbreaks , tbreaks = rye_v1$tenix_breaks )
#debugonce(plot_optigs_agp)
#plot_optigs_agp( data = rye_v4p3$optig_mapping , hicmap = rye_v4p3_hic_map_v1 , selchr = "chr1R" , range = c(20383839,41383839) )


#integrates with TRITEX assembly objects
#plots optig alignments in against the whole AGP in a multipage PDF
multipage_optig_agp <- function( data , hicmap , windowsize = 30e6 , stepsize = 28e6 , fname = "agp_optigs.pdf" , ... ){
  pdf(fname,width=30,height=7)
  #ask <- character()
  
  setkey(hicmap$agp,scaffold)
  setkey(data,scaffold)
  assembly <- merge(hicmap$agp[scaffold != "gap",.(scaffold , agp_start , agp_end , agp_chr)],data,all=T)
  
  totagplength <- assembly[,.SD[1],scaffold][,.(scaffold,agp_start,agp_end)][,sum(agp_end-agp_start)+(100*.N)]
  agp_sofar <- 0
  for(chr in unique(assembly$agp_chr %>% sort)){
    totrange <- c( min(assembly$agp_start) , max(assembly$agp_end) )
    starts <- seq(totrange[1]-1,totrange[2]+1,stepsize)
    for (s in starts){
      cat("Calculating Chromosome ", chr ," from ",s/1e6 %>% round," to ",(s+windowsize)/1e6 %>% round,"Mb (approx ", round((agp_sofar/totagplength)*100),"%)\n")
      #ask <- readline(prompt="Press [enter] to continue, Q to quit, p to plot ...")
      #if(ask=="Q") {break}
      #if(ask=="p") {}
      plot_optigs_agp( data = data , hicmap = hicmap , selchr = chr , range = c(s,s+windowsize) , ... ) %>% print
      #plot_optigs_agp( data = data , hicmap = hicmap , selchr = chr , range = c(s,s+windowsize) , obreaks = rye_v1$optig_opbreaks , tbreaks = rye_v1$tenix_breaks)
      agp_sofar <- agp_sofar + stepsize
    }
    #if(ask=="Q") {break}
  }
  dev.off()
}
# debugonce(plot_optigs_agp)
# debugonce(multipage_optig_agp)
# multipage_optig_agp(data = rye_v1$optig_mapping ,  hicmap = rye_v1_hic_map_v1 , fname="agp_optigs_allchrs_breaks.pdf" , obreaks = rye_v1$optig_opbreaks , tbreaks = rye_v1$tenix_breaks)

#integrates with TRITEX assembly objects
#takes a manually written table of superscaffolds and adds a superscaffold object to the assembly
superscaffold <- function( assembly , superscaff_list , species , cores = 21 ){
  
  info <- assembly$info #done ... more to be added in reanchoring step
  optigs <- assembly$optig_mapping #done
  cssaln <- assembly$cssaln #done
  fpairs <- assembly$fpairs #done
  molecules <- assembly$molecules
  if(!is.null(assembly$popseq$scaffold)){
    transpose_popseq <- TRUE
    popseq <- assembly$popseq
  }
  
  cat("Superscaffolding: Munging superscaffold info\n")
  
  #superscaffold object, ssl, will be the bible for superscaffolding throughout, and be attached to the object forthwith. Since we are overwriting some orig_scaffold-to-broken-scaffold conversions, we also save these in the table
  ssl <- copy( superscaff_list )
  ssl <- info[,.(scaffold,scaffold_length = as.numeric(length),orig_scaffold,orig_start)][ssl , on="scaffold" ]
  setorder(ssl,superscaffold,ss_order)
  ssl[ scaffold=="gap" ,scaffold_length := ss_gaplength]
  ssl[, ss_end := cumsum(scaffold_length), by="superscaffold"]
  ssl[ , ss_start := c(1,ss_end[1:(.N-1)]+1) , by="superscaffold"]
  sslen <- ssl[, .( ss_length = sum(scaffold_length) ),by="superscaffold"] #sslen is aggregated info per superscaffold. Currently just a length ... 
  ssl <- sslen[ssl,on="superscaffold"]
  
  
  
  info_same <- info[!scaffold %in% ssl$scaffold][,.(scaffold,length,mr,mri,mr_10x,superscaffold=FALSE)]
  info_change <- sslen[ , .(scaffold=superscaffold,length=ss_length,superscaffold=TRUE,mr=NA_real_,mri=NA_real_,mr_10x=NA_real_) ]
  
  
  new_assembly <- list( ss_info = bind_rows( ssl , info[!scaffold %in% ssl$scaffold][,.(superscaffold = NA , ss_length = length , ss_order = NA , ss_orientation = NA , ss_gaplength = NA , ss_end = length , ss_start = 1 , scaffold , orig_scaffold , orig_start )] ) )
  new_assembly$info <- bind_rows(info_same,info_change)
  
  
  
  
  cat("Superscaffolding: Transpose cssaln\n")
  
  copy(cssaln) -> z
  
  z_same <- z[!scaffold %in% ssl$scaffold]
  z_change <- z[scaffold %in% ssl$scaffold]
  z_change <- ssl[,.(scaffold,superscaffold,ss_start,ss_end,ss_orientation,ss_length)][z_change,on="scaffold"]
  z_change[ , prev_scaffold_length := orig_scaffold_length ]
  z_change[ , orig_scaffold_length := NULL ]
  z_change[ , prev_pos := pos ]
  z_change[ , orig_pos := NULL ]
  z_change[ , pos := ifelse(ss_orientation>0 , (ss_start-1) + prev_pos , (ss_end+1) - prev_pos ) ]
  z_change[ , orig_scaffold := NULL]
  z_change[ , prev_scaffold := scaffold]
  z_change[ , scaffold := superscaffold]
  z_change[ , superscaffold := NULL][ , superscaffold := TRUE]
  
  z_same[ , prev_scaffold := orig_scaffold]
  z_same[ , orig_scaffold := NULL]
  z_same[ , prev_scaffold_length := orig_scaffold_length ]
  z_same[ , orig_scaffold_length := NULL ]
  z_same[ , prev_pos := orig_pos ]
  z_same[ , orig_pos := NULL ]
  z_same[ , superscaffold := FALSE]
  
  z <- bind_rows(z_same , z_change )
  new_assembly$cssaln <- z
  
  
  
  cat("Superscaffolding: Transpose optigs\n")
  z <- copy(optigs)
  z_same <- z[!scaffold %in% ssl$scaffold]
  z_change <- z[scaffold %in% ssl$scaffold]
  
  #stick in superscaffold for scaffold
  
  z_change <- ssl[,.(scaffold,superscaffold,ss_start,ss_end,ss_orientation,ss_length)][z_change,on="scaffold"]
  z_change[, prev_scaffold := scaffold ]
  z_change[, scaffold := superscaffold ]
  z_change[, superscaffold := NULL]
  z_change[ , prev_cmap_length := cmap_length ]
  z_change[ ref_query=="ref" , cmap_length := ss_length ]
  z_change[ , prev_cmap_pos := cmap_pos ]
  z_change[ ref_query=="ref" , cmap_pos := ifelse(ss_orientation>0 , (ss_start-1) + prev_cmap_pos , (ss_end+1) - prev_cmap_pos ) ]
  
  z_same[, prev_scaffold := scaffold ]
  z_same[ , prev_cmap_length := cmap_length ]
  z_same[ , prev_cmap_pos := cmap_pos ]
  
  new_assembly$optig_mapping <- bind_rows( z_change , z_same )
  
  
  
  
  cat("Superscaffolding: Transpose fpairs\n")
  
  z <- copy(fpairs)
  z[,link_idx := 1:.N]
  z1 <- z[,.(scaffold1,chr1,pos1,orig_scaffold1,orig_pos1,link_idx)]
  
  z1_change <- z1[scaffold1 %in% ssl$scaffold][,scaffold := scaffold1]
  z1_change <- ssl[,.(scaffold,superscaffold,ss_start,ss_end,ss_orientation)][z1_change,on="scaffold"]
  z1_change[,prev_scaffold1 := scaffold1]
  z1_change[,scaffold1 := superscaffold]
  z1_change[,prev_pos1 := pos1 ]
  z1_change[,pos1 := ifelse(ss_orientation>0 , (ss_start-1) + prev_pos1 , (ss_end+1) - prev_pos1 )]
  z1_change[,superscaffold := NULL]
  z1_change <- z1_change[,.(scaffold1,pos1,orig_scaffold1,orig_pos1,prev_scaffold1,prev_pos1,link_idx)]
  #
  z1_same <- z1[!scaffold1 %in% ssl$scaffold]
  z1_same[,prev_scaffold1 := scaffold1]
  z1_same[,prev_pos1 := pos1 ]
  z1_same <- z1_same[,.(scaffold1,pos1,orig_scaffold1,orig_pos1,prev_scaffold1,prev_pos1,link_idx)]
  
  z1 <- bind_rows(z1_same,z1_change)
  
  
  z2 <- z[,.(scaffold2,chr2,pos2,orig_scaffold2,orig_pos2,link_idx)]
  z2_change <- z2[scaffold2 %in% ssl$scaffold][,scaffold := scaffold2]
  z2_change <- ssl[,.(scaffold,superscaffold,ss_start,ss_end,ss_orientation)][z2_change,on="scaffold"]
  z2_change[,prev_scaffold2 := scaffold2]
  z2_change[,scaffold2 := superscaffold]
  z2_change[,prev_pos2 := pos2 ]
  z2_change[,pos2 := ifelse(ss_orientation>0 , (ss_start-1) + prev_pos2 , (ss_end+1) - prev_pos2 )]
  z2_change[,superscaffold := NULL]
  z2_change <- z2_change[,.(scaffold2,pos2,orig_scaffold2,orig_pos2,prev_scaffold2,prev_pos2,link_idx)]
  #
  z2_same <- z2[!scaffold2 %in% ssl$scaffold]
  z2_same[,prev_scaffold2 := scaffold2]
  z2_same[,prev_pos2 := pos2 ]
  z2_same <- z2_same[,.(scaffold2,pos2,orig_scaffold2,orig_pos2,prev_scaffold2,prev_pos2,link_idx)]
  
  z2 <- bind_rows(z2_same,z2_change)
  
  
  z <- z1[z2,on="link_idx"]
  

  new_assembly$fpairs <- z
  
  
  
  
  cat("Superscaffolding: Transpose molecules\n")
  
  z <- copy(molecules)
  
  z_same <- z[!scaffold %in% ssl$scaffold]
  z_change <- z[scaffold %in% ssl$scaffold]
  
  z_change <- ssl[,.(scaffold,superscaffold,ss_start,ss_end,ss_orientation)][z_change,on="scaffold"]
  
  z_change[,prev_scaffold := scaffold]
  z_change[,scaffold := superscaffold]
  z_change[,prev_start := start]
  z_change[,prev_end := end]
  z_change[,start := ifelse(ss_orientation>0 , (ss_start-1) + prev_start , (ss_end+1) - prev_end )]
  z_change[,end := ifelse(ss_orientation>0 , (ss_start-1) + prev_end , (ss_end+1) - prev_start )]
  z_change <- z_change[,.(orig_scaffold,scaffold,start,end,barcode,npairs,sample,length,orig_start,orig_end,prev_scaffold,prev_start,prev_end)]
  
  z_same[,prev_scaffold := scaffold]
  z_same[,prev_start := start]
  z_same[,prev_end := end]
  
  z <- bind_rows(z_same,z_change)
  
  new_assembly$molecules <- z
  
  if (transpose_popseq) {
    cat("Superscaffolding: Looks like pseudo-popseq with real css; transposing 'popseq' object\n")
    
    z <- copy(popseq)
    
    z_change <- z[scaffold %in% ssl$scaffold]
    if(nrow(z_change) > 1){
      z_change <- ssl[,.(scaffold,superscaffold,ss_start,ss_end,ss_orientation,ss_length)][z_change,on="scaffold"]
      z_change[,prev_scaffold := scaffold ]
      z_change[,scaffold := superscaffold]
      z_change[,prev_scaffold_length := scaffold_length ]
      z_change[,scaffold_length := ss_length]
      z_change[,prev_pos := pos ]
      z_change[,pos := ifelse(ss_orientation>0 , (ss_start-1) + prev_pos , (ss_end+1) - prev_pos )]
      z_change[,superscaffold := NULL]
      z_change <- z_change[,.(scaffold_length,scaffold,css_contig,pos,popseq_alphachr,popseq_chr,popseq_cM,css_contig_length,prev_scaffold,prev_scaffold_length,prev_pos,superscaffold = TRUE)]
    }
    
    z_same <- z[!scaffold %in% ssl$scaffold]
    z_same[,prev_scaffold := scaffold ]
    z_same[,prev_scaffold_length := scaffold_length ]
    z_same[,prev_pos := pos ]
    z_same[,superscaffold := FALSE]
    
    z <- bind_rows( z_same , z_change )
    
    new_assembly$popseq <- z
    
    
  } else {
    stop("Superscaffolding function not (yet?) implemented for popseq data reliant on cssaln for positional info (as inherited from MMs pipeline).\n")
  }
  
  
  
  cat("Superscaffolding: (Re)anchor scaffolds\n")
  anchor_scaffolds(new_assembly, popseq=new_assembly$popseq, species=species ) -> new_assembly
  
  
  #recalc hic_coverage, 10x cov, mr_10x, mr, (see break routine end)
  
  info <- copy(new_assembly$info)
  
  if("mr_10x" %in% names(new_assembly$info)){
    new_assembly$info[, mr_10x := NULL]
  }
  
  if("mr" %in% names(new_assembly$info)){
    new_assembly$info[, mr := NULL]
    new_assembly$info[, mri := NULL]
  }
  
  if("cov" %in% names(assembly) & nrow(fpairs) > 0){
    cat("Superscaffolding: (Re)calculating Hi-C coverage\n")
    cov_change <- add_hic_cov(new_assembly, scaffolds=new_assembly$info[superscaffold==TRUE]$scaffold, binsize=assembly$binsize, minNbin=assembly$minNbin, innerDist=assembly$innerDist, cores=cores)
    cov_same <- assembly$cov[!scaffold %in% new_assembly$info[superscaffold==TRUE]$scaffold]
    if(nrow(cov_change$cov) > 0){
      new_assembly$cov <- bind_rows( cov_same , cov_change$cov )
    } else {
      new_assembly$cov <- cov_same
    }
    
    info_same <- new_assembly$info[superscaffold==FALSE]
    info_change <- cov_change$info
    new_assembly$info <- bind_rows( info_same , info_change )
    
  } else {
    assembly_new$cov <- data.table()
  }
  
  if("molecule_cov" %in% names(assembly) & nrow(molecules) > 0){
    cat("Superscaffolding: (Re)calculating 10X molecule coverage\n")
    mcov_change <- add_molecule_cov(new_assembly, scaffolds=info[superscaffold==TRUE]$scaffold, binsize=assembly$mol_binsize, cores=cores)
    
    info_same <- info[superscaffold==FALSE]
    new_assembly$info <- bind_rows( info_same[, names(mcov_change$info), with=F], mcov_change$info)
    
    new_assembly$mol_binsize <- assembly$mol_binsize
    
    
    mcov_same <- assembly$molecule_cov[!info[superscaffold==TRUE]$scaffold, on="scaffold"]
    if(nrow(mcov_change$molecule_cov) > 0){
      new_assembly$molecule_cov <- rbind(mcov_same, mcov_change$molecule_cov)
    } else {
      new_assembly$molecule_cov <- mcov_same
    }
  } else {
    new_assembly$molecule_cov <- data.table()
  }
  
  new_assembly$binsize <- assembly$binsize
  new_assembly$innerDist <- assembly$innerDist
  new_assembly$minNbin <- assembly$minNbin
  
  
  
  new_assembly
}
#superscaffold( assembly = rye_v4t , superscaff_list = sscl , species = "rye" )