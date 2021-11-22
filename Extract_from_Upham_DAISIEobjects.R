## Extract colonization and branching times from the Malagasy clades
## from the posterior distribution of trees of Upham et al.2019
# Script by Luis Valente

rm(list=ls())

library(ape)
library(phytools)
library(stringr)


# # "Main" dataset settings
# #Upham_tree_type<-"DNA_only"
# ## SET COLONISATION SCENARIO 1 or 2
# Colonisation_scenario<-1
# # SET EXTINCTION SCENARIO "low_impact" or "high_impact"
# Extinction_scenario<-'high_impact'
# ## Set tree type "posterior" or "maxcred"
# posterior_or_maxcred<-"maxcred"
# # Set ISLAND AGE
# island_age<-88
# # Set MAINLAND POOL SIZE
# M<-1000
# ## Proportion of bat species in the mainland pool
# prop_type2_pool<-0.22


#####  SETTINGS
# SET WHICH TREE TYPE TO USE: "DNA_only" or "Complete" trees
Upham_tree_type<-"DNA_only"
## SET COLONISATION SCENARIO 1 or 2
Colonisation_scenario<-1
# SET EXTINCTION SCENARIO "low_impact" or "high_impact"
Extinction_scenario<-'high_impact'
## Set tree type "posterior" or "maxcred"
posterior_or_maxcred<-"maxcred"
# Set ISLAND AGE
island_age<-88
# Set MAINLAND POOL SIZE
M<-1000
## Proportion of bat species in the mainland pool
prop_type2_pool<-0.22


## Load the mammal trees
if(posterior_or_maxcred=="posterior"){
  if(Upham_tree_type=="Complete")
  {mtrees<-read.nexus("Upham_1000_ND_Completed5911.trees")}
  if(Upham_tree_type=="DNA_only")
  {mtrees<-read.nexus("Upham_1000_ND_DNAonly4098.trees")}
}

if(posterior_or_maxcred=="maxcred"){
  if(Upham_tree_type=="Complete")
  {mtre_save<-read.nexus("MaxCred_Upham_VertLife_ND_Completed5911.tre")}
  if(Upham_tree_type=="DNA_only")
  {mtre_save<-read.nexus("MaxCred_Upham_VertLife_ND_DNAonly4098.tre")}
  
  mtrees<-list()
  mtrees[[1]]<-mtre_save
  
}

## Function to extract tip ages for single taxa
get_leaf_age<-function(tree,leaf_name){
  the_node<-which(tree$tip.label==leaf_name)
  return(mtree$edge.length[which(mtree$edge[,2]==the_node)])}

## Function to extract stem and branching times for cases where closest relative to the clade
## can vary from tree to tree
## and for cases where the ingroup is not monophyletic in 100% of the trees
## Code for this function modified from code available online by Liam Revell (http://blog.phytools.org/)
get_stem_and_branching_times<-function(tree,ingroup) {
  H<-nodeHeights(tree)
  H<-max(H)-H
  subclade<-keep.tip(tree,ingroup)
  crown_btimes<-as.vector(branching.times(subclade))
  crown_btimes<-sort(crown_btimes,decreasing = T)
  crown_node<-findMRCA(tree,ingroup)
  stem<-H[tree$edge==phytools:::getAncestors(tree,crown_node,"parent")][1]
  btimes<-c(stem,crown_btimes)
  return(btimes)
}



## START Loop to go through the trees to extract times

btimes_list<-list()

for (i in 1:length(mtrees)) {
  
  mtree<-mtrees[[i]]
  
  
  ## Singleton colonisations (BATS) present in Upham et al. complete
  btimes_Chaerephon_atsinanana<-get_leaf_age(mtree,"Chaerephon_atsinanana_MOLOSSIDAE_CHIROPTERA")
  btimes_Coleura_kibomalandy<-get_leaf_age(mtree,"Coleura_kibomalandy_EMBALLONURIDAE_CHIROPTERA")
  btimes_Eidolon_dupreanum<-get_leaf_age(mtree,"Eidolon_dupreanum_PTEROPODIDAE_CHIROPTERA")
  btimes_Laephotis_matroka<-get_leaf_age(mtree,"Neoromicia_matroka_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Mops_leucostigma<-get_leaf_age(mtree,"Mops_leucostigma_MOLOSSIDAE_CHIROPTERA")
  btimes_Mops_midas<-get_leaf_age(mtree,"Mops_midas_MOLOSSIDAE_CHIROPTERA")
  btimes_Mormopterus_jugularis<-get_leaf_age(mtree,"Mormopterus_jugularis_MOLOSSIDAE_CHIROPTERA")
  btimes_Myotis_goudoti<-get_leaf_age(mtree,"Myotis_goudoti_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Neoromicia_bemainty<-get_leaf_age(mtree,"Hypsugo_bemainty_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Otomops_madagascariensis<-get_leaf_age(mtree,"Otomops_madagascariensis_MOLOSSIDAE_CHIROPTERA")
  btimes_Pipistrellus_hesperidus<-get_leaf_age(mtree,"Pipistrellus_hesperidus_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Pipistrellus_raceyi<-get_leaf_age(mtree,"Pipistrellus_raceyi_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Pteropus_rufus<-get_leaf_age(mtree,"Pteropus_rufus_PTEROPODIDAE_CHIROPTERA")
  btimes_Rousettus_madagascariensis<-get_leaf_age(mtree,"Rousettus_madagascariensis_PTEROPODIDAE_CHIROPTERA")
  btimes_Scotophilus_robustus<-get_leaf_age(mtree,"Scotophilus_robustus_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Tadarida_fulminans<-get_leaf_age(mtree,"Tadarida_fulminans_MOLOSSIDAE_CHIROPTERA")
  btimes_Chaerephon_jobimena<-get_leaf_age(mtree,"Tadarida_jobimena_MOLOSSIDAE_CHIROPTERA")
  btimes_Taphozous_mauritianus<-get_leaf_age(mtree,"Taphozous_mauritianus_EMBALLONURIDAE_CHIROPTERA")
  btimes_Triaenops_menamena<-get_leaf_age(mtree,"Triaenops_menamena_HIPPOSIDERIDAE_CHIROPTERA")
  
  ## Singleton bat colonisations not present in Upham complete nor DNA_only
  ## Chaerephon_leucogaster absent from both, so we use crown age of the family as MaxAge
  ingroup<-str_subset(mtree$tip.label,'MOLOSSIDAE')
  ## Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA (South America) comes out in an odd place on the tree,
  # so it is removed, otherwise we would get the incorrect crown age for the family
  ingroup<-ingroup[-which(ingroup=='Nyctinomops_kalinowskii_MOLOSSIDAE_CHIROPTERA')]
  btimes_Chaerephon_leucogaster<-get_stem_and_branching_times(mtree,ingroup)
  btimes_Chaerephon_leucogaster<-btimes_Chaerephon_leucogaster[2]
  
  
  ## Singleton bat colonisation not present in Upham DNA_only tree (but present in Complete)
  if(Upham_tree_type=="DNA_only")
  {
    ## Pipistrellus raceyi absent from DNA_only tree so we use crown age of the family Vespertilionidae
    # as MaxAge (Pippistrellus genus is not monophyletic on the DNA only trees)
    ingroup<-str_subset(mtree$tip.label,'VESPERTILIONIDAE')
    btimes_Pipistrellus_raceyi<-get_stem_and_branching_times(mtree,ingroup)
    btimes_Pipistrellus_raceyi<-btimes_Pipistrellus_raceyi[2] }
  
  #### Bat clades with more than one species
  
  ## Laephotis_malagasyensis_robertsi
  ingroup<-c("Neoromicia_malagasyensis_VESPERTILIONIDAE_CHIROPTERA","Neoromicia_robertsi_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Laephotis_malagasyensis_robertsi<-get_stem_and_branching_times(mtree,ingroup)
  
  
  ### Macronycteris
  # Macronycteris Scenario 1 - two colonisations (besaoka+commersoni together, besaoka added as missing to commersoni btimes)
  if(Colonisation_scenario==1){
    btimes_Macronycteris_besaoka_commersoni<-get_leaf_age(mtree,"Hipposideros_commersoni_HIPPOSIDERIDAE_CHIROPTERA")
    ## Macronycteris cryptovalorona absent from both trees, so we use crown age of the family as MaxAge
    ingroup<-str_subset(mtree$tip.label,'HIPPOSIDERIDAE')
    btimes_Macronycteris_cryptovalorona<-get_stem_and_branching_times(mtree,ingroup)
    btimes_Macronycteris_cryptovalorona<-btimes_Macronycteris_cryptovalorona[2]
  }
  # Macronycteris Scenario 2 - three colonisations
  if(Colonisation_scenario==2){
    btimes_Macronycteris_commersoni<-get_leaf_age(mtree,"Hipposideros_commersoni_HIPPOSIDERIDAE_CHIROPTERA")
    
    ## Macronycteris besaoka absent from both trees, so we use crown age of the family as MaxAge
    ingroup<-str_subset(mtree$tip.label,'HIPPOSIDERIDAE')
    btimes_Macronycteris_besaoka<-get_stem_and_branching_times(mtree,ingroup)
    btimes_Macronycteris_besaoka<-btimes_Macronycteris_besaoka[2]
    
    ## Macronycteris cryptovalorona absent from both trees, so we use crown age of the family as MaxAge (same as above)
    ingroup<-str_subset(mtree$tip.label,'HIPPOSIDERIDAE')
    btimes_Macronycteris_cryptovalorona<-get_stem_and_branching_times(mtree,ingroup)
    btimes_Macronycteris_cryptovalorona<-btimes_Macronycteris_cryptovalorona[2] }
  
  
  ## Miniopterus
  ### Resolution is very bad for this radiation so it is not appropriate to extract branching times,
  ## instead we use the Crown age of the genus as a max colonisation time and add all species as missing
  ## We chose crown and not stem because Miniopterus of Madagascar are not basal in the genus
  # Miniopterus Scenario 1 - One single colonisation
  # Miniopterus Scenario 2 - Four colonisations
  if(Colonisation_scenario==1){
    ingroup<-str_subset(mtree$tip.label,'Miniopterus')
    subclade<-keep.tip(mtree,ingroup)
    btimes_Miniopterus_all<-max(branching.times(subclade)) }
  
  # Miniopterus Scenario 2 - four colonisations
  if(Colonisation_scenario==2){
    
    ingroup<-str_subset(mtree$tip.label,'Miniopterus')
    subclade<-keep.tip(mtree,ingroup)
    btimes_Miniopterus_1<-max(branching.times(subclade))
    
    btimes_Miniopterus_griveaudi<-get_leaf_age(mtree,"Miniopterus_griveaudi_VESPERTILIONIDAE_CHIROPTERA")
    btimes_Miniopterus_mahafaliensis<-get_leaf_age(mtree,"Miniopterus_mahafaliensis_VESPERTILIONIDAE_CHIROPTERA")
    btimes_Miniopterus_sororculus<-get_leaf_age(mtree,"Miniopterus_sororculus_VESPERTILIONIDAE_CHIROPTERA")
    
    if(Upham_tree_type=="DNA_only") {
      ## Use crown age of Miniopterus as MaxAge
      btimes_Miniopterus_mahafaliensis<-btimes_Miniopterus_1
      btimes_Miniopterus_sororculus<-btimes_Miniopterus_1}
    
  }
  
  
  ## Myzopoda
  ingroup<-str_subset(mtree$tip.label,'Myzopoda')
  btimes_Myzopoda<-get_stem_and_branching_times(mtree,ingroup)
  
  ## Paremballonura
  ingroup<-c("Emballonura_tiavato_EMBALLONURIDAE_CHIROPTERA","Emballonura_atrata_EMBALLONURIDAE_CHIROPTERA")
  btimes_Paremballonura<-get_stem_and_branching_times(mtree,ingroup)
  
  ## Paratriaenops
  ## This assumes that Paratriaenops pauliani (endemic to Aldabra) colonised Aldabra from Madagascar
  ingroup<-c("Paratriaenops_auritus_HIPPOSIDERIDAE_CHIROPTERA",
             "Paratriaenops_furculus_HIPPOSIDERIDAE_CHIROPTERA")
  btimes_Paratriaenops<-get_stem_and_branching_times(mtree,ingroup)
  
  ## Scotophilus_maro+tandr
  ingroup<-c("Scotophilus_marovaza_VESPERTILIONIDAE_CHIROPTERA",
             "Scotophilus_tandrefana_VESPERTILIONIDAE_CHIROPTERA")
  btimes_Scotophilus_maro_tandr<-get_stem_and_branching_times(mtree,ingroup)
  
  
  ##### LAND MAMMALS
  
  ### CARNIVORA/ EUPLERIDAE
  ingroup<-str_subset(mtree$tip.label,'EUPLERIDAE')
  ## remove Galidictis grandidieri as it is now considered subspecies of fasciata
  # (this taxon is not present in the DNA-only tree anyway)
  ingroup<-str_subset(ingroup,"Galidictis_grandidieri_EUPLERIDAE_CARNIVORA",negate = T)
  
  ## remove Salanoia_durrelli as it is has been synonimized with Salanoia_concolor
  ingroup<-str_subset(ingroup,"Salanoia_durrelli_EUPLERIDAE_CARNIVORA",negate = T)
  
  
  if(Upham_tree_type=='Complete'){
    if(Extinction_scenario=='low_impact'){
      # From Complete trees, low human impact, branching time of Cryptoprocta spelea should be removed
      ingroup<-str_subset(ingroup,"Cryptoprocta_spelea_EUPLERIDAE_CARNIVORA",negate = T)}
  }
  btimes_Eupleridae<-get_stem_and_branching_times(mtree,ingroup)
  
  
  ### TENRECS
  #### NOTE: 3 species of otter shrews were assigned TENRECIDAE in the tip labels
  ## of Upham et al, but are now considered to be a different family
  # and are not found in Madagascar,
  ## so these names were first changed in all trees to remove
  # TENRECIDAE from their tip labels:
  # "Potamogale_velox_AFROSORICIDA
  # "Micropotamogale_ruwenzorii_AFROSORICIDA"
  # "Micropotamogale_lamottei_AFROSORICIDA"
  ## Tenrec + Bibymalagasy Scenario 1 - Single colonisation
  if(Colonisation_scenario==1){
    ingroup<-str_subset(mtree$tip.label,'TENRECIDAE')
    btimes_Afrosoricida<-get_stem_and_branching_times(mtree,ingroup) }
  
  ## Tenrec + Bibymalagasy Scenario 2 - Bibymalagasy separate colonisation of unknown time
  if(Colonisation_scenario==2){
    ingroup<-str_subset(mtree$tip.label,'TENRECIDAE')
    btimes_Tenrecs<-get_stem_and_branching_times(mtree,ingroup)
    
    ## For Bibymalagasia use stem age of Afrosoricida as MaxAge
    ingroup<-str_subset(mtree$tip.label,'AFROSORICIDA')
    btimes_Bibymalagasia<-get_stem_and_branching_times(mtree,ingroup)
    btimes_Bibymalagasia<-btimes_Bibymalagasia[1]
  }
  
  
  ### LEMURS
  ingroup<-str_subset(mtree$tip.label,'LEMURS')
  btimes_Lemurs<-get_stem_and_branching_times(mtree,ingroup)
  
  
  ### NESOMYIDAE
  ## All the names are typed here here because NESOMYIDAE also includes a large group of species
  ## from outside Madagascar so string-matching NESOMYIDAE against Upham would not work
  if(Upham_tree_type=="DNA_only")
  {
    ingroup<-c("Brachytarsomys_albicauda_NESOMYIDAE_RODENTIA",
               "Brachyuromys_betsileoensis_NESOMYIDAE_RODENTIA",
               "Brachyuromys_ramirohitra_NESOMYIDAE_RODENTIA",
               "Eliurus_antsingy_NESOMYIDAE_RODENTIA",
               "Eliurus_carletoni_NESOMYIDAE_RODENTIA",
               "Eliurus_ellermani_NESOMYIDAE_RODENTIA",
               "Eliurus_grandidieri_NESOMYIDAE_RODENTIA",
               "Eliurus_majori_NESOMYIDAE_RODENTIA",
               "Eliurus_minor_NESOMYIDAE_RODENTIA",
               "Eliurus_myoxinus_NESOMYIDAE_RODENTIA",
               "Eliurus_tanala_NESOMYIDAE_RODENTIA",
               "Eliurus_webbi_NESOMYIDAE_RODENTIA",
               "Gymnuromys_roberti_NESOMYIDAE_RODENTIA",
               "Hypogeomys_antimena_NESOMYIDAE_RODENTIA",
               "Macrotarsomys_bastardi_NESOMYIDAE_RODENTIA",
               "Macrotarsomys_ingens_NESOMYIDAE_RODENTIA",
               "Monticolomys_koopmani_NESOMYIDAE_RODENTIA",
               "Nesomys_audeberti_NESOMYIDAE_RODENTIA",
               "Nesomys_rufus_NESOMYIDAE_RODENTIA",
               "Voalavo_gymnocaudus_NESOMYIDAE_RODENTIA")
  }
  
  if(Upham_tree_type== "Complete")
  {
    ingroup<-c("Brachytarsomys_albicauda_NESOMYIDAE_RODENTIA",
               "Brachytarsomys_villosa_NESOMYIDAE_RODENTIA",
               "Brachyuromys_betsileoensis_NESOMYIDAE_RODENTIA",
               "Brachyuromys_ramirohitra_NESOMYIDAE_RODENTIA",
               "Eliurus_antsingy_NESOMYIDAE_RODENTIA",
               "Eliurus_carletoni_NESOMYIDAE_RODENTIA",
               "Eliurus_danieli_NESOMYIDAE_RODENTIA",
               "Eliurus_ellermani_NESOMYIDAE_RODENTIA",
               "Eliurus_grandidieri_NESOMYIDAE_RODENTIA",
               "Eliurus_majori_NESOMYIDAE_RODENTIA",
               "Eliurus_minor_NESOMYIDAE_RODENTIA",
               "Eliurus_myoxinus_NESOMYIDAE_RODENTIA",
               "Eliurus_penicillatus_NESOMYIDAE_RODENTIA",
               "Eliurus_petteri_NESOMYIDAE_RODENTIA",
               "Eliurus_tanala_NESOMYIDAE_RODENTIA",
               "Eliurus_webbi_NESOMYIDAE_RODENTIA",
               "Gymnuromys_roberti_NESOMYIDAE_RODENTIA",
               "Hypogeomys_antimena_NESOMYIDAE_RODENTIA",
               "Macrotarsomys_bastardi_NESOMYIDAE_RODENTIA",
               "Macrotarsomys_ingens_NESOMYIDAE_RODENTIA",
               "Macrotarsomys_petteri_NESOMYIDAE_RODENTIA",
               "Monticolomys_koopmani_NESOMYIDAE_RODENTIA",
               "Nesomys_audeberti_NESOMYIDAE_RODENTIA",
               "Nesomys_lambertoni_NESOMYIDAE_RODENTIA",
               "Nesomys_rufus_NESOMYIDAE_RODENTIA",
               "Voalavo_antsahabensis_NESOMYIDAE_RODENTIA",
               "Voalavo_gymnocaudus_NESOMYIDAE_RODENTIA")
  }
  ## Some non-Madagascar species often appear within the Madagascar clade,
  ## including in the MCC tree from Upham, but with extremely low support:
  ## Dendroprionomys_rousseloti, Prionomys_batesi, Megadendromus_nikolausi
  ## the following code already excludes their branching times in case they are
  ## retrieved within the clade for a given tree
  btimes_Nesomyinae<-get_stem_and_branching_times(mtree,ingroup)
  
  
  
  
  ## Hippos
  ### Resolution is  bad for hippos, so it is best not to extract branching times,
  ## instead we use the stem age of the genus Hippopotamus as a max colonisation time and add all species as missing
  ## Hippos Scenario 1 - single colonisation
  if(Colonisation_scenario==1){
    ingroup<-str_subset(mtree$tip.label,'Hippopotamus')
    btimes_Hippos<-get_stem_and_branching_times(mtree,ingroup)
    btimes_Hippos<-btimes_Hippos[1]
    
    ## None of the Madagascar hippos are in the DNA only tree, so use age of hippo genus as max
    if(Upham_tree_type=="DNA_only") {
      btimes_Hippos<-get_leaf_age(mtree,"Hippopotamus_amphibius_HIPPOPOTAMIDAE_CETARTIODACTYLA")}
    
  }
  
  ## Hippos Scenario 2 - two colonisations, with laloumena separate
  if(Colonisation_scenario==2){
    ingroup<-str_subset(mtree$tip.label,'Hippopotamus')
    btimes_Hipp<-get_stem_and_branching_times(mtree,ingroup)
    # repeat 2 types because max age is the same for both independent colonisations as it is the stem of Hippopotamus
    btimes_Hippos_mad_lem<-btimes_Hipp[1]
    btimes_Hippos_laloumena<-btimes_Hipp[1]
    
    if(Upham_tree_type=="DNA_only") {
      btimes_Hippos_mad_lem<-get_leaf_age(mtree,"Hippopotamus_amphibius_HIPPOPOTAMIDAE_CETARTIODACTYLA")
      btimes_Hippos_laloumena<-get_leaf_age(mtree,"Hippopotamus_amphibius_HIPPOPOTAMIDAE_CETARTIODACTYLA")}
    
  }
  
  
  
  
  ## COLONISATION SCENARIO 1
  if(Colonisation_scenario==1){
    btimes_list[[i]]<-
      list(
        Afrosoricida=btimes_Afrosoricida,
        Chaerephon_atsinanana=btimes_Chaerephon_atsinanana,
        Chaerephon_jobimena=btimes_Chaerephon_jobimena,
        Chaerephon_leucogaster=btimes_Chaerephon_leucogaster,
        Coleura_kibomalandy=btimes_Coleura_kibomalandy,
        Eidolon_dupreanum=btimes_Eidolon_dupreanum,
        Eupleridae=btimes_Eupleridae,
        Hippos=btimes_Hippos,
        Laephotis_matroka=btimes_Laephotis_matroka,
        Laephotis_malagasyensis_robertsi=btimes_Laephotis_malagasyensis_robertsi,
        Lemurs=btimes_Lemurs,
        Macronycteris_besaoka_commersoni=btimes_Macronycteris_besaoka_commersoni,
        Macronycteris_cryptovalorona=btimes_Macronycteris_cryptovalorona,
        Miniopterus_all=btimes_Miniopterus_all,
        Mops_leucostigma=btimes_Mops_leucostigma,
        Mops_midas=btimes_Mops_midas,
        Mormopterus_jugularis=btimes_Mormopterus_jugularis,
        Myotis_goudoti=btimes_Myotis_goudoti,
        Myzopoda=btimes_Myzopoda,
        Neoromicia_bemainty=btimes_Neoromicia_bemainty,
        Nesomyinae=btimes_Nesomyinae,
        Otomops_madagascariensis=btimes_Otomops_madagascariensis,
        Paratriaenops=btimes_Paratriaenops,
        Paremballonura=btimes_Paremballonura,
        Pipistrellus_hesperidus=btimes_Pipistrellus_hesperidus,
        Pipistrellus_raceyi=btimes_Pipistrellus_raceyi,
        Pteropus_rufus=btimes_Pteropus_rufus,
        Rousettus_madagascariensis=btimes_Rousettus_madagascariensis,
        Scotophilus_maro_tandr=btimes_Scotophilus_maro_tandr,
        Scotophilus_robustus=btimes_Scotophilus_robustus,
        Tadarida_fulminans=btimes_Tadarida_fulminans,
        Taphozous_mauritianus=btimes_Taphozous_mauritianus,
        Triaenops_menamena=btimes_Triaenops_menamena)}
  
  ## COLONISATION SCENARIO 2
  if(Colonisation_scenario==2){
    btimes_list[[i]]<-
      list(
        Bibymalagasia=btimes_Bibymalagasia,
        Chaerephon_atsinanana=btimes_Chaerephon_atsinanana,
        Chaerephon_jobimena=btimes_Chaerephon_jobimena,
        Chaerephon_leucogaster=btimes_Chaerephon_leucogaster,
        Coleura_kibomalandy=btimes_Coleura_kibomalandy,
        Eidolon_dupreanum=btimes_Eidolon_dupreanum,
        Eupleridae=btimes_Eupleridae,
        Hippos_mad_lem=btimes_Hippos_mad_lem,
        Hippos_laloumena=btimes_Hippos_laloumena,
        Laephotis_matroka=btimes_Laephotis_matroka,
        Laephotis_malagasyensis_robertsi=btimes_Laephotis_malagasyensis_robertsi,
        Lemurs=btimes_Lemurs,
        Macronycteris_besaoka=btimes_Macronycteris_besaoka,
        Macronycteris_commersoni=btimes_Macronycteris_commersoni,
        Macronycteris_cryptovalorona=btimes_Macronycteris_cryptovalorona,
        Miniopterus_1=btimes_Miniopterus_1,
        Miniopterus_griveaudi=btimes_Miniopterus_griveaudi,
        Miniopterus_mahafaliensis=btimes_Miniopterus_mahafaliensis,
        Miniopterus_sororculus=btimes_Miniopterus_sororculus,
        Mops_leucostigma=btimes_Mops_leucostigma,
        Mops_midas=btimes_Mops_midas,
        Mormopterus_jugularis=btimes_Mormopterus_jugularis,
        Myotis_goudoti=btimes_Myotis_goudoti,
        Myzopoda=btimes_Myzopoda,
        Neoromicia_bemainty=btimes_Neoromicia_bemainty,
        Nesomyinae=btimes_Nesomyinae,
        Otomops_madagascariensis=btimes_Otomops_madagascariensis,
        Paratriaenops=btimes_Paratriaenops,
        Paremballonura=btimes_Paremballonura,
        Pipistrellus_hesperidus=btimes_Pipistrellus_hesperidus,
        Pipistrellus_raceyi=btimes_Pipistrellus_raceyi,
        Pteropus_rufus=btimes_Pteropus_rufus,
        Rousettus_madagascariensis=btimes_Rousettus_madagascariensis,
        Scotophilus_maro_tandr=btimes_Scotophilus_maro_tandr,
        Scotophilus_robustus=btimes_Scotophilus_robustus,
        Tadarida_fulminans=btimes_Tadarida_fulminans,
        Taphozous_mauritianus=btimes_Taphozous_mauritianus,
        Tenrecs=btimes_Tenrecs,
        Triaenops_menamena=btimes_Triaenops_menamena)}
  
  
  print(paste("Extracted from tree",i,"out of", length(mtrees)))
}

library(DAISIE)
library(readxl)

## LOAD btimes_list and select type2 species

if(Colonisation_scenario==1){
  ## load table with Clade names, Status and Missing_species for each clade
  mtable<-read_excel("/Users/luis/Desktop/MADAGASCAR/PHYLOGENETIC_DATA/Colonisations_Scenario1.xlsx")
  
  
  
  type2_bat_name_list<-c("Chaerephon_atsinanana",
                         "Chaerephon_jobimena",
                         "Chaerephon_leucogaster",
                         "Coleura_kibomalandy",
                         "Eidolon_dupreanum",
                         "Laephotis_matroka",
                         "Macronycteris_besaoka_commersoni",
                         "Macronycteris_cryptovalorona",
                         "Mops_leucostigma",
                         "Mops_midas",
                         "Mormopterus_jugularis",
                         "Myotis_goudoti",
                         "Neoromicia_bemainty",
                         "Otomops_madagascariensis",
                         "Pipistrellus_hesperidus",
                         "Pipistrellus_raceyi",
                         "Pteropus_rufus",
                         "Rousettus_madagascariensis",
                         "Scotophilus_robustus",
                         "Tadarida_fulminans",
                         "Taphozous_mauritianus",
                         "Triaenops_menamena",
                         "Laephotis_malagasyensis_robertsi",
                         "Myzopoda",
                         "Paremballonura",
                         "Scotophilus_maro_tandr",
                         "Paratriaenops",
                         "Miniopterus_all")
  
}


if(Colonisation_scenario==2){
  ## load table with Clade names, Status and Missing_species for each clade
  mtable<-read_excel("/Users/luis/Desktop/MADAGASCAR/PHYLOGENETIC_DATA/Colonisations_Scenario2.xlsx")
  
  if(Extinction_scenario=='low_impact'){
    mtable<-mtable[-which(mtable[,'Clade_name']=='Bibymalagasia'),]
    mtable<-mtable[-which(mtable[,'Clade_name']=='Macronycteris_besaoka'),]
    mtable<-mtable[-which(mtable[,'Clade_name']=='Hippos_laloumena'),]
  }
  
  
  type2_bat_name_list<-c("Chaerephon_atsinanana",
                         "Chaerephon_jobimena",
                         "Chaerephon_leucogaster",
                         "Coleura_kibomalandy",
                         "Eidolon_dupreanum",
                         "Laephotis_matroka",
                         "Macronycteris_besaoka",
                         "Macronycteris_commersoni",
                         "Macronycteris_cryptovalorona",
                         "Miniopterus_griveaudi",
                         "Miniopterus_mahafaliensis",
                         "Miniopterus_sororculus",
                         "Mops_leucostigma",
                         "Mops_midas",
                         "Mormopterus_jugularis",
                         "Myotis_goudoti",
                         "Neoromicia_bemainty",
                         "Otomops_madagascariensis",
                         "Pipistrellus_hesperidus",
                         "Pipistrellus_raceyi",
                         "Pteropus_rufus",
                         "Rousettus_madagascariensis",
                         "Scotophilus_robustus",
                         "Tadarida_fulminans",
                         "Taphozous_mauritianus",
                         "Triaenops_menamena",
                         "Laephotis_malagasyensis_robertsi",
                         "Myzopoda",
                         "Paremballonura",
                         "Scotophilus_maro_tandr",
                         "Paratriaenops",
                         "Miniopterus_1")
  
  # Remove Macronycteris besaoka from bat list for low human impact scenario
  if(Extinction_scenario=='low_impact'){
    type2_bat_name_list<-type2_bat_name_list[-7]
  }
  
  
  
}

## Find order of btimes to match same order as the clade names in the data prep table

clade_ordering<-match(mtable$Clade_name,names(btimes_list[[1]]))


## TYPE 1 ONLY
## Loop through the different trees to create separate DAISIE data objects
Madagascar_list<-list()
for (i in 1:length(btimes_list)){
  ##  order the list in case the table and btimes are not in the same order
  ordered_btimes_list<-btimes_list[[i]][clade_ordering]
  ## convert the branching times into a factor that can be added to a data frame
  btimes_as_factor =  data.frame(col1 = unlist(lapply(ordered_btimes_list, paste, collapse = ", ")))
  
  if(Upham_tree_type=="Complete") { 
    if(Extinction_scenario=='high_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_COMPLETE","Miss_spec_to_add_Complete")],Branching_times=btimes_as_factor$col1)
    }
    if(Extinction_scenario=='low_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_COMPLETE","LOWIMPACT_Miss_spec_to_add_Complete")],Branching_times=btimes_as_factor$col1)
    }
  }
  
  if(Upham_tree_type=="DNA_only") { 
    if(Extinction_scenario=='high_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_DNAonly","Miss_spec_to_add_DNA_only")],Branching_times=btimes_as_factor$col1)
    }
    if(Extinction_scenario=='low_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_DNAonly","LOWIMPACT_Miss_spec_to_add_DNA_only")],Branching_times=btimes_as_factor$col1)
    }
  }
  colnames(the_mtable)<-c("Clade_name","Status","Missing_species", "Branching_times")
  
  ## I added this suppressWarnings because DAISIE_dataprep doesn't like the format of the branching times
  # with NAs somehow (even though NAs are correctly interpreted), but I didn't know how to fix it
  # In any case the output is correct
  Madagascar_list[[i]]<-suppressWarnings(DAISIE_dataprep(datatable = the_mtable, island_age=island_age,M=M))
  print(paste("Prepared type1 DAISIE data from tree",i,"out of", length(btimes_list)))
}

#### TYPE 2
Madagascar_list_type2<-list()
for (i in 1:length(btimes_list)){
  ##  order the list in case the table and btimes are not in the same order
  ordered_btimes_list<-btimes_list[[i]][clade_ordering]
  ## convert the branching times into a factor that can be added to a data frame
  btimes_as_factor =  data.frame(col1 = unlist(lapply(ordered_btimes_list, paste, collapse = ", ")))
  
  if(Upham_tree_type=="Complete") {   
    if(Extinction_scenario=='high_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_COMPLETE","Miss_spec_to_add_Complete")],Branching_times=btimes_as_factor$col1)
    }
    if(Extinction_scenario=='low_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_COMPLETE","LOWIMPACT_Miss_spec_to_add_Complete")],Branching_times=btimes_as_factor$col1)
    }
  }
  
  if(Upham_tree_type=="DNA_only") { 
    if(Extinction_scenario=='high_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_DNAonly","Miss_spec_to_add_DNA_only")],Branching_times=btimes_as_factor$col1)
    }
    if(Extinction_scenario=='low_impact'){
      the_mtable<-cbind(mtable[,c("Clade_name","DAISIE_status_DNAonly","LOWIMPACT_Miss_spec_to_add_DNA_only")],Branching_times=btimes_as_factor$col1)
    }
  }
  colnames(the_mtable)<-c("Clade_name","Status","Missing_species", "Branching_times")
  
  ## I added this suppressWarnings because DAISIE_dataprep doesn't accept the format of the branching times
  # with NAs for some reason (even though NAs are correctly interpreted), but I did not know how to fix it
  # In any case the output is correct
  Madagascar_list_type2[[i]]<-suppressWarnings(DAISIE_dataprep(datatable = the_mtable, island_age=island_age,M=M,
                                                               number_clade_types = 2, 
                                                               list_type2_clades = type2_bat_name_list,
                                                               prop_type2_pool = prop_type2_pool))
  print(paste("Prepared type1&2 DAISIE data from tree",i,"out of", length(btimes_list)))
  
}

if(posterior_or_maxcred=="maxcred"){
  save(Madagascar_list,file=paste("MAD_Maxcred_",Upham_tree_type,"_Mpool",M,"_IslandAge",island_age,"_Scenario",Colonisation_scenario,"_",Extinction_scenario,".Rdata",sep=""))
  save(Madagascar_list_type2,file=paste("MAD_Maxcred_",Upham_tree_type,"_Mpool",M,"_IslandAge",island_age,"_Scenario",Colonisation_scenario,"_",Extinction_scenario,"_ptype",prop_type2_pool,".Rdata",sep=""))}


if(posterior_or_maxcred=="posterior"){
  save(Madagascar_list,file=paste("MAD_1000_",Upham_tree_type,"_Mpool",M,"_IslandAge",island_age,"_Scenario",Colonisation_scenario,"_",Extinction_scenario,".Rdata",sep=""))
  save(Madagascar_list_type2,file=paste("MAD_1000_",Upham_tree_type,"_Mpool",M,"_IslandAge",island_age,"_Scenario",Colonisation_scenario,"_",Extinction_scenario,"_ptype",prop_type2_pool,".Rdata",sep=""))}
