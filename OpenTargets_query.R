#install.packages("ghql")

library(dplyr)
library(ghql)
library(jsonlite)
require(biomaRt)

#convert variants from rsID to chromosome_position_ref_alt
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

snps <- unique(GOI_GWAS_100kb$SNP) #this is where I had the rsID from PGC3 stored


snp_biomart <- getBM(attributes=c(
  "refsnp_id", "chr_name", "chrom_start", "chrom_end", "chrom_strand",
  "allele",  "allele_1", "minor_allele", "associated_variant_risk_allele" ),
  filters="snp_filter", values=snps,
  mart=ensembl, uniqueRows=TRUE)

snp_biomart <- na_if(snp_biomart, '')
snp_biomart$allele_1 <- coalesce(snp_biomart$allele_1, snp_biomart$associated_variant_risk_allele)

snp_pos <- snp_biomart %>% select(chr_name, chrom_start, allele_1, minor_allele, refsnp_id)
snp_pos <- snp_pos %>%  filter(allele_1!=minor_allele) %>%  filter(chr_name %in% seq(1,22)) %>% distinct()
snp_pos_id <- str_c(snp_pos$chr_name, '_', snp_pos$chrom_start, '_', snp_pos$allele_1, '_', snp_pos$minor_allele)
snp_pos_id <- as.data.frame(cbind(snp_pos$refsnp_id, snp_pos_id))

#bipolar disorder PGC3 has identifier EFO_0000289 

identifier <- "EFO_0000289"


## Set up to query Open Targets Genetics API
otg_cli <- GraphqlClient$new(url = "https://api.genetics.opentargets.org/graphql")
otg_qry <- Query$new()

## Query for GWAS study locus details
otg_qry$query('l2g_query', 'query l2gQuery($studyId: String!, $variantId: String!){
    studyInfo(studyId: $studyId){
    numAssocLoci
    ancestryInitial
    nTotal
    nCases
    pubAuthor
  }
  studyLocus2GeneTable(studyId: $studyId, variantId: $variantId){
    rows {
      gene {
        id
        symbol
      }
      hasColoc
      yProbaModel
      distanceToLocus
    }
  }
}')

## Execute the query for one locus
variables <- list(studyId = "GCST012465", variantId = snp_pos_id$snp_pos_id[1])
result <- fromJSON(otg_cli$exec(otg_qry$queries$l2g_query, variables, flatten = TRUE))$data

#loop over all loci
l2g_loci <- list()
for( l in snp_pos_id$snp_pos_id){
  variables <- list(studyId = "GCST012465", variantId = l)
  result <- fromJSON(otg_cli$exec(otg_qry$queries$l2g_query, variables, flatten = TRUE))$data
  table <- result$studyLocus2GeneTable
  l2g_loci[l] <- table
}

l2g_loci <- l2g_loci[sapply(l2g_loci, is.data.frame)] #keep the ones where we got data (returned dataframes)


l2g_loci_df <- bind_rows(l2g_loci, .id= "SNP")

rs_not_scored <- snp_pos_id$V1[!snp_pos_id$snp_pos_id %in% l2g_loci_df$SNP] #find the rsID for the snps that was not successful in l2g scoring. We probably have the wrong formatted alleles for these
#look these up manually to get format

snp_not_scored <- c("1_163776152_C_T", "15_90883330_G_A", "15_42612706_A_C",  "11_66557112_T_C", "8_34294974_G_A", 
                    "6_98117335_C_T", "6_166581772_C_G", "14_99252882_A_G", "2_192873610_T_A", "22_40757875_T_C",
                    "9_37090541_T_TTAC", "16_9136959_G_A", "7_132185838_A_C")

snp_not_scored <- cbind(rs_not_scored, snp_not_scored)

#loop over remaining loci
l2g_loci_rest <- list()
for( l in snp_not_scored){
  variables <- list(studyId = "GCST012465", variantId = l)
  result <- fromJSON(otg_cli$exec(otg_qry$queries$l2g_query, variables, flatten = TRUE))$data
  table <- result$studyLocus2GeneTable
  l2g_loci_rest[l] <- table
}

l2g_loci_rest <- l2g_loci_rest[sapply(l2g_loci_rest, is.data.frame)] #keep the ones where we got data (returned dataframes). Missing intergenic loci

l2g_loci_rest_df <- bind_rows(l2g_loci, .id= "SNP")

#add to the previous round

l2g_loci_all <- bind_rows(l2g_loci_df , l2g_loci_rest_df)
l2g_loci_all <- l2g_loci_all %>% left_join(snp_pos_id, by = c("SNP"="snp_pos_id")) #insert rsID
l2g_loci_all <- l2g_loci_all%>% rename("rsID" = V1) %>%  distinct()


l2g_loci_all <- l2g_loci_all %>% left_join(select(GOI_GWAS_100kb, Locus, SNP, CHR, P, OR, Source, Nearest_gene), 
                                           by= c("rsID" ="SNP" )) #insert more information on locus, position etc

l2g_loci_all <- l2g_loci_all %>% rename("l2g"= yProbaModel)
l2g_loci_all$symbol <- l2g_loci_all$gene$symbol
l2g_loci_all$gene_id <- l2g_loci_all$gene$id

l2g_loci_all <- l2g_loci_all %>%  distinct()

l2g_loci_all <- l2g_loci_all %>% select(-gene,  -hasColoc)

#get the gene with highest l2g score from each locus
genes_max_l2g_locus <- l2g_loci_all %>%  
  group_by(Locus) %>% 
  mutate(max_l2g= max(l2g)) %>% 
  filter(max_l2g== l2g) %>%  
  distinct() %>%  
  arrange(desc(l2g)) %>% 
  select(Locus, symbol, l2g)

#get all genes with l2g scor >= 0.50
genes_l2g_0.5_locus <- 
  l2g_loci_all %>% 
  filter(l2g>= 0.5)  %>%  
  distinct() %>%  
  arrange(desc(l2g)) %>% 
  select(Locus, symbol, l2g)

dir <- "/Users/asbjorh/PhD/RNAseq_temporal/genes_of_interest/OpenTargets/"
write_tsv(genes_max_l2g_locus , file = paste0(dir, "/genes_max_l2g_locus.txt"))
write_tsv(genes_l2g_0.5_locus , file= paste0(dir, "/genes_l2g_0.5_locus.txt"))
write_tsv(l2g_loci_all, file= paste0(dir, "/genes_l2g_all_locus.txt"))


#look up the remaining manually

missing <- GOI_GWAS_100kb %>% 
  filter(!Locus %in% l2g_loci_all$Locus ) %>% 
  filter(Nearest_gene!= "intergenic") %>% 
  select(SNP, Locus, CHR, `Gene start (bp)`) %>%  distinct() 


#get snp position
snp_pos <- 
  getBM(attributes = c('chr_name', 'chrom_start', 'refsnp_id'),
        filters = 'snp_filter' ,
        values = GOI_GWAS_100kb$SNP,
        mart=ensembl )

snp_pos <- as_tibble(snp_pos) %>% filter(!str_detect(chr_name, c('PATCH', 'CHR'))) #remove the unwanted hits (non-canonical chr)
colnames(snp_pos) <- c("CHR", "pos" ,"SNP") 



#use the web interface to look up rsID, download Gene prioritisation using locus-to-gene pipeline from Mullins et al 2021.
f <- list.files(path = "/Users/asbjorh/PhD/RNAseq_temporal/genes_of_interest/OpenTargets/missing_in_query/",  full.names  = TRUE)

dat <- lapply(f, 
       function(i){          #get the data from each sample
         x = read_tsv(i, col_names = TRUE)
         x = x[, c(1,2,3,8)] 
         # add a column to say which file they're from
         x$file = i
         # Return your data
         x
       }
)


add <-
  do.call("rbind.data.frame", dat) #make a dataframe out of the imported datasets (combine by )

colnames(add) <- c("gene.symbol", "gene.id", "l2g", "distanceToLocus", "file") 
add$file <- #remove unnescessary prefix text from sample names
  str_remove(add$file, "/Users/asbjorh/PhD/RNAseq_temporal/genes_of_interest/OpenTargets/missing_in_query/") 

add$file <- #remove unnescessary suffix text from sample names
  str_remove(add$file, ".tsv") 

add <- add %>% separate( col = file, into = c(1,2,3,"SNP"), sep = "-")

add <- add %>%  select(SNP, gene.id, gene.symbol, l2g, distanceToLocus)

add <- add %>% 
  separate( col = SNP, into = c("CHR","pos","A1","A2"), sep = "_", remove = FALSE) 

add$pos <- as.integer(add$pos) 
add$CHR <- as.integer(add$CHR) 

snp_pos$CHR <- as.integer(snp_pos$CHR)
add <- add %>%
  left_join(snp_pos, by= c("CHR", "pos")) %>%  
  rename("SNP" = SNP.x, "rsID"= SNP.y)


#in need of rsID
last_rsID <- cbind(c("11_66087090_A_G", "15_42610048_A_G", "20_49432969_T_C"), 
                   c("rs489337", "rs1197546", "rs237475"))
last_rsID <- as.data.frame(last_rsID)
colnames(last_rsID) <- c("SNP", "rsID")


add <- add %>%
  left_join(last_rsID, by= c("SNP")) 

add <- add %>% mutate(rsID= coalesce(add$rsID.x, add$rsID.y))
add <- add %>% select(-rsID.x, -rsID.y )
add <- add %>% select(-A1, -A2 , -pos)

add <- add %>% rename("symbol" = gene.symbol)
add <- add %>% left_join(select(GOI_GWAS_100kb, Locus, SNP,  P, OR, Source, Nearest_gene), 
                        by= c("rsID"= "SNP")) %>% distinct()

#add the manually cured scores
l2g_loci_all$gene.id <- l2g_loci_all$gene$id
l2g_loci_all$gene.symbol <- l2g_loci_all$gene$symbol

l2g_loci_all_2 <- bind_rows(l2g_loci_all, add)
l2g_loci_all_2 <- l2g_loci_all_2 %>% select(SNP, symbol, gene.id, l2g, distanceToLocus, rsID, Locus, CHR, P, OR, Source, Nearest_gene)

l2g_loci_all <- l2g_loci_all_2

l2g_loci_all <- 
  l2g_loci_all %>% 
  mutate( rsID= str_replace(rsID, "rs489337", "rs475805")) %>% 
  mutate( SNP= str_replace(SNP, "11_66087090_A_G", "11_66081267_G_A")) 

l2g_loci_all <- 
  l2g_loci_all %>% 
  mutate( rsID= str_replace(rsID, "rs237475", "rs237460")) %>% 
  mutate(Locus = case_when(rsID == "rs237460" ~ 62))

l2g_loci_all <- 
  l2g_loci_all %>% 
  mutate( rsID= str_replace(rsID, "rs1197546", "rs4447398")) 
  
l2g_loci_all <- 
  l2g_loci_all %>% select(-Locus) %>%  left_join(select(GOI_GWAS_100kb, SNP, Locus), by= c("rsID" = "SNP")) 
  


#find loci missing high L2G
loci_missing <- setdiff(seq(1:64), genes_l2g_0.5_locus$Locus)
#get the nearest genes for these loci
loci_missing_nearest <- GOI_GWAS_100kb %>% 
  filter(Locus %in% loci_missing) %>% 
  select(Nearest_gene, Locus) %>% 
  filter(Nearest_gene != "intergenic") %>% 
  distinct()

write_tsv(loci_missing_nearest, file= paste0(dir, "loci_missing_nearest.txt"))

