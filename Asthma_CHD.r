library(TwoSampleMR)
library(readr)
library(MRPRESSO)
library(ieugwasr)
library(MendelianRandomization)
library(MRlap)
library(meta)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.NCBI.GRCh38)


##### Asthma meta-analysis instrumental variables #####

Asthma <- read.table(gzfile("Asthma_Bothsex_eas_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"), sep='\t', header = F)
Asthma <- Asthma[,-c(14:17)]
colnames(Asthma) <- c('CHR', 'POS', 'REF', 'ALT', 'rsid', 'all_meta_AF', 'inv_var_meta_beta', 'inv_var_meta_sebeta',
              'inv_var_meta_p', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
Asthma <- Asthma[Asthma$inv_var_meta_p<5e-8, ]

### clumping ###
dat <- data.frame(rsid=Asthma$rsid, pval=Asthma$inv_var_meta_p)
retained_snps <- ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, bfile = 'D:/LD_reference/EAS',
                          plink_bin = 'D:/LD_reference/plink.exe')
Asthma <- Asthma[Asthma$rsid%in%retained_snps$rsid, ]

### Exclude SNPs with minor allele frequency <0.01 ###
Asthma <- Asthma[(Asthma$all_meta_AF>0.01)&(Asthma$all_meta_AF<0.99), ]

write.csv(Asthma, 'Asthma_IVs_multiple_biobank.csv', quote = F, row.names = F)

Asthma <- read.csv('Asthma_IVs_multiple_biobank.csv', header = T, stringsAsFactors = F)
Asthma <- format_data(Asthma, chr_col = 'CHR', pos_col = 'POS', effect_allele_col = 'ALT',
                      other_allele_col = 'REF', snp_col = 'rsid', eaf_col = 'all_meta_AF',
                      beta_col = 'inv_var_meta_beta', se_col = 'inv_var_meta_sebeta', pval_col = 'inv_var_meta_p')
### estimate pseudo r2 of asthma
Asthma$samplesize.exposure <- 322655
add_rsq(Asthma)

### mean F-statistics ###

mean((Asthma$beta.exposure)^2/(Asthma$se.exposure)^2)

#### BBJ CAD Outcome ####
BBJ_CAD <- extract_outcome_data(snps=Asthma$SNP, outcomes = 'bbj-a-159', proxies = T)
Asthma.BBJ_CAD <- harmonise_data(Asthma, BBJ_CAD)
write.csv(Asthma.BBJ_CAD, 'Asthma_to_BBJ_CAD.csv', quote = F, row.names = F)

Asthma.BBJ_CAD <- read.csv('Asthma_to_BBJ_CAD.csv', header = T, stringsAsFactors = F)
ivw_res.BBJ_CAD <- mr(Asthma.BBJ_CAD, method_list = c('mr_ivw'))
ivw_res.BBJ_CAD$b <- 0.693*ivw_res.BBJ_CAD$b
ivw_res.BBJ_CAD$se <- 0.693*ivw_res.BBJ_CAD$se
generate_odds_ratios(ivw_res.BBJ_CAD)

#### UKB East Asian CAD Outcome ####

UKB_CAD <- read.table(gzfile("EA_cad.txt.gz"), sep=',', header = T)
UKB_CAD <- UKB_CAD[UKB_CAD$ID%in%Asthma$SNP, ]
UKB_CAD$BETA <- log(UKB_CAD$OR)
UKB_CAD <- format_data(UKB_CAD, type='outcome', snp_col = 'ID', chr_col = 'X.CHROM',
                       pos_col = 'POS', beta_col = 'BETA', se_col = 'SE', pval_col = 'P',
                       effect_allele_col = 'A1', other_allele_col = 'REF', eaf_col = 'ALT_FREQS')
Asthma.UKB_CAD <- harmonise_data(Asthma, UKB_CAD)
write.csv(Asthma.UKB_CAD, 'Asthma_to_UKB_CAD.csv', quote = F, row.names = F)

Asthma.UKB_CAD <- read.csv('Asthma_to_UKB_CAD.csv', header = T, stringsAsFactors = F)
ivw_res.UKB_CAD <- mr(Asthma.UKB_CAD, method_list = c('mr_ivw'))
ivw_res.UKB_CAD$b <- 0.693*ivw_res.UKB_CAD$b
ivw_res.UKB_CAD$se <- 0.693*ivw_res.UKB_CAD$se
generate_odds_ratios(ivw_res.UKB_CAD)


#### ivw-meta estimates ####
beta <- c(ivw_res.BBJ_CAD$b, ivw_res.UKB_CAD$b)
se <- c(ivw_res.BBJ_CAD$se, ivw_res.UKB_CAD$se)
dat <- data.frame(ids=c('BBJ', 'UKB'), beta=beta, se=se)
metagen(data=dat, TE=beta, seTE=se, overall=T, sm='OR', random = F)



##### MR-lap method #####


##### BBJ CAD data processing ####

BBJ_CAD <- read.table(gzfile("CAD.auto.rsq07.mac10.txt.gz"), sep=' ', header = T)
BBJ_CAD <- BBJ_CAD[,-c(3,6,13,14,15,16,17,18,19,20)]
write.table(BBJ_CAD, 'CAD_GWAS_BBJ.txt', quote=F, row.names = F)

BBJ_CAD <- read.table('CAD_GWAS_BBJ.txt', header=T, stringsAsFactors = F)
colnames(BBJ_CAD) <- c('CHR', 'POS', 'Other_allele', 'Effect_allele', 'EAF',
                   'N', 'Beta', 'SE', 'T-stat', 'P_value')
format_sumstats(path = BBJ_CAD, ref_genome = "GRCh37", dbSNP=144, save_path = "CAD.BBJ.tsv.gz")


##### UKB CAD data processing ####
UKB_CAD <- read.table(gzfile("EA_cad.txt.gz"), sep=',', header = T)
UKB_CAD <- UKB_CAD[,-c(5,7)]
UKB_CAD$BETA <- log(UKB_CAD$OR)
UKB_CAD <- UKB_CAD[,-which(colnames(UKB_CAD)=='OR')]
colnames(UKB_CAD) <- c('SNP', 'CHR', 'POS', 'Other_allele', 'Effect_allele', 'N',
                       'SE', 'T-stat', 'P_value', 'EAF', 'Beta')
UKB_CAD <- UKB_CAD[,-which(colnames(UKB_CAD)=='T-stat')]
UKB_CAD$N <- 89+2210
format_sumstats(path = UKB_CAD, ref_genome = "GRCh37", dbSNP=144, save_path = "CAD.UKB.tsv.gz")

##### Asthma data processing ####
Asthma <- read.table(gzfile("Asthma_Bothsex_eas_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"), sep='\t', header = F)
Asthma <- Asthma[,-c(14:17)]
colnames(Asthma) <- c('CHR', 'POS', 'Other_allele', 'effect_allele', 'SNP', 'EAF', 'Beta', 'SE',
                      'P_value', 'inv_var_het_p', 'direction', 'N_case', 'N_ctrl')
Asthma <- Asthma[(Asthma$CHR)!=23, ]
Asthma$N <- Asthma$N_case + Asthma$N_ctrl
Asthma <- Asthma[,-c(10,11,12,13)]
format_sumstats(path = Asthma, ref_genome = "GRCh38", dbSNP=144, save_path = "Asthma.EA.tsv.gz")

#### because MRlap function exclude the SNPs located within the MHC region
#### so we modified the MRlap function to retain the MHC region
#### the new function called MRlap_EAS.r

source('MRlap_EAS.r')

### BBJ CAD results ###
Asthma <- read.table('Asthma.EA.tsv.gz', header = T, stringsAsFactors = F)
CAD <- read.table('CAD.BBJ.tsv.gz', header = T, stringsAsFactors = F)
CAD <- CAD[,-c(10)]
colnames(CAD) <- c('snp', 'chr', 'pos', 'a1', 'a2', 'freq', 'N', 'beta', 'se', 'P')

Asthma <- Asthma[Asthma$SNP%in%CAD$snp, ]
CAD <- CAD[CAD$snp%in%Asthma$SNP, ]

pos.hg19 <- data.frame(SNP=CAD$snp, Chr_hg19=CAD$chr, pos_hg19=CAD$pos)
Asthma <- merge(Asthma, pos.hg19, by='SNP')
Asthma <- Asthma[,-c(2,3)]
colnames(Asthma) <- c('snp',  'a1', 'a2', 'freq', 'beta', 'se', 'P', 'N', 'chr', 'pos')

res1 <- MRlap(exposure = Asthma,
              exposure_name = 'Asthma',
              outcome = CAD,
              outcome_name = 'CAD',
              ld = 'D:/LD_reference/eas_w_ld_chr',
              hm3 = 'D:/LD_reference/w_hm3.noMHC.snplist',
              MR_pruning_LD = 0.001,
              MR_pruning_dist = 10000,
              MR_reverse = 0.001)
res1$MRcorrection$corrected_effect <- 0.693*res1$MRcorrection$corrected_effect
res1$MRcorrection$corrected_effect_se <- 0.693*res1$MRcorrection$corrected_effect_se
exp(res1$MRcorrection$corrected_effect)
lo_ci <- res1$MRcorrection$corrected_effect - qnorm(0.975)*res1$MRcorrection$corrected_effect_se
up_ci <- res1$MRcorrection$corrected_effect + qnorm(0.975)*res1$MRcorrection$corrected_effect_se
exp(c(lo_ci, up_ci))

### UKB CAD results ###
Asthma <- read.table('Asthma.EA.tsv.gz', header = T, stringsAsFactors = F)
CAD <- read.table('CAD.UKB.tsv.gz', header = T, stringsAsFactors = F)
colnames(CAD) <- c('snp', 'chr', 'pos', 'a1', 'a2', 'N', 'se', 'P', 'freq', 'beta')

Asthma <- Asthma[Asthma$SNP%in%CAD$snp, ]
CAD <- CAD[CAD$snp%in%Asthma$SNP, ]

pos.hg19 <- data.frame(SNP=CAD$snp, Chr_hg19=CAD$chr, pos_hg19=CAD$pos)
Asthma <- merge(Asthma, pos.hg19, by='SNP')
Asthma <- Asthma[,-c(2,3)]
colnames(Asthma) <- c('snp',  'a1', 'a2', 'freq', 'beta', 'se', 'P', 'N', 'chr', 'pos')

res2 <- MRlap(exposure = Asthma,
              exposure_name = 'Asthma',
              outcome = CAD,
              outcome_name = 'CAD',
              MR_threshold = 5e-8,
              ld = 'D:/LD_reference/eas_w_ld_chr',
              hm3 = 'D:/LD_reference/w_hm3.noMHC.snplist',
              MR_pruning_LD = 0.001,
              MR_pruning_dist = 10000,
              MR_reverse = 0.01)
res2$MRcorrection$corrected_effect <- 0.693*res2$MRcorrection$corrected_effect
res2$MRcorrection$corrected_effect_se <- 0.693*res2$MRcorrection$corrected_effect_se
exp(res2$MRcorrection$corrected_effect)
lo_ci <- res2$MRcorrection$corrected_effect - qnorm(0.975)*res2$MRcorrection$corrected_effect_se
up_ci <- res2$MRcorrection$corrected_effect + qnorm(0.975)*res2$MRcorrection$corrected_effect_se
exp(c(lo_ci, up_ci))

#### ivw-meta estimates ####
beta <- c(res1$MRcorrection$corrected_effect, res2$MRcorrection$corrected_effect)
se <- c(res1$MRcorrection$corrected_effect_se, res2$MRcorrection$corrected_effect_se)
dat <- data.frame(ids=c('BBJ', 'UKB'), beta=beta, se=se)
metagen(data=dat, TE=beta, seTE=se, overall=T, sm='OR', random = F)
