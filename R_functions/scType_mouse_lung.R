###############################################################################
# Code from Claude on 14/10/2026, to convert the scType ref 
# ("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx")
# from human to mouse symbols, and add genes specific to lung tissues
#
# scType with a mouse marker database (lung)
#   1. convert ScTypeDB_full (human) -> mouse symbols (babelgene + manual map)
#   2. add mouse lung-specific cell types (tissueType = "Lung_mouse")
#   3. build gene sets WITHOUT gene_sets_prepare() so mouse case is preserved
###############################################################################
library(dplyr)
library(babelgene)
library(openxlsx)

# Download the full database
db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
temp_file <- tempfile(fileext = ".xlsx")
download.file(db_url, destfile = temp_file, mode = "wb")

# Read the database
db <- read.xlsx(temp_file)

# -----------------------------------------------------------------------------
# 1. Conversion
# -----------------------------------------------------------------------------

# Protein / surface-marker names, old aliases and human genes with a known mouse
# counterpart that babelgene cannot map automatically.
# (names = uppercase token as found in the db, values = mouse symbol(s))
# Everything NOT listed here and NOT mapped by babelgene is dropped. That
# includes pseudogenes, lncRNAs (RP11-*, CTD-*, LINC*, AC0*...), snoRNAs and
# human-specific genes (GNLY, S100A12, CD1A/CD1C, IL32, CCL14, FCN3, PSG*, CGB,
# TLR10, AZU1, PAEP, MYEOV, ...).
manual_map <- c(
  # immunoglobulins / TCR
  "IGD"="Ighd", "IGHD"="Ighd", "IGJ"="Jchain", "TRAC"="Trac", "TRDC"="Trdc",
  "CD3"="Cd3e,Cd3d,Cd3g", "CD3Z"="Cd247", "CD8"="Cd8a,Cd8b1",
  # CD nomenclature -> gene symbols
  "CD10"="Mme", "CD11B"="Itgam", "CD11C"="Itgax", "CD13"="Anpep",
  "CD15"="Fut4", "CD16"="Fcgr3,Fcgr4", "FCGR3"="Fcgr3,Fcgr4", "CD20"="Ms4a1",
  "CD25"="Il2ra", "CD29"="Itgb1", "CD30"="Tnfrsf8", "CD31"="Pecam1",
  "CD32"="Fcgr2b", "CD41"="Itga2b", "CD41A"="Itga2b", "CD42A"="Gp9",
  "CD42B"="Gp1ba", "CD43"="Spn", "CD45"="Ptprc", "CD45RA"="Ptprc",
  "CD45RO"="Ptprc", "CD56"="Ncam1", "CD61"="Itgb3", "CD62L"="Sell",
  "CD64"="Fcgr1", "CD73"="Nt5e", "CD90"="Thy1", "CD94"="Klrd1", "CD95"="Fas",
  "CD99"="Cd99", "CD105"="Eng", "CD110"="Mpl", "CD115"="Csf1r",
  "CD116"="Csf2ra", "CD117"="Kit", "CD122"="Il2rb", "CD123"="Il3ra",
  "CD125"="Il5ra", "CD127"="Il7r", "CD133"="Prom1", "CD135"="Flt3",
  "CD138"="Sdc1", "CD146"="Mcam", "CD150"="Slamf1", "CD161"="Klrb1c",
  "CD166"="Alcam", "CD193"="Ccr3", "CD203C"="Enpp3", "CD206"="Mrc1",
  "MRC"="Mrc1", "CD235A"="Gypa", "CD271"="Ngfr", "CD314"="Klrk1",
  "NKG2D"="Klrk1", "KLRC4-KLRK1"="Klrk1", "NKP46"="Ncr1", "CD49F"="Itga6",
  "SCA-1"="Ly6a", "FC-EPSILON RI-ALPHA"="Fcer1a",
  "NESTIN"="Nes", "OCT4"="Pou5f1", "KERATIN-19"="Krt19", "KERATIN-7"="Krt7",
  # MHC class II (functional equivalents)
  "HLA-DR"="H2-Aa,H2-Ab1,H2-Eb1", "HLA-DPB1"="H2-Ab1", "HLA-DPA1"="H2-Aa",
  # human genes with a different mouse symbol / functional analogue
  "CXCL8"="Cxcl1,Cxcl2",              # functional analogues, no Cxcl8 in mouse
  "C10ORF10"="Depp1", "SDPR"="Cavin2", "IFI15"="Isg15", "CLEC2B"="Clec2d",
  "CLEC141"="Clec14a", "G6PC"="G6pc1", "C6ORF25"="G6b", "FAM19A3"="Tafa3",
  "GPR97"="Adgrg3", "GPR116"="Adgrf5", "GPR123"="Adgra1", "BAI1"="Adgrb1",
  "WDR65"="Cfap57", "FAM26F"="Calhm6", "FAM132B"="Erfe", "AMICA1"="Jaml",
  "IFNG-AS1"="Tmevpg1", "WT1-AS"="Wt1os", "CXORF36"="Dipk2b",
  "MT1G"="Mt1", "PPAPDC1A"="Plpp4", "MTL5"="Tesmin", "C9ORF117"="Cfap157",
  "MFSD4"="Mfsd4a", "NOV"="Ccn3", "FYB"="Fyb1", "MIR145"="Mir145a",
  "NGN2"="Neurog2", "GCAP2"="Guca1b", "OPMW"="Opn1mw", "LIM1"="Lhx1",
  "SUR2B"="Abcc9", "GLYT1"="Slc6a9", "VACHT"="Slc18a3", "VGAT"="Slc32a1",
  "SLIT"="Slit1,Slit2,Slit3", "NR1I4"="Nr1i3", "NMT"="Nmt1", "UNC13"="Unc13a",
  # same gene in mouse, babelgene just has no 1:1 entry
  "ARSE"="Arse", "ELANE1"="Elane", "NGP"="Ngp", "CHIL3"="Chil3",
  "TCN1"="Tcn1", "GCSAML"="Gcsaml", "MEG3"="Meg3", "ABCA8A"="Abca8a",
  "OG9X"="Og9x", "FGF15"="Fgf15", "CAR8"="Car8", "KRT3"="Krt3",
  "SLCO1A4"="Slco1a4", "ZPF58"="Zfp58", "ATP5MD"="Atp5md", "TTC40"="Ttc40",
  "MIR144"="Mir144", "SPRR2F"="Sprr2f", "APCDD1L"="Apcdd1l",
  "MIR143HG"="Mir143hg", "EFCAB1"="Efcab1", "SOX2-OT"="Sox2ot",
  "CEMP1"="Cemp1", "GPR1"="Gpr1", "MIR218-1"="Mir218-1", "MIR194-2"="Mir194-2",
  "CD163L1"="Cd163l1", "FRMD1"="Frmd1",
  # symbols that were already mouse genes in the db (just restore the case)
  "GM10754"="Gm10754", "GM20687"="Gm20687", "GM13861"="Gm13861",
  "GM4876"="Gm4876", "IFI203"="Ifi203", "TRP73"="Trp73", "LHX1OS"="Lhx1os",
  "OTX2OS1"="Otx2os1", "LY6A"="Ly6a", "ABCB1A"="Abcb1a", "ZFP366"="Zfp366",
  "KITL"="Kitl", "AOX3"="Aox3", "ZFHX2OS"="Zfhx2os", "DLX6OS1"="Dlx6os1",
  "DLX1AS"="Dlx1as", "CYP2J12"="Cyp2j12", "UGT8A"="Ugt8a",
  "SLC6A20A"="Slc6a20a", "AIRN"="Airn", "DCDC2A"="Dcdc2a",
  "TBX3OS1"="Tbx3os1", "CLEC4A1"="Clec4a1", "MS4A6D"="Ms4a6d",
  "LYZ2"="Lyz2", "MS4A6C"="Ms4a6c"
)

split_genes <- function(x) {
  if (is.na(x)) return(character(0))
  g <- toupper(trimws(strsplit(x, ",")[[1]]))
  g[g != ""]
}

all_genes <- unique(unlist(lapply(c(db$geneSymbolmore1, db$geneSymbolmore2), split_genes)))
ort <- orthologs(genes = all_genes, species = "mouse", human = TRUE)
map <- split(ort$symbol, ort$human_symbol)

convert_str <- function(x) {
  g <- split_genes(x)
  if (!length(g)) return(NA_character_)
  out <- vapply(g, function(s) {
    if (s %in% names(manual_map)) manual_map[[s]]
    else if (s %in% names(map))   paste(map[[s]], collapse = ",")
    else ""
  }, character(1), USE.NAMES = FALSE)
  m <- unique(unlist(strsplit(out, ",")))
  m <- m[m != ""]
  if (!length(m)) NA_character_ else paste(m, collapse = ",")
}

db_mouse <- db
db_mouse$geneSymbolmore1 <- vapply(db$geneSymbolmore1, convert_str, character(1), USE.NAMES = FALSE)
db_mouse$geneSymbolmore2 <- vapply(db$geneSymbolmore2, convert_str, character(1), USE.NAMES = FALSE)

# what is still not converted (should now be only lncRNAs, pseudogenes, etc.)
still_unmapped <- setdiff(all_genes, c(names(map), names(manual_map)))
head(still_unmapped, 50)

write.xlsx(db_mouse, "ScTypeDB_mouse_wolung.xlsx")

# -----------------------------------------------------------------------------
# 2. OPTIONAL - Mouse lung-specific cell types from Claude
#    (negative markers are deliberately few: they subtract from the score)
# -----------------------------------------------------------------------------
mk <- function(cell, pos, neg = NA_character_, short = cell) {
  data.frame(tissueType = "Lung_mouse", cellName = cell,
             geneSymbolmore1 = pos, geneSymbolmore2 = neg,
             shortName = short, stringsAsFactors = FALSE)
}

lung_rows <- bind_rows(
  ## ---- Epithelium --------------------------------------------------------
  mk("Alveolar type 1 cells (AT1)", "Ager,Hopx,Pdpn,Aqp5,Rtkn2,Emp2,Clic5,Akap5,Cav1,Col4a3", "Sftpc,Lamp3", "Lung_AT1"),
  mk("Alveolar type 2 cells (AT2)", "Sftpc,Sftpb,Sftpa1,Sftpd,Lamp3,Abca3,Slc34a2,Napsa,Lyz1,Cxcl15,Etv5,Cebpa,Scd1,Lpcat1,Chil1,Muc1", "Ager,Hopx", "Lung_AT2"),
  mk("Transitional AT2 / DATP", "Krt8,Cldn4,Sfn,Lgals3,Cdkn1a,Krt19,Krt7", NA, "Lung_TransitionalAT2"),
  mk("Club cells", "Scgb1a1,Scgb3a2,Cyp2f2,Cbr2,Cyp4b1,Sec14l3,Scgb3a1", "Foxj1", "Lung_Club"),
  mk("Ciliated cells", "Foxj1,Tppp3,Dynlrb2,Ccdc153,Ccdc39,Cdhr3,Rsph1,Pifo,Sntn,Cfap126,Capsl,Mlf1,Dnah5,Tubb4b,Ccdc113", "Scgb1a1", "Lung_Ciliated"),
  mk("Goblet cells", "Muc5ac,Muc5b,Agr2,Spdef,Foxa3,Tff2,Bpifb1,Gp2,Reg3g,Cxcl17", "Foxj1", "Lung_Goblet"),
  mk("Basal cells", "Trp63,Krt5,Krt14,Krt15,Krt17,Ngfr,Aqp3,Itga6,Bcam", "Scgb1a1", "Lung_Basal"),
  mk("Pulmonary neuroendocrine cells", "Ascl1,Calca,Chga,Chgb,Syp,Uchl1,Scg5,Grp,Neurod1", NA, "Lung_PNEC"),
  mk("Pulmonary ionocytes", "Foxi1,Cftr,Ascl3,Atp6v1b1", NA, "Lung_Ionocyte"),
  mk("Tuft cells", "Pou2f3,Dclk1,Trpm5,Avil,Gnat3,Alox5,Il25,Ptgs1", NA, "Lung_Tuft"),

  ## ---- Endothelium -------------------------------------------------------
  mk("Endothelial cells (pan)", "Pecam1,Cdh5,Cldn5,Kdr,Tie1,Emcn,Eng,Flt1", NA, "Lung_EC"),
  mk("General capillary EC (gCap)", "Gpihbp1,Plvap,Kit,Cd93,Cd36,Ptprb", "Car4,Ednrb", "Lung_gCap"),
  mk("Aerocytes (aCap)", "Car4,Ednrb,Apln,Tbx2,Hpgd,Prx,Sostdc1,Emcn", "Gpihbp1,Plvap", "Lung_aCap"),
  mk("Arterial EC", "Gja5,Gja4,Bmx,Dkk2,Efnb2,Sox17,Fbln5,Sema3g,Hey1", NA, "Lung_ArterialEC"),
  mk("Venous EC", "Vwf,Nr2f2,Ephb4,Ackr1,Slc6a2,Selp,Sele,Vcam1,Cpe,Prss23", NA, "Lung_VenousEC"),
  mk("Lymphatic EC", "Prox1,Lyve1,Pdpn,Flt4,Ccl21a,Mmrn1,Reln,Tbx1", NA, "Lung_LymphaticEC"),

  ## ---- Mesenchyme --------------------------------------------------------
  mk("Alveolar fibroblasts", "Npnt,Ces1d,Scube2,Tcf21,Gpc3,Limch1", "Pi16", "Lung_Alv.fibroblast"),
  mk("Adventitial fibroblasts", "Pi16,Dpt,Sfrp2,Ly6a,Mfap5,Igfbp6,Serpinf1,Col1a1", "Npnt", "Lung_Adv.fibroblast"),
  mk("Peribronchial fibroblasts", "Hhip,Aspn,Fgf18,Lgr6", NA, "Lung_Peribronchial_fib."),
  mk("Lipofibroblasts", "Tcf21,Plin2,Lpl,Fgf10", NA, "Lung_Lipofibroblast"),
  mk("Myofibroblasts", "Acta2,Tagln,Pdgfra,Eln,Tgfbi", NA, "Lung_Myofibroblast"),
  mk("Fibrotic fibroblasts (Cthrc1+)", "Cthrc1,Postn,Tnc,Ddah1,Col1a1,Col3a1,Fn1,Acta2", NA, "Lung_Fibrotic_fib."),
  mk("Pericytes", "Pdgfrb,Cspg4,Rgs5,Kcnj8,Abcc9,Higd1b,Cox4i2,Notch3,Mcam", NA, "Lung_Pericyte"),
  mk("Vascular smooth muscle", "Acta2,Myh11,Cnn1,Tagln,Lmod1,Myocd", NA, "Lung_VascularSMC"),
  mk("Airway smooth muscle", "Acta2,Myh11,Actg2,Des,Cnn1,Tagln", NA, "Lung_AirwaySMC"),
  mk("Mesothelial cells", "Msln,Upk3b,Wt1,Lrrn4,Upk1b,Krt19,Rspo1", NA, "Lung_Mesothelial"),
  mk("Chondrocytes (trachea)", "Col2a1,Acan,Sox9,Col9a1,Col9a2,Col11a1", NA, "Lung_Chondrocyte"),

  ## ---- Myeloid -----------------------------------------------------------
  mk("Alveolar macrophages", "Siglecf,Ear1,Ear2,Chil3,Itgax,Marco,Pparg,Car4,Krt79,Fabp1,Fabp5,Mctp1,Atp6v0d2", "Cx3cr1,Ccr2", "Lung_AM"),
  mk("Interstitial macrophages", "C1qa,C1qb,C1qc,Cx3cr1,Lyve1,Folr2,Pf4,Cd68,Mertk,Fcgr1,Mrc1,Selenop", "Siglecf", "Lung_IM"),
  mk("Spp1+ profibrotic macrophages", "Spp1,Trem2,Gpnmb,Cd9,Lgals3,Ctsd,Fabp5", NA, "Lung_Spp1+_Mac"),
  mk("Classical monocytes", "Ly6c2,Ccr2,Plac8,Sell,Cd14", NA, "Lung_Ly6C-hi_Mono"),
  mk("Non-classical monocytes", "Nr4a1,Ace,Cx3cr1,Cd300e,Treml4,Fcgr4,Adgre4", "Ccr2,Sell", "Lung_Ly6C-lo_Mono"),
  mk("Neutrophils", "Ly6g,S100a8,S100a9,Retnlg,Cxcr2,Csf3r,Mmp9,Ngp,Camp,Lcn2,Ltf,Mmp8", NA, "Lung_Neutrophil"),
  mk("Eosinophils", "Siglecf,Ccr3,Il5ra,Epx,Prg2,Ear6,Gata1", "Itgax", "Lung_Eosinophil"),
  mk("Basophils", "Mcpt8,Cd200r3,Ms4a2,Fcer1a,Il4", NA, "Lung_Basophil"),
  mk("Mast cells", "Kit,Fcer1a,Cpa3,Mcpt4,Tpsb2,Tpsab1,Ms4a2,Hdc,Cma1", NA, "Lung_Mast"),
  mk("cDC1", "Xcr1,Clec9a,Batf3,Irf8,Cadm1,Itgae,Tlr3,Cd24a", NA, "Lung_cDC1"),
  mk("cDC2", "Sirpa,Cd209a,Irf4,Itgam,Clec4a4,Clec10a", NA, "Lung_cDC2"),
  mk("Migratory DC (mregDC)", "Ccr7,Fscn1,Ccl22,Lamp3,Cd40,Cd83,Il12b,Cd274,Pdcd1lg2,Cd200", "Sftpc", "Lung_mregDC"),
  mk("Plasmacytoid DC", "Siglech,Bst2,Ly6d,Tcf4,Irf7,Cox6a2,Ccr9,Klk1,Runx2,Pacsin1", NA, "Lung_pDC"),

  ## ---- Lymphoid ----------------------------------------------------------
  mk("B cells (mouse)", "Cd19,Ms4a1,Cd79a,Cd79b,Pax5,Cd22,Ighm,Ighd,Cr2,Fcer2a,Ebf1,Bank1,Cd72", "Cd3e", "B"),
  mk("Plasma cells (mouse)", "Jchain,Sdc1,Xbp1,Prdm1,Mzb1,Tnfrsf17,Derl3,Slamf7,Ssr4,Pou2af1,Ighg1,Igha", "Ms4a1,Cd19", "Plasma"),
  mk("T cells (pan)", "Cd3e,Cd3d,Cd3g,Cd247,Trac,Trbc1,Trbc2,Lck,Cd2,Cd5,Cd6", NA, "T"),
  mk("CD4 T cells", "Cd4,Il7r,Tcf7,Cd40lg,Icos,Cd28", "Cd8a,Cd8b1", "CD4 T"),
  mk("Regulatory T cells", "Foxp3,Il2ra,Ikzf2,Ctla4,Tnfrsf4,Tnfrsf18,Izumo1r,Cd4", "Cd8a", "Treg"),
  mk("Th1 cells", "Tbx21,Ifng,Cxcr3,Il12rb2,Stat4,Cd4", NA, "Th1"),
  mk("Th2 cells", "Gata3,Il4,Il5,Il13,Ccr4,Cd4", NA, "Th2"),
  mk("Th17 cells", "Rorc,Il17a,Il17f,Il23r,Ccr6,Il22,Il1r1,Cd4", NA, "Th17"),
  mk("T follicular helper cells", "Cxcr5,Pdcd1,Bcl6,Icos,Maf,Cd4", NA, "Tfh"),
  mk("Naive CD8 T cells (mouse)", "Cd8a,Cd8b1,Ccr7,Sell,Tcf7,Lef1,Cd3e", "Cd4", "Naive CD8"),
  mk("Effector CD8 T cells", "Cd8a,Cd8b1,Gzma,Gzmb,Prf1,Ifng,Nkg7,Klrg1,Ccl5,Cd3e", "Cd4", "Eff CD8"),
  mk("Tissue-resident memory CD8 T cells", "Cd8a,Cd8b1,Cd69,Itgae,Itga1,Cxcr6,Rgs1,Cd3e", "Cd4", "Trm CD8"),
  mk("Exhausted CD8 T cells", "Cd8a,Pdcd1,Havcr2,Lag3,Tigit,Tox,Entpd1,Cd3e", "Cd4", "Exh CD8"),
  mk("Gamma-delta T cells", "Trdc,Tcrg-C1,Tcrg-C2,Cd3e,Sox13,Blk,Cd163l1,Il17a,Rorc,Ccr6,Il1r1,Maf", "Cd4,Cd8a", "gd T"),
  mk("NK cells (mouse)", "Ncr1,Klrb1c,Klrk1,Nkg7,Gzma,Prf1,Klrd1,Eomes,Itga2,Klre1", "Cd3e", "NK"),
  mk("iNKT cells", "Zbtb16,Cd3e,Klrb1c,Cd44,Icos,Il4", NA, "iNKT"),
  mk("ILC2", "Il1rl1,Gata3,Il5,Il13,Klrg1,Il17rb,Areg,Arg1,Rora,Il7r,Icos", "Cd3e", "ILC2"),
  mk("ILC3", "Rorc,Il22,Il23r,Kit,Ccr6,Nrp1,Il7r", "Cd3e", "ILC3"),

  ## ---- Other -------------------------------------------------------------
  mk("Megakaryocytes / platelets", "Pf4,Ppbp,Itga2b,Gp9,Gp1bb,Mpl,Gp5", NA, "Mk/Platelet"),
  mk("Erythroid cells", "Hba-a1,Hba-a2,Hbb-bs,Hbb-bt,Alas2,Slc4a1,Gypa,Snca", NA, "Erythroid")
)

db_all <- bind_rows(db_mouse, lung_rows)
write.xlsx(db_all, "ScTypeDB_mouse_withlung.xlsx")

# -----------------------------------------------------------------------------
# 3. Gene sets WITHOUT gene_sets_prepare() from scType  (keeps mouse case, no HGNC check)
# -----------------------------------------------------------------------------
# prepare_gene_sets_mouse <- function(db, tissue) {
#   d <- db[db$tissueType %in% tissue, ]
#   parse <- function(x) lapply(x, function(s) {
#     if (is.na(s) || s == "") character(0) else unique(trimws(strsplit(s, ",")[[1]]))
#   })
#   gs  <- parse(d$geneSymbolmore1); names(gs)  <- d$cellName
#   gs2 <- parse(d$geneSymbolmore2); names(gs2) <- d$cellName
#   list(gs_positive = gs, gs_negative = gs2)
# }
# 
# gs_list <- prepare_gene_sets_mouse(db_all, unique(db_all$tissueType))
# gs_list$gs_positive[["Alveolar type 2 cells (AT2)"]]   # should be Sftpc, Sftpb, ...

# -----------------------------------------------------------------------------
# 4. Scoring (note gene_names_to_uppercase = FALSE)
# -----------------------------------------------------------------------------
# library(Seurat); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# es.max <- sctype_score(scRNAseqData = seu[["RNA"]]@scale.data, scaled = TRUE,
#                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative,
#                        gene_names_to_uppercase = FALSE)
#
# Tip: scale.data only contains the genes that were scaled (variable features by
# default). Run ScaleData(seu, features = rownames(seu)) first, otherwise many
# markers (e.g. Siglecf, Foxj1) may be missing from the scoring matrix.




###############################################################################
# Using the file "Cell_marker_Mouse.xlsx" downloaded 14/10/2026 from 
# http://www.bio-bigdata.center/CellMarker_download.html
# and converting it to a file compatible with scType
###############################################################################

# Load required libraries
library(dplyr)

# Read the input file
input_file <- "Cell_marker_Mouse.xlsx"
data <- read.xlsx(input_file)

# Clean the data: remove rows with missing gene symbols or tissue types
data_clean <- data %>%
  filter(!is.na(Symbol) & !is.na(tissue_type)) %>%
  mutate(Symbol = trimws(Symbol))  # Trim whitespace

# Group by tissue_type and aggregate gene symbols
scType_format <- data_clean %>%
  group_by(tissue_type) %>%
  summarise(
    tissueType = first(tissue_type),
    cellName = first(tissue_type), 
    geneSymbolmore1 = paste(Symbol, collapse = ","),
    geneSymbolmore2 = "",  # There is no such classification in the CellMarker file
    shortName = first(tissue_type)
  ) %>%
  ungroup() %>%
  select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2, shortName)

# Write the output file
output_file <- "cellmarker_sctype_format.xlsx"
write.xlsx(scType_format, output_file, col_names = TRUE)

# Print confirmation
head(scType_format)
dim(scType_format) # 271 cell types