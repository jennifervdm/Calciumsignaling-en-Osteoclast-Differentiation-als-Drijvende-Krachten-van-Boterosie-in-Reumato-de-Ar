#----
# Mapping
# Deze data is gemaakt met de paired end sequencing methode (fw en rv)
setwd("C:/Users/Jenni/OneDrive - NHL Stenden/j2/p4/transcriptomics/casus reuma/")
getwd()
# Package voor Bioconductor, voor bio data zoals RNA- en DNA-seq
install.packages('BiocManager')
BiocManager::install('Rsubread')
library(BiocManager)
library(Rsubread)

# Referentiegenoom downloaden, verplaats GCF file naar working directory
# Refgen kan via NCBI of Ensembl gedownload worden 
# Gebruik hetzelfde refgen voor de GTF gen-notaties

# Indexeren (download in je wd locatie)
# Dit is het verkleinen van de sequenties om groepen bp's te vergelijken (snel)
buildindex(basename = 'ref_human',
           reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
           memory = 4000,
           indexSplit = TRUE)

# Mappen data
# Downloaden alle samples
align.SRR4785819 <- align(index = "ref_human",
                         readfile1 = "SRR4785819_1_subset40k.fastq",
                         readfile2 = "SRR4785819_2_subset40k.fastq",
                         output_file = "SRR4785819.BAM")

align.SRR4785820 <- align(index = "ref_human",
                         readfile1 = "SRR4785820_1_subset40k.fastq",
                         readfile2 = "SRR4785820_2_subset40k.fastq",
                         output_file = "SRR4785820.BAM")

align.SRR4785828 <- align(index = "ref_human",
                         readfile1 = "SRR4785828_1_subset40k.fastq",
                         readfile2 = "SRR4785828_2_subset40k.fastq",
                         output_file = "SRR4785828.BAM")

align.SRR4785831 <- align(index = "ref_human",
                         readfile1 = "SRR4785831_1_subset40k.fastq",
                         readfile2 = "SRR4785831_2_subset40k.fastq",
                         output_file = "SRR4785831.BAM")

align.SRR4785979 <- align(index = "ref_human",
                          readfile1 = "SRR4785979_1_subset40k.fastq",
                          readfile2 = "SRR4785979_2_subset40k.fastq",
                          output_file = "SRR4785979.BAM")

align.SRR4785980_ <- align(index = "ref_human",
                          readfile1 = "SRR4785980_1_subset40k.fastq",
                          readfile2 = "SRR4785980_2_subset40k.fastq",
                          output_file = "SRR4785980.BAM")


align.SRR4785986 <- align(index = "ref_human",
                          readfile1 = "SRR4785986_1_subset40k.fastq",
                          readfile2 = "SRR4785986_2_subset40k.fastq",
                          output_file = "SRR4785986.BAM")

align.SRR4785988 <- align(index = "ref_human",
                          readfile1 = "SRR4785988_1_subset40k.fastq",
                          readfile2 = "SRR4785988_2_subset40k.fastq",
                          output_file = "SRR4785988.BAM")

# Laad Rsamtools voor sorteren en indexeren (dowloaden indien nodig)
# Dit is nodig om de data te kunnen zien in een viewer
BiocManager::install('Rsamtools')
library(Rsamtools)

# Bestandsnamen voor matrix van de monsters
samples <- c('SRR4785819',
             'SRR4785820',
             'SRR4785828',
             'SRR4785831',
             'SRR4785979',
             'SRR4785980',
             'SRR4785986',
             'SRR4785988')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})
# Indexeer de gesorteerde BAM-file
lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})

#----
# Count Matrix
library(Rsubread)
setwd("C:/Users/Jenni/OneDrive - NHL Stenden/j2/p4/transcriptomics/casus reuma")
getwd()

# GTF gen-notaties downloaden
# Zorg dat GTF uit dezelfde database komt als GCF van eerder, anders kunnen er kleine variaties zijn
# Vergelijken BAM met GTF in een lijst (meerdere matrixen)
# Wat doet featureCounts en zijn zijn functies?
??featureCounts
# Nu gaan we hetzelfde doen met alle samples in 1 lijst
# Eerst een matrix maken voor alle samples
all.samples = c('Data_RA_raw/BAM/SRR4785819.BAM',
                'Data_RA_raw/BAM/SRR4785820.BAM',
                'Data_RA_raw/BAM/SRR4785828.BAM',
                'Data_RA_raw/BAM/SRR4785831.BAM',
                'Data_RA_raw/BAM/SRR4785979.BAM',
                'Data_RA_raw/BAM/SRR4785980.BAM',
                'Data_RA_raw/BAM/SRR4785986.BAM',
                'Data_RA_raw/BAM/SRR4785988.BAM')

# Dan de matrix invullen in featureCounts (pairedEnd TRUE/FALSE)
count_matrix <- featureCounts(
  files = all.samples,
  annot.ext = "Data_RA_raw/genomic_human.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.featureType = "gene",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)
# Bekijk de matrix
str(count_matrix)
# Bekijk de annotation
head(count_matrix$annotation)
# We hebben alleen de counts nodig
counts = count_matrix$counts
# Check counts (kan ook in nieuw tabblad)
head(counts)
write.csv(counts, "human_countmatrix.csv")
# In dit geval wordt de volledige versie gestuurd door een docent

#----
# Statistiek en analyse
# Inladen tekstbestand
setwd("C:/Users/Jenni/OneDrive - NHL Stenden/j2/p4/transcriptomics/casus reuma/github/data_processed/")
getwd()
file.exists("count_matrix_RA.txt")
count_matrix_RA = read.delim("count_matrix_RA.txt")
head(count_matrix_RA)
str(count_matrix_RA)

# Download packages DESeq2 en EnhancedVolcano
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

# Inladen packages
library(DESeq2)
library(EnhancedVolcano)

# Metadata
# Maak een vector met categoriën van de behandelingen uit je dataset
treatment_RA= c("Control", 
                "Control", 
                "Control", 
                "Control",
                "Rheumatoid arthritis (established)",
                "Rheumatoid arthritis (established)",
                "Rheumatoid arthritis (established)", 
                "Rheumatoid arthritis (established)")
treatment_RA
# Maak een tabel voor de Metadataframe
treatment_table_RA = data.frame(treatment_RA)
View(treatment_table_RA)
# Download de Metadata
write.csv2(treatment_table_RA, "Metadata_RA.csv")

# Verander de kolomnamen van de .BAM files, zodat het geen .BAM meer heet
colnames(count_matrix_RA) = c("Control1", 
                              "Control2", 
                              "Control3", 
                              "Control4",
                              "RA1",
                              "RA2",
                              "RA3", 
                              "RA4")
head(count_matrix_RA)
#Verander de rijnamen van de behandeling tabel
rownames(treatment_table_RA) = c("Control1", 
                                 "Control2", 
                                 "Control3", 
                                 "Control4",
                                 "RA1",
                                 "RA2",
                                 "RA3", 
                                 "RA4")
head(treatment_table_RA)

# Statistiek analyses
# Maak de DESeqDataSet aan
dds_RA = DESeqDataSetFromMatrix(countData = count_matrix_RA,
                                colData = treatment_table_RA,
                                design = ~treatment_RA)
# Maak de DESeq dataset
dds_RA = DESeq(dds_RA)
# Voer de aanalyse uit voor de Volcano plot
resultaten_RA = results(dds_RA)
# Download de resultaten
write.table(resultaten_RA, file = "Differential_Expression_Resultaten_RA.cvs", row.names = TRUE, col.names = TRUE)
head(resultaten_RA)
# Bekijk hoeveel genen veranderd zijn
sum(resultaten_RA$padj < 0.05 & resultaten_RA$log2FoldChange > 1, na.rm = TRUE)
# ^Dit is 2085
sum(resultaten_RA$padj < 0.05 & resultaten_RA$log2FoldChange < -1, na.rm = TRUE)
# ^Dit is 2487
# Welke genen zijn opvallend?
hoogste_fold_change_RA = resultaten_RA[order(resultaten_RA$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change_RA = resultaten_RA[order(resultaten_RA$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde_RA = resultaten_RA[order(resultaten_RA$padj, decreasing = FALSE), ]
# Bekijk de volgende resultaten
head(hoogste_fold_change_RA)
head(laagste_fold_change_RA)
head(laagste_p_waarde_RA)

# Maak de Volcano plot
VolcanoPlot_RA = EnhancedVolcano(resultaten_RA,
                lab = rownames(resultaten_RA),
                x = 'log2FoldChange',
                y = 'padj')
# Verander de limieten van data en ruimte tussen labels
VolcanoPlot_RA + scale_x_continuous(
  limits = c(-13, 13),
  breaks = seq(-14, 14, by = 2))
# Download de Volcano plot
dev.copy(png, 'Volcanoplot_RA.png', 
         width = 10,
         height = 10,
         units = 'in',
         res = 500)
# Sluit dev om de afbeelding in de wd te kunnen zien
dev.off()


#---
# GO-analyse
# Voer nu de GO-analyse (enriched gene ontology terms) uit
# Hierbij kijk je naar de functies van genen 
# Tutorial: https://cloud.wikis.utexas.edu/wiki/spaces/bioiteam/pages/47732482/GO+Enrichment+using+goseq
source("http://bioconductor.org/biocLite.R")
install.packages("BiocManager")
BiocManager::install(c(
  "goseq",
  "org.Hs.eg.db",
  "AnnotationDbi"))
BiocManager::install(c(
  "DESeq2"))

library(goseq)
library(org.Hs.eg.db)
library(AnnotationDbi)
# Bekijk type ID codes, we hebben SYMBOL (ID: 5S-rRNA), Entrez (ID: 7157) is nodig
head(resultaten_RA)
# Bekijk de data of de gen namen op kolommen staan
colnames(resultaten_RA)
# Maak een vector voor de genen, nu wordt gekeken welke wel significant zijn (1 vs 0)
gene.vector <- as.integer(resultaten_RA$padj < 0.05)
names(gene.vector) <- rownames(resultaten_RA)

# Verwijder NA
keep <- !is.na(resultaten_RA$padj)
gene.vector <- gene.vector[keep]
names(gene.vector) <- rownames(resultaten_RA)[keep]

# Check keys
keys <- names(gene.vector)
head(keys)
length(keys)

# Verander ID codes van SYMBOL naar ENTREZ
mapping <- mapIds(
  org.Hs.eg.db,
  keys = keys,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first")

# Verander de ID codes nu naar een ENTREZ ID
names(gene.vector) <- mapping
gene.vector <- gene.vector[!is.na(names(gene.vector))]

# GOseq uitvoeren
pwf <- nullp(gene.vector, "hg38", "knownGene")
# Bekijk welke processen veel voorkomen in de lijst
GO.wall <- goseq(pwf, "hg38", "knownGene")
head(GO.wall)

GO_overrep_p05_vector <- GO.wall$category[GO.wall$over_represented_pvalue<.05]
GO_overrep_p05 <- GO.wall[GO.wall$category %in% GO_overrep_p05_vector, ]
View(GO_overrep_p05)

library(GO.db)
capture.output(for(go in GO_overrep_p05_vector[1:1517]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="GO-analyse_sig.txt")

# De data weergeeft category, over_represented_pvalue en nog een aantal sets
# Vooral de biologische processen zijn interessant, ook de p-waarde kan handig zijn
# In deze analyse vind je voornamelijk pathways die met het immuunsysteem te maken hebben
# Het meest interessante voor dit onderzoek zijn juist wat minder onderzochte processen


#----
# Pathway analyse
# Download clusterProfiler, pathview en KEGGREST
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGGREST")

# Inladen packages
library(clusterProfiler)
library(pathview)
library(KEGGREST)
library(org.Hs.eg.db)
library(AnnotationDbi)
# clusterProfiler kan evt. visualisatie geven aan GO-analyse https://doi.org/10.1089/omi.2011.0118

# KEGG pathway-analyse, dit kan per gen worden uitgevoerd
# Hierbij worden de KEGG-pathways bekeken en vergeleken met de resultaten
# Zo wordt duidelijk welke pathways betrokken zijn bij bepaalde genen
# Op de website kegg.jp kunnen pathways gemapt worden met genen uit de Volcano plot
# Website: https://www.kegg.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=hsa

# Maak de significante dataset van P-waarde <0.05 en log2FoldChange van > 1 & < -1
sig_genen <- rownames(
  resultaten_RA[
    !is.na(resultaten_RA$padj) &
      resultaten_RA$padj < 0.05 &
      abs(resultaten_RA$log2FoldChange) > 1,])
# Check data
length(sig_genen)
head(sig_genen)

# Convert SYMBOL naar Entrez
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = sig_genen,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first")
entrez_ids <- unique(na.omit(entrez_ids))

# Voer KEGG uit
kegg_resultaten_RA <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa")

# Maak een dataset van deze resultaten
as.data.frame(kegg_resultaten_RA)
head(as.data.frame(kegg_resultaten_RA))

# Visualisatie voor de meest voorkomende significante pathways
dotplot(kegg_resultaten_RA)
barplot(kegg_resultaten_RA)

# Maak de vector nodig voor de analyse
pathview_vector <- resultaten_RA$log2FoldChange
names(pathview_vector) <- row.names(resultaten_RA)
head(pathview_vector)

# Doe de analyse
pathview(
  gene.data = pathview_vector,
  pathway.id = "hsa04020",
  species = "hsa",
  gene.idtype = "SYMBOL",
  limit = list(gene = 5))