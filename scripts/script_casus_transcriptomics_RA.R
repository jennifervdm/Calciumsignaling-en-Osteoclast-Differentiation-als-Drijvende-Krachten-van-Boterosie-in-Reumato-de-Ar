#----
# Mapping
# Deze data is gemaakt met de paired end sequencing methode (fw en rv)
setwd("C:/Users/Jenni/OneDrive - NHL Stenden/j2/p4/transcriptomics/casus reuma/Data_RA_raw/")
getwd()
# Package voor Bioconductor, voor bio data zoals RNA- en DNA-seq
install.packages('BiocManager')
library(BiocManager)
BiocManager::install('Rsubread')
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
setwd("C:/Users/Jenni/OneDrive - NHL Stenden/j2/p4/transcriptomics/casus reuma/Data_RA_raw/")
getwd()
# GTF gen-notaties downloaden
# Zorg dat GTF uit dezelfde database komt als GCF van eerder, anders kunnen er kleine variaties zijn
# Vergelijken BAM met GTF in een lijst (meerdere matrixen)
# Wat doet featureCounts en zijn zijn functies?
?featureCounts
# Nu gaan we hetzelfde doen met alle samples in 1 lijst
# Eerst een matrix maken voor alle samples
all.samples = c('SRR4785819.BAM',
                'SRR4785820.BAM',
                'SRR4785828.BAM',
                'SRR4785831.BAM',
                'SRR4785979.BAM',
                'SRR4785980.BAM',
                'SRR4785986.BAM',
                'SRR4785988.BAM')
# Dan de matrix invullen in featureCounts (pairedEnd TRUE/FALSE)
count_matrix <- featureCounts(
  files = all.samples,
  annot.ext = "genomic_human.gtf",
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
# Nu wordt de volledige versie gestuurd door een docent
#----
# Statistiek en analyse
# Inladen tekstbestand
setwd("C:/Users/Jenni/OneDrive - NHL Stenden/j2/p4/transcriptomics/casus reuma/")
getwd()
file.exists("count_matrix_RA.txt")
count_matrix_RA = read.delim("count_matrix_RA.txt")
head(count_matrix_RA)
str(count_matrix_RA)

# Download packages DESeq2, KEGGREST, EnhancedVolcano en pathview
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGGREST")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")

# Inladen packages
library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)

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
# Maak een tabel voor de dataframe
treatment_table_RA = data.frame(treatment_RA)
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
# Sla de resultaten op
write.table(resultaten_RA, file = "Resultaten_RA.csv", row.names = TRUE, col.names = TRUE)
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
# Sluit dev om de afbeelding in downloads te kunnen zien
dev.off()


# Nu een plot met interessante genen die niet te maken hebben met T-cellen
Interessante_Genen = subset(
  resultaten_RA)
# Filter de genen die met T-cellen te maken hebben
Tcel_Genen = c("CD3D", "CD3E", "CD4", "CD8A", "CD28", "CXCR1", "HLA-V", "RAB3IL1", "SRGN", "BCL2A1", "PTGFR", "ADAMDEC1")
# Maak de nieuwe dataset
Interessante_Genen = Interessante_Genen[
  !(rownames(Interessante_Genen) %in% Tcel_Genen), ]
# Maak de Volcano Plot met interessante genen
VolcanoPlot_RA_Interessante_Genen = EnhancedVolcano(Interessante_Genen,
                                 lab = rownames(Interessante_Genen),
                                 x = 'log2FoldChange',
                                 y = 'padj')

VolcanoPlot_RA_Interessante_Genen + scale_x_continuous(
  limits = c(-13, 13),
  breaks = seq(-14, 14, by = 2))
# Download de Volcano plot
dev.copy(png, 'Volcanoplot_RA_Interessante_Genen.png', 
         width = 10,
         height = 10,
         units = 'in',
         res = 500)
# Sluit dev om de afbeelding in downloads te kunnen zien
dev.off()






BiocManager::install("clusterProfiler")
library(clusterProfiler)

kegg = enrichKEGG(
  gene = entrez_ids,
  organism = "hsa"
)



# Voer nu de pathway analyse uit
# Hierbij worden de KEGG-pathways bekeken en vergeleken met de resultaten
# Zo wordt duidelijk welke pathways betrokken zijn bij bepaalde genen
# Op de website kegg.jp kunnen pathways gemapt worden met genen uit de Volcano plot
# Website: https://www.kegg.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=hsa


















# Voer nu de GO-analyse (gene ontology analyse) uit
# Hierbij kijk je naar de functies van genen 
# Je hebt een aantal packages nodig
#
library(BiocManager)
library(tidyverse)
library(dplyr)
BiocManager::install("goseq")
browseVignettes("goseq")
library(goseq)
BiocManager::install("geneLenDataBase")
library(geneLenDataBase)
BiocManager::install(org.Dm.eg.db)
library(org.Dm.eg.db)
# Overige packages
library(Rsubread)
library(Rsamtools)
library(DESeq2)
library(pathview)

# Je hebt een vector met alle genen en differentieel tot expressie gebrachte genen
# Hiervoor moet eerst je data gelezen worden
gene.vector <- as.integer(assayed.genes %in% de.genes)
names(gene.vector) <- assayed.genes
head(gene.vector)


#
supportedOrganisms() %>% filter(str_detect(Genome, "hg19"))


# Rijnamen opslaan, zodat je die kan gebruiken voor genen
rijnamen = rownames(resultaten_RA)
# Maak dataframe
dataframe.resultaten_RA = as.data.frame(resultaten_RA)
head(dataframe.resultaten_RA)


# Maak een lijst van differentieel tot expressie gebrachte genen
sigData <- as.integer(!is.na(shrinkLvV$FDR) & shrinkLvV$FDR < 0.05)
names(sigData) <- shrinkLvV$GeneID

pwf <- nullp(sigData, "hg19", "geneSymbol", bias.data = shrinkLvV$medianTxLength)

goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:BP"))

goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")




# Count matrix + resultaten + 
# Uiteindelijk een Go analyse, die moet zelf uitgezocht worden
# Zelf kiezen welke biologissche pathway? ahv Go analyse kiezen
# volcano, go analyse, pathview/deseq2(?)
# Kies biologische pathway, bespreek in discussie