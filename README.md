# Calciumsignaling en Osteoclast Differentiation als Drijvende Krachten van Boterosie in Reumatoïde Artritis: een transcriptomics analyse
## 📁 Inhoud/structuur
- `assets` - Overige documenten voor de opmaak van deze pagina.
- `bronnen` - Gebruikte bronnen en AI gebruik.
- `data_processed` – Bewerkte data voor helder overzicht van resultaten en uitvoeren van analyses. 
- `data_raw` – Ruwe data afkomstig van 8 individuen, waarvan 4 RA hebben, en 4 gezond zijn.  
- `data_stewardship` - Aantoning van de competentie Beheren niveau I.
- `scripts` – Scripts waarin de data geanalyseerd wordt.
- `resultaten` - Figuren zoals grafieken en tabellen.
- `README.md` - Het document om de tekst hier te genereren.

---

## Inleiding

Reumatoïde artritis (RA) is een veelvoorkomende autoimmuunziekte. In 2024 waren er zo’n 200.000 Nederlanders met deze aandoening [(Reumatoïde artritis (RA) | Leeftijd en geslacht | Volksgezondheid en Zorg, z.d.)](<https://www.vzinfo.nl/reumatoide-artritis-ra/leeftijd-en-geslacht>). RA wordt gekenmerkt door synoviale hyperplasie met pannusvorming [(Chetina & Markova, 2019)](<https://doi.org/10.1134/S1990750819010049>). Synovium, een type slijmvlies raakt ontstoken, wat lijdt tot de erosie van bot- en kraakbeenweefsel wegens overname van pannusweefsel [(Gravallese & Monach, 2015)](<https://doi.org/10.1016/B978-0-323-09138-1.00094-2>). Momenteel focust onderzoek zich op de genen die met T-cellen te maken hebben. [(Padyukov, 2022)](<https://doi.org/10.1007/s00281-022-00912-0>). Er zijn al kleine klinische studies bezig met een mogelijk geneesmiddel, CAR T-celtherapie, wat veelbelovende behandelingsresultaten bied voor diverse autoimmuunziekten, waaronder RA [(Freeley, 2025)](<https://doi.org/10.1007/s12016-025-09113-7>). 

Transcriptomics is een veelgebruikte studie binnen genetisch en medisch onderzoek. Het transcriptoom, wat hierbij onderzocht wordt, geeft informatie over hoeveelheid expressie in alle genen. Met analyses van de functies van genen met hoge expressie, kan een beeld worden geschetst van de basale processen die plaatsvinden bij zieke individuen, in vergelijking met gezonde personen [(Monzó et al., 2025)](<https://doi.org/10.1038/s41576-025-00828-z>).

In dit onderzoek is ingedoken op de processen betrokken bij RA, door transcriptiefactoren te onderzoeken aan de hand van een VolcanoPlot, GO-analyse en KEGG-pathway analyse. De nadruk is gelegd op genen en gene ontologies waar minder onderzoek naar gedaan is om een dieper begrip te krijgen van diverse processen en pathways die betrokken zijn bij RA, in dit geval is focus gelegd op boterosie, één van de destructievere symptomen. Een interessante pathway hiervoor is osteoclast differentiation, wat botvernieuwing kan beïnvloeden [(Filgueira, 2010)](<https://www.sciencedirect.com/science/article/pii/B9780123748959000050>).

## Methoden

Er is ingezoomd op de transcriptomics van RA, om genen op te sporen met onderzoekspotentie [(figuur 1)](assets/Flowchart_Methode_RA.png). De dataset is afkomstig van het onderzoek van [Platzer et al. (2019)](<https://doi.org/10.1371/journal.pone.0219698>) [(tabel 1)](data_processed/Metadata_RA.cvs), en is uitgewerkt in R V4.5.2 in dit [script](scripts/script_casus_transcriptomics_RA.R).

*Tabel 1. Patiënten uit onderzochte dataset, afkomstig van het onderzoek van [Platzer et al. (2019)](<https://doi.org/10.1371/journal.pone.0219698>). De ruwe sequencing data in FASTQ bestanden is afkomstig van 8 vrouwen, 4 met RA (leeftijden 54-66), vastgesteld voor >12 maanden en positief getest op autoantistoffen ACPA. En een negatief geteste controlegroep van 4 (leeftijden 15-42).*
|     ID      | Age |   Sex   |               Status                |
| ----------- | --- | ------- | ----------------------------------- |
| SRR4785819  | 31  | female  |               Normal                |
| SRR4785820  | 15  | female  |               Normal                |
| SRR4785828  | 31  | female  |               Normal                |
| SRR4785831  | 42  | female  |               Normal                |
| SRR4785979  | 54  | female  | Rheumatoid arthritis (established)  |
| SRR4785980  | 55  | female  | Rheumatoid arthritis (established)  |
| SRR4785986  | 60  | female  | Rheumatoid arthritis (established)  |
| SRR4785988  | 59  | female  | Rheumatoid arthritis (established)  |

Bij data mapping wordt een index gebouwd met [`BiocManager`](<https://bioconductor.org/install/>) V1.30.27/V3.22 package [`Rsubread`](<https://bioconductor.org/packages/Rsubread/>) V2.24.0. Hierin worden FASTA bestanden van het RefSeq referentiegenoom [GRCh38.p14](<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/>) afkomstig van de NCBI genome database gemaakt. Het alignen tot BAM bestanden is in paired-end. [`Rsamtools`](<https://bioconductor.org/packages//release/bioc/html/Rsamtools.html>) V2.26.0 sorteert en indexeert de BAM files. De Count Matrix en BAM files matrix wordt gemaakt aan de hand van een GTF annotatiebestand, nogmaals het NCBI GRCh38.p14 genoom.

Een VolcanoPlot wordt gemaakt voor genexpressie bepaling. Hiervoor zijn de [`BiocManager`](<https://bioconductor.org/install/>) packages [`DESeq2`](<https://doi.org/10.3791/62528>) V1.50.2 (Liu et al., 2021), en [`EnhancedVolcano`](<http://bioconductor.org/packages/EnhancedVolcano/>) V1.28.2 nodig. Een differentiële analyse wordt over de data uitgevoerd. De VolcanoPlot bevat de log₂ fold change en gecorrigeerde P-waarde <0.05.

Bepalen van verschillende Gene Ontologies gaat via [`goseq`](<http://bioconductor.org/packages/goseq/>) V1.62.0, [`org.Hs.eg.db`](<http://bioconductor.org/packages/org.Hs.eg.db/>), en [`AnnotationDbi`](<http://bioconductor.org/packages/AnnotationDbi/>). Hier wordt een gecorrigeerde P-waarde van <0.05 aangehouden en het hg38 genoom wordt gebruikt. De ontologie [osteoclast differentiation](<https://amigo.geneontology.org/amigo/term/GO:0030316>) is onderzocht.

KEGG pathway analyses vereisen [`clusterProfiler`](<https://doi.org/10.1089/omi.2011.0118>) V4.18.4 (Yu et al., 2012), [`pathview`](<https://bioconductor.org/packages/release/bioc/html/pathview.html>) V1.50.0, [`KEGGREST`](<http://bioconductor.org/packages/KEGGREST/>) V1.50.0, [`org.Hs.eg.db`](<http://bioconductor.org/packages/org.Hs.eg.db/>) V3.22.0, en [`AnnotationDbi`](<http://bioconductor.org/packages/AnnotationDbi/>) V1.72.0. Enriched KEGG-analyses weergeven genen die meer of minder tot expressie komen in bepaalde pathways, en de pathview weergeeft een pathway-figuur, uitgevoerd op het gen van interesse [(CAMK4)](<https://www.kegg.jp/entry/hsa:CAMK4>) met een log2FoldChange vector.

<p align="center">
  <img src="assets/Flowchart_Methode_RA.png "alt="Flowchart Methode" width="600"/>
</p>

*Figuur 1. Flowchart van de gebruikte methode. Een referentiegenoom wordt geïndexeerd. RNA seq data en index worden gemapt naar BAM files. De reads van deze files worden vergeleken met de GTF annotatie en geteld. In de data-analyse worden een Volcano plot gemaakt, een GO-analyse, en KEGG Pathway-analyse uitgevoerd.*

## 📊 Resultaten

In de transcriptoomanalyse was de differentiële genexpressie tussen de condities van RA (n=4) en controle (n=4) onderzocht. Een Volcano plot is gemaakt om alle significante genen met differentiële expressie te weergeven, een GO-analyse voor significante biologische processen, en KEGG pathway analyse voor betrokken pathways bij handmatig geselecteerde genen.

### Volcano Plot

De Volcano plot weergaf 2085 opgereguleerde en 2487 neergereguleerde [(figuur 2)](resultaten/Volcanoplot_RA.png). Van de 29407 genen waren 4572 significant, P-waarde <0.05 en een log₂ fold change van > 1 en < -1 waren meegenomen. De meeste genen met expressie waren gerelateerd aan B-cellen, T-cellen en macrofagen. Er is een gen gekozen in verband met boterosie. Het gen [CAMK4](<https://www.kegg.jp/entry/hsa:CAMK4>) (Ca2+/calmodulin-dependent protein kinase) is onderzocht wegens betrekking tot osteoclast differentiation. De log₂ fold change van dit gen was opgereguleerd tot 3.31 met een gecorrigeerde P-waarde van 3.47e-05 [(tabel 2)](resultaten/CAMK4_differentiële-analyse.cvs).

<p align="center">
  <img src="resultaten/Volcanoplot_RA.png "alt="Resultaten Volcanoplot" width="600"/>
</p>

*Figuur 2. Volcano plot van genen met differentiële genexpressie bij RA. De data vergelijkt gezonde individuen (n=4) met RA-patienten (n=4). Genen n=29407, de P-waarde van <0.05 is meegenomen. De x-range loopt van -14 naar 14, aangezien alle data binnen deze punten ligt. Er is gekozen elke 2 waarden op de x-as aan te geven voor overzicht.*

*Tabel 2. Resultaten van de differentiële analyse van het gen CAMK4.*
|       | baseMean | log2FoldChange |   lfcSE   |  stat   |   pvalue    |     padj    |
| ----- | -------- | -------------- | --------- | ------- | ----------- | ----------- |
| CAMK4 | 236.8115 |    3.30909     |  0.689901 | 4.79647 | 1.61485e-06 | 3.46967e-05 |

### Gene Ontology Analyse

De Gene Ontology (GO)-analyse weergaf verschillende processen van RA-patiënten. De analyse toonde een significante verrijking binnen het proces [osteoclast differentiation](https://amigo.geneontology.org/amigo/term/GO:0030316) (0.005981202; overrep. p=0.005981202, underrep. p=0.9967119). 42 van de 105 pathways-gerelateerde genen waren differentieel tot expressie gekomen. De meeste genen waren niet significant binnen de Volcano Plot, en hierdoor minder relevant. CAMK4 was wel interessant wegens significante overrepresentatie.

### KEGG Pathway Analyse

[CAMK4](<https://www.kegg.jp/entry/hsa:CAMK4>) is geanalyseerd met een KEGG-pathway analyse. Het is betrokken bij calcium signaling, cAMP signaling en osteoclast differentiation. Aan de hand hiervan is een calcium signaling pathway analyse uitgevoerd [(figuur 3)](resultaten/hsa04020_pathview_RA.png). Veel betrokken genen bij de pathway waren niet significant binnen de Volcano plot, RYR en de CAMK familie wel. Binnen de pathway had CAMK een onderrepresentatie.

<p align="center">
  <img src="resultaten/hsa04020_pathview_RA.png" alt="Resultaten Pathway Analyse" width="600"/>
</p>

*Figuur 3. KEGG-pathway analyse van de calcium signaling pathway. Opgereguleerde genen zijn rood en neergereguleerde genen zijn groen.*

## Conclusie

De focus in het RA onderzoeksveld ligt momenteel op CAR T-cel therapie testen voor behandelingen  [(Freeley, 2025)](<https://doi.org/10.1007/s12016-025-09113-7>). Tot een grootschalig onderzoek uitgevoerd is, blijft het verstandig te investeren in alternatief onderzoek. De focus van dit onderzoek ligt op genen gerelateerd aan het destructieve RA symptoom boterosie. Toenemende kennis over boterosie kan een keerpunt vormen voor begrip van achterliggende processen van RA.

Met de Volcano Plot, KEGG-analyse en pathway analyse kan een relatie aangetoond worden met het functioneren van cellulaire processen van botten, zoals osteoclast differentiation. Binnen deze drie analyses is het gen CAMK4 significant gevonden, samen met een rol in de calcium signalering pathway en osteoclast differentiation lijkt het gen relevant te zijn. Het CAMK4 gen bewijst een grote rol te spelen in autoimmuunziekten, zoals SLE, EAE en EAP, door inbalans van T-cellen. Dit gen heeft ook invloed op parenchymale celmassa modulatie, zoals celproliferatie, apoptose en celbeschadiging [(Xu et al., 2024)](<https://www.sciencedirect.com/science/article/pii/S0006295224001874>). Ook speelt het een rol in tumorvorming, wat een link heeft met het CAR T-cel therapie, oorspronkelijk een kankermedicijn [(Fu et al., 2026)](<https://doi.org/10.1016/j.cellsig.2026.112557>). Dit suggereert dat het gen, indien disfunctioneel, invloed heeft op de werkingen die in beide tumoren en boterosie aanrichten. Hiermee heeft het gen mogelijk interactie met CAR T-cel therapie.

Uiteraard zal in verder vervolgonderzoek verdiept moeten worden in het CAMK4 gen om aantoning te kunnen doen van werking in het complete proces, er moet verder onderzoek gedaan worden naar de interactie van het gen en het osteoclast differentiation proces, en receptie naar CAR T-cel therapie. Ook vinden en uitspitten van andere significante genen zoals SNX10 die gerelateerd zijn aan osteoclast differentiatie wordt geadviseerd.

## Bronnen

AmiGo2 Osteoclast Differentiation. (z.d.). Geraadpleegd 11 juli 2026, van <https://amigo.geneontology.org/amigo/term/GO:0030316>

AnnotationDbi. (z.d.). Bioconductor. Geraadpleegd 29 mei 2026, van <http://bioconductor.org/packages/AnnotationDbi/>

Bioconductor—Install. (z.d.). Geraadpleegd 29 mei 2026, van <https://bioconductor.org/install/>

Chetina, E. V., & Markova, G. A. (2019). Prospects for the Use of Gene Expression Analysis in Rheumatology. Biochemistry (Moscow), Supplement Series B: Biomedical Chemistry, 13(1), 13-25. <https://doi.org/10.1134/S1990750819010049>

EnhancedVolcano. (z.d.). Bioconductor. Geraadpleegd 29 mei 2026, van <http://bioconductor.org/packages/EnhancedVolcano/>

Filgueira, L. (2010). Chapter 5—Osteoclast Differentiation and Function. In D. Heymann (Red.), Bone Cancer (pp. 59-66). Academic Press. <https://doi.org/10.1016/B978-0-12-374895-9.00005-0>

Freeley, M. (2025). CAR T Cell Therapy for Rheumatoid Arthritis. Clinical Reviews in Allergy & Immunology, 68(1), 100. <https://doi.org/10.1007/s12016-025-09113-7>

Fu, R., Wang, Y., Wang, W., Song, Q., Zhang, H., & Yang, L. (2026). The CAMK protein family: Pivotal roles in calcium-dependent disease pathways and emerging therapeutic strategies. Cellular Signalling, 144, 112557. <https://doi.org/10.1016/j.cellsig.2026.112557>

Goseq. (z.d.). Bioconductor. Geraadpleegd 29 mei 2026, van <http://bioconductor.org/packages/goseq/>

Gravallese, E. M., & Monach, P. A. (2015). The rheumatoid joint: Synovitis and tissue destruction. In Rheumatology (pp. 768-784). Mosby. <https://doi.org/10.1016/B978-0-323-09138-1.00094-2>

Homo sapiens genome assembly GRCh38.p14. (z.d.). NCBI. Geraadpleegd 28 mei 2026, van <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/>

KEGGREST. (z.d.). Bioconductor. Geraadpleegd 29 mei 2026, van <http://bioconductor.org/packages/KEGGREST/>

KEGG T01001: CAMK4. (z.d.). Geraadpleegd 13 juli 2026, van <https://www.kegg.jp/entry/hsa:CAMK4>

Liu, S., Wang, Z., Zhu, R., Wang, F., Cheng, Y., & Liu, Y. (2021). Three Differential Expression Analysis Methods for RNA Sequencing: Limma, EdgeR, DESeq2. Journal of Visualized Experiments: JoVE, (175). <https://doi.org/10.3791/62528>

Monzó, C., Liu, T., & Conesa, A. (2025). Transcriptomics in the era of long-read sequencing. Nature Reviews. Genetics, 26(10), 681-701. <https://doi.org/10.1038/s41576-025-00828-z>

Org.Hs.eg.db. (z.d.). Bioconductor. Geraadpleegd 29 mei 2026, van <http://bioconductor.org/packages/org.Hs.eg.db/>

Padyukov, L. (2022). Genetics of rheumatoid arthritis. Seminars in Immunopathology, 44(1), 47-62. <https://doi.org/10.1007/s00281-022-00912-0>

Pathview. (z.d.). Bioconductor. Geraadpleegd 29 mei 2026, van <http://bioconductor.org/packages/pathview/>

Platzer, A., Nussbaumer, T., Karonitsch, T., Smolen, J. S., & Aletaha, D. (2019). Analysis of gene expression in rheumatoid arthritis and related conditions offers insights into sex-bias, gene biotypes and co-expression patterns. PLoS ONE 14 (7). <https://doi.org/10.1371/journal.pone.0219698>

Reumatoïde artritis (RA) | Leeftijd en geslacht | Volksgezondheid en Zorg. (z.d.). VZinfo. Geraadpleegd 26 mei 2026, van <https://www.vzinfo.nl/reumatoide-artritis-ra/leeftijd-en-geslacht>

Rsamtools. (z.d.). Bioconductor. Geraadpleegd 10 juli 2026, van <http://bioconductor.org/packages/Rsamtools/>

Rsubread. (z.d.). Bioconductor. Geraadpleegd 29 mei 2026, van <http://bioconductor.org/packages/Rsubread/>

Xu, H., Yong, L., Gao, X., Chen, Y., Wang, Y., Wang, F., & Hou, X. (2024). CaMK4: Structure, physiological functions, and therapeutic potential. Biochemical Pharmacology, 224, 116204. <https://doi.org/10.1016/j.bcp.2024.116204>

Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). clusterProfiler: An R package for comparing biological themes among gene clusters. Omics: A Journal of Integrative Biology, 16(5), 284-287. <https://doi.org/10.1089/omi.2011.0118>