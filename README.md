# Transcriptomics van reumatoïde artritis
## 📁 Inhoud/structuur
- `assets` - Overige documenten voor de opmaak van deze pagina.
- `bronnen` - Gebruikte bronnen en AI gebruik.
- `data processed` – Bewerkte data voor helder overzicht van resultaten en uitvoeren van analyses. 
- `data raw` – Ruwe data afkomstig van 8 individuen, waarvan 4 RA hebben, en 4 gezond zijn.  
- `data_stewardship` - Aantoning van de competentie Beheren niveau I 
- `scripts` – Scripts waarin de data geanalyseerd wordt.
- `resultaten` - Figuren zoals grafieken en tabellen.
- `README.md` - Het document om de tekst hier te genereren.

---

## Inleiding

Reumatoïde artritis (RA) is een veelvoorkomende autoimmuunziekte. In 2024 waren er zo’n 200.000 Nederlanders met deze aandoening [(Reumatoïde artritis (RA) | Leeftijd en geslacht | Volksgezondheid en Zorg, z.d.)](#vzinfo2026). RA wordt gekenmerkt door synoviale hyperplasie met pannusvorming [(Chetina & Markova, 2019)](#chetina2019). Synovium, een type slijmvlies raakt ontstoken, wat lijdt tot de erosie van bot- en kraakbeenweefsel wegens vorming en overname van pannusweefsel [(Gravallese & Monach, 2015)](#gravallese2015). Momenteel focust onderzoek zich op de genen die met T-cellen te maken hebben [(Padyukov, 2022)](#padyukov2022). Er zijn al kleine klinische studies bezig met een mogelijk geneesmiddel, CAR T-celtherapie [(Freeley, 2025)](#freeley2025). 

Transcriptomics is een veelgebruikte studie binnen genetisch en medisch onderzoek. Het transcriptoom, wat hierbij onderzocht wordt, geeft informatie over hoeveelheid expressie in alle genen. Naar het kijken van de functies van genen met hoge expressie, kan een beeld worden geschetst van de basis processen die plaatsvinden bij zieke individuen, in vergelijking met gezonde personen [(Monzó et al., 2025)](#monzó2025)(Monzó et al., 2025)(#monzó2025).

In dit onderzoek wordt een transcriptomics analyse uitgevoerd. De nadruk wordt gelegd op genen waar minder onderzoek naar gedaan is om een dieper begrip te krijgen van alle processen en pathways die betrokken zijn bij RA, met behulp van een VolcanoPlot, GO-analyse en KEGG pathway-analyse.

## Methoden

Er wordt ingezoomd op de transcriptomics van RA, om genen op te sporen met onderzoekspotentie [(figuur 1)](). De ruwe sequencing data in FASTQ bestanden is afkomstig van 8 vrouwen, 4 met RA (leeftijden 54-66), vastgesteld voor >12 maanden en positief getest op autoantistoffen ACPA. En een negatief geteste controlegroep van 4 (leeftijden 15-42). De verkregen data is uitgewerkt in R V4.5.2 in dit [script](scripts/script_casus_transcriptomics_RA.R).

Bij data mapping wordt een index gebouwd met `BiocManager` V1.30.27 package `Rsubread` V2.24.0. Hierin worden FASTA bestanden van het RefSeq referentiegenoom [GRCh38.p14](<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/>) van de NCBI genome database gemaakt. Het alignen tot BAM bestanden is in paired-end. Rsamtools V2.26.0 sorteert en indexeert de BAM files. De Count Matrix en BAM files matrix wordt gemaakt aan de hand van een GTF annotatiebestand, nogmaals het NCBI GRCh38.p14 genoom. De counts worden geïsoleerd en als csv opgeslagen.

Een VolcanoPlot wordt gemaakt voor genexpressie bepaling. Hiervoor zijn de `BiocManager` packages `DESeq2` V1.50.2 en `EnhancedVolcano` V1.28.2 nodig. Een DESeq dataset wordt aangemaakt. De VolcanoPlot bevat de log2FoldChange en gecorrigeerde P-waarde <0.05.

Pathway analyses vereisen `pathview` V1.50.0 en `KEGGREST` V1.50.0. Een enriched KEGG stelt het juiste organisme vast en de pathview wordt uitgevoerd op het gen van interesse met een log2FoldChange vector.

Bepalen van verschillende Gene Ontologies gaat via `goseq` V1.62.0, `org.Hs.eg.db` V3.22.0 en `AnnotationDbi` V1.72.0. Hier wordt een gecorrigeerde P-waarde van <0.05 aangehouden en het hg38 genoom wordt gebruikt.


<p align="center">
  <img src="assets/Flowchart_Methode_RA.png" alt="Flowchart Methode" width="600"/>
</p>
*Figuur 1. Flowchart van de gebruikte methode. Een referentiegenoom wordt geïndexeerd. RNA seq data en Index worden gemapt naar BAM files. De reads van deze files worden samen met de GTF annotatie geteld.*

## 📊 Resultaten

-

-
<p align="center">
  <img src="resultaten/Volcanoplot_RA.png" alt="Resultaten Volcanoplot" width="600"/>
</p>
*Figuur 2. Volcano plot van de verschillen in gesequencete genen van gezonde individuen (n=4) tegenover RA patienten (n=4). Genen n=29407, p-waarde is meegenomen. De x-range loopt van -14 naar 14, aangezien alle data binnen deze punten ligt. Er is gekozen elke 2 waarden op de x-as aan te geven voor overzicht.*



<p align="center">
  <img src="resultaten/hsa04662.pathview_RA.png" alt="Resultaten Pathway Analyse" width="600"/>
</p>
Figuur 3. KEGG pathway van het ... gen.

## Conclusie

-

-

-


## Bronnen
**NOG NIET OP VOLGORDE**

<a id="chetina2019"></a> Chetina, E. V., & Markova, G. A. (2019). Prospects for the Use of Gene Expression Analysis in Rheumatology. Biochemistry (Moscow), Supplement Series B: Biomedical Chemistry, 13(1), 13-25. <https://doi.org/10.1134/S1990750819010049>

<a id="freeley2025"></a> Freeley, M. (2025). CAR T Cell Therapy for Rheumatoid Arthritis. Clinical Reviews in Allergy & Immunology, 68(1), 100. <https://doi.org/10.1007/s12016-025-09113-7>

<a id="gravallese2015"></a> Gravallese, E. M., & Monach, P. A. (2015). The rheumatoid joint: Synovitis and tissue destruction. In Rheumatology (pp. 768-784). Mosby. <https://doi.org/10.1016/B978-0-323-09138-1.00094-2>

<a id="monzó2025"></a> Monzó, C., Liu, T., & Conesa, A. (2025). Transcriptomics in the era of long-read sequencing. Nature Reviews. Genetics, 26(10), 681-701. <https://doi.org/10.1038/s41576-025-00828-z>

<a id="vzinfo2026"></a> Reumatoïde artritis (RA) | Leeftijd en geslacht | Volksgezondheid en Zorg. (z.d.). VZinfo. Geraadpleegd 26 mei 2026, van <https://www.vzinfo.nl/reumatoide-artritis-ra/leeftijd-en-geslacht>

<a id="padyukov2022"></a> Padyukov, L. (2022). Genetics of rheumatoid arthritis. Seminars in Immunopathology, 44(1), 47-62. <https://doi.org/10.1007/s00281-022-00912-0>

Homo sapiens genome assembly GRCh38.p14. (z.d.). NCBI. Geraadpleegd 28 mei 2026, van <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/>

<a id="???"></a>

<a id="???"></a>

<a id="???"></a>