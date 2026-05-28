# Transcriptomics van reumatoïde artritis
## 📁 Inhoud/structuur

- `data/raw/` – Ruwe data afkomstig van 8 individuen, waarvan 4 RA hebben, en 4 gezond zijn.  
- `data/processed` - Verwerkte datasets gegenereerd met scripts 
- `scripts/` – Scripts waarin de data geanalyseerd wordt
- `resultaten/` - Figuren zoals grafieken en tabellen
- `bronnen/` - Gebruikte bronnen en AI gebruik
- `README.md` - Het document om de tekst hier te genereren
- `assets/` - Overige documenten voor de opmaak van deze pagina
- `data_stewardship/` - Aantoning van de competentie Beheren niveau I 

---

## Inleiding

Reumatoïde artritis (RA) is een veelvoorkomende autoimmuunziekte. In 2024 waren er zo’n 200.000 Nederlanders met deze aandoening [(Reumatoïde artritis (RA) | Leeftijd en geslacht | Volksgezondheid en Zorg, z.d.)](#vzinfo2026). RA wordt gekenmerkt door synoviale hyperplasie met pannusvorming [(Chetina & Markova, 2019)](#chetina2019)[(Chetina & Markova, 2019)](#chetina2019). Synovium, een type slijmvlies raakt ontstoken, wat lijdt tot de erosie van bot- en kraakbeenweefsel wegens vorming en overname van pannusweefsel [(Gravallese & Monach, 2015)](#gravallese2015). Momenteel focust onderzoek zich op de genen die met T-cellen te maken hebben [(Padyukov, 2022)](#padyukov2022). Er zijn al kleine klinische studies bezig met een mogelijk geneesmiddel, CAR T-celtherapie [(Freeley, 2025)](#freeley2025). 

Transcriptomics is een veelgebruikte studie binnen genetisch en medisch onderzoek. Het transcriptoom, wat hierbij onderzocht wordt, geeft informatie over hoeveelheid expressie in alle genen. Naar het kijken van de functies van genen met hoge expressie, kan een beeld worden geschetst van de basis processen die plaatsvinden bij zieke individuen, in vergelijking met gezonde personen [(Monzó et al., 2025)](#monzó2025)(Monzó et al., 2025)(#monzó2025).

In dit onderzoek wordt een transcriptomics analyse uitgevoerd. De nadruk wordt gelegd op genen waar minder onderzoek naar gedaan is om een dieper begrip te krijgen van alle processen en pathways die betrokken zijn bij RA, met behulp van een VolcanoPlot, GO-analyse en KEGG pathway-analyse.

## Methoden
Referentiegenoom:opzoeken en indexeren (FASTA: fastq. en gtf.)???

KEGG pathway
Genexpressie bepalen aan de hand van de DESeq2 V*versie* package in R
Bepalen verschillende gene ontologies in R






-

-

## 📊 Resultaten

-
`hoi :+) voor code blokken`
-
alt = “Figuur 1. Volcano plot van de verschillen in gesequencete genen van gezonde individuen (n=4) tegenover RA patienten (n=4). Genen n=29407, p-waarde is meegenomen. De x-range loopt van -14 naar 14, aangezien alle data binnen deze punten ligt. Er is gekozen elke 2 waarden op de x-as aan te geven voor overzicht.”

![img](assets/Volcanoplot_RA.png)


## Conclusie

-

## Bronnen
NOG NIET OP VOLGORDE
<a id="chetina2019"></a> Chetina, E. V., & Markova, G. A. (2019). Prospects for the Use of Gene Expression Analysis in Rheumatology. Biochemistry (Moscow), Supplement Series B: Biomedical Chemistry, 13(1), 13-25. <https://doi.org/10.1134/S1990750819010049>
<a id="freeley2025"></a> Freeley, M. (2025). CAR T Cell Therapy for Rheumatoid Arthritis. Clinical Reviews in Allergy & Immunology, 68(1), 100. <https://doi.org/10.1007/s12016-025-09113-7>
<a id="gravallese2015"></a> Gravallese, E. M., & Monach, P. A. (2015). The rheumatoid joint: Synovitis and tissue destruction. In Rheumatology (pp. 768-784). Mosby. <https://doi.org/10.1016/B978-0-323-09138-1.00094-2>
<a id="monzó2025"></a> Monzó, C., Liu, T., & Conesa, A. (2025). Transcriptomics in the era of long-read sequencing. Nature Reviews. Genetics, 26(10), 681-701. <https://doi.org/10.1038/s41576-025-00828-z>
<a id="vzinfo2026"></a> Reumatoïde artritis (RA) | Leeftijd en geslacht | Volksgezondheid en Zorg. (z.d.). VZinfo. Geraadpleegd 26 mei 2026, van <https://www.vzinfo.nl/reumatoide-artritis-ra/leeftijd-en-geslacht>
<a id="padyukov2022"></a> Padyukov, L. (2022). Genetics of rheumatoid arthritis. Seminars in Immunopathology, 44(1), 47-62. <https://doi.org/10.1007/s00281-022-00912-0>

