# VPaC

**V**ariant **Pa**thogenicity **C**lassification

Random Forest - based pathogenicity classifier of human genetic variation. Uses curated set of ClinVar and solved UK10K eye mendelian disorders with a broad set of *in silico*, population, gene conseequence, and conservation based metrics and annotations to build highly specific variant classification score. 

There are three main distinctions in how this algorithm is built and assessed compared to popular pathogenic scoring systems like REVEL, CADD, and FATHMM.

1. A much deeper set of putative non-pathogenic mutations
2. A highly rich set of annotations ranging from the scores like CADD and FATHMM, to non-coding annotations like LINSIGHT, to *in silico* consequences from VEP, and to population allelle frequency metrics. 
3. Asssement of algorithms under more *realistic* situations that genetic clinics will face for mendelian disorders for singletons and trios. 