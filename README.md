# Distill

Variant Pathogenicity Classification in rare mendelian disease (no there is no crazy acronyn involved). 

Random Forest + xgboost + Deep LSTM - based pathogenicity classifier of human genetic variation. Uses curated set of ClinVar, solved UK10K eye mendelian disorders, and gnomAD DNA variants with a broad set of *in silico*, population, gene conseequence, and conservation based metrics and annotations to build highly **specific** variant classification score. 

There are a few minor and one major main distinctions in how this algorithm is built and assessed compared to popular pathogenic scoring systems like REVEL, CADD, and FATHMM.

Minor differences:
- Deep set of putative non-pathogenic mutations, including real WGS from rare mendelian disease
- Rich set of annotations ranging from the scores like CADD and FATHMM, to non-coding annotations like LINSIGHT, to *in silico* consequences from VEP, and to population allelle frequency metrics.
- Ensembl model from three high performing frameworks (Random Forest, xgboost, Keras)

Major difference:
- Training of algorithms under more *realistic* situations that genetic clinics will face for mendelian disorders for singletons and trios. 
  - 250:1 not_pathogenic:pathogenic ratio
