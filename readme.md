Readme
========================================================

Date: lipiec 03 2015

### List of files


```
## Michal, you slacker! File(s): cv_analysis.R, cv_data.RData, cv_encodings2.R, cv_encodings.R, encoded_amyloids2.R, encoded_amyloids2.RData, encoded_amyloids.R, encoded_amyloids.RData, figure, final_plot_dat.RData, frame_summarize_functions.RData, multigram_quipt.R, ngram_analysis2.R, ngram_analysis_frame_size.R, only6.csv, pair_dat.RData, performance_plot.R, presentation.bbl, presentation_cp1250.Rnw, presentation.pdf, presentation.Rnw, presentation.tex, properties_graphs.R, property_app, prop_names.RData, randomForestsb2.R, randomForestsb3.R, randomForestsb.R, report2.html, report2.RData, report2.Rmd, report3.html, report3.RData, report3.Rmd, results.RData, SP.png, tableA.csv do not have a proper description.
```

<!-- html table generated in R 3.2.1 by xtable 1.7-4 package -->
<!-- Fri Jul  3 13:45:20 2015 -->
<table border=1>
<tr> <th> File name </th> <th> Description </th> <th> Sourced files </th>  </tr>
  <tr> <td> aa_encodings2.R </td> <td> Amino acid encodings used in the report3. Uses AAIndex database. </td> <td>  </td> </tr>
  <tr> <td> aa_encodings.R </td> <td> Amino acid encodings used in the report2. Uses some older variant AAIndex database (?). </td> <td>  </td> </tr>
  <tr> <td> AA_index_mk2.csv </td> <td> Properties from AAIndex annotated by MK. </td> <td>  </td> </tr>
  <tr> <td> aa_nprop.RData </td> <td> Normalized values of properties from AAIndex. </td> <td>  </td> </tr>
  <tr> <td> aa_properties.R </td> <td> Extract from seqinr package AAIndex data, annotate it partially and save. </td> <td>  </td> </tr>
  <tr> <td> all_sizes.csv </td> <td> Created by ngram_analysis.R, contains important n-grams and their properties in amyloids and non-amyloids. </td> <td>  </td> </tr>
  <tr> <td> amyloid_data.csv </td> <td> Data set from Gasior, Kotulska 2015. </td> <td>  </td> </tr>
  <tr> <td> amyloid_data_extraction.R </td> <td> Extraction of data from amyloid_data.csv, some machine learning. Done for presentation.Rnw </td> <td> randomForestsb2.R </td> </tr>
  <tr> <td> amyloid_fold_res2.RData </td> <td> Results of cross-validation of the amyloid data done for report3. </td> <td>  </td> </tr>
  <tr> <td> amyloid_fold_res.RData </td> <td> Results of cross-validation of the amyloid data done for report2. </td> <td>  </td> </tr>
  <tr> <td> amyloid_neg.fasta </td> <td> Negative (non-amyloid) sequences. </td> <td>  </td> </tr>
  <tr> <td> amyloid_pos.fasta </td> <td> Positive (amyloid) sequences. </td> <td>  </td> </tr>
  <tr> <td> amyloid_present.bib </td> <td> Bibliography for amyloid presentation. </td> <td>  </td> </tr>
  <tr> <td> amyloid_seqs.RData </td> <td> Created by amyloid_data_extraction.R and used in presentation. </td> <td>  </td> </tr>
  <tr> <td> analyze_ngrams_function.R </td> <td> A function that extracts important n-grams and calculate their frequency. </td> <td>  </td> </tr>
  <tr> <td> aoc.R </td> <td> Areas of the competence for classifiers. </td> <td>  </td> </tr>
  <tr> <td> aoc.RData </td> <td> Areas of the competence for classifiers - data for presentation. </td> <td>  </td> </tr>
  <tr> <td> bigger6.csv </td> <td> Important n-grams extracted from sequences longer than 6. </td> <td>  </td> </tr>
  <tr> <td> bigram_analysis.R </td> <td> Analysis of bigrams for presentation. </td> <td>  </td> </tr>
  <tr> <td> bigram_analysis.RData </td> <td> Data from analysis of bigrams for presentation. </td> <td>  </td> </tr>
  <tr> <td> casualRF.R </td> <td> Random forests on amyloid sequences without encoding. </td> <td> read_data.R </td> </tr>
  <tr> <td> ngram_analysis.R </td> <td> Analysis of n-grams extracted from amyloids. The analysis was performed separately for several groups of peptides. </td> <td> read_data.R, analyze_ngrams_function.R </td> </tr>
  <tr> <td> read_data.R </td> <td> Read data from .fasta files </td> <td>  </td> </tr>
  <tr> <td> readme.md </td> <td> Parsed readme.Rmd </td> <td>  </td> </tr>
  <tr> <td> readme.Rmd </td> <td> A template for automatically generated readme. Date and list of file included. </td> <td>  </td> </tr>
  <tr> <td> report1.html </td> <td> Parsed report1.Rmd </td> <td>  </td> </tr>
  <tr> <td> report1.Rmd </td> <td> First report containing results from ngram_analysis.R </td> <td>  </td> </tr>
   </table>

### Repository rules
1. Each file must have a desctiption in readme using the R code in *readme.Rmd*.
2. Use comments frequently.
