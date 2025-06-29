Benchmark Datasets
================

| Dataset Name | Domain | Description | Samples | Features | Groups | Dataset Reference | Main Experimental Reference | Methods applied† | Best Method (Reference) | Best Metrics (Reference) | Best Method (Experiments) | Best Metrics (Experiments) | Dataset available at |
|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
| **meats** | Food | Homogenized raw meat samples from 5 animal species (Beef, Chicken, Lamb, Pork, Turkey), measured by near infrared spectroscopy | 231 | 1,050 | 5 | McElhinney et al. (1999) | Anderlucci et al. (2022) | AffPr, GMM, HDDC-B, HDDC-C, KM, PAM, RP, RPEClu, SpeCl, Ward | RPEClu | ARI=0.32 | IS + RP (GMM) + LCE | ARI=0.69 | R package `RPEClust` |
| **ALLAML** | Genomics | Acute lymphoblastic leukemia (ALL) and acute myeloid leukemia (AML) mass spectrometry data | 72 | 7,129 | 2 | Fodor (1997) | Guo et al. (2018) | KM + varied UFS | KM + DGUFS | ACC=78.4, NMI=19.1 | N/A | N/A | <https://jundongl.github.io/scikit-feature/datasets.html> |
| CNS tumors | Genomics | Embryonal CNS tumors: medulloblastomas, PNET, rhabdoid tumors, gliomas, normal cerebellum | 48 | TBD | 5 | Pomeroy et al. (2002) | Monti et al. (2003) | CC, HC, SOM | HC + CC | RAND=0.549 | N/A | N/A |  |
| **leukemia** | Genomics | Bone marrow samples from AML, T-lineage ALL, and B-lineage ALL leukemia patients | 38 | 7,129 | 3 | Golub et al. (1999) | Monti et al. (2003) | CC, HC, SOM | SOM + CC | RAND=0.721 | N/A | N/A | R package `golubEsets` |
| **lung_cancer** | Genomics | Human lung carcinomas: adenocarcinoma, squamous cell carcinoma, carcinoid, normal lung | 203 | 3,313 | 4 | Bhattacharjee et al. (2001) | Monti et al. (2003) | CC, HC, SOM | HC + CC | ARI=0.310 | IS + RP (GMM) + LCE | ARI=0.41 | <https://homes.di.unimi.it/valenti/DATA/MICROARRAY-DATA/Bhattacharjee-Lung> |
| **lymphoma** | Genomics | Lymphoma gene expression: DLBCL, FL, CLL | 62 | 4,026 | 3 | Alizadeh et al. (2000) | Anderlucci et al. (2022) | AffPr, GMM, HDDC-B, HDDC-C, KM, PAM, RP, RPEClu, SpeCl, Ward | RPEClu | ARI=1.00 | N/A | N/A | R package `spls` |
| Normal tissues | Genomics | 13 normal tissue types: breast, prostate, lung, colon, etc. | 99 | TBD | 13 | Ramaswamy et al. (2001) | Monti et al. (2003) | CC, HC, SOM | HC + CC | RAND=0.457 | N/A | N/A |  |
| Novartis multi-tissue | Genomics | Tissue samples from 4 cancer types: breast, prostate, lung, colon | 103 | TBD | 4 | Su et al. (2002) | Monti et al. (2003) | CC, HC, SOM | HC + CC | RAND=0.921 | N/A | N/A |  |
| **prostate_ge** | Genomics | Prostate cancer microarray dataset (tumor and normal tissue) | 102 | 5,966 | 2 | Singh et al. (2002) | Guo et al. (2018) | KM + varied UFS | KM + DGUFS | ACC=65.3, NMI=6.59 | IS + RP (KM) + LCA | ACC=75.5, NMI=19.7 | <https://jundongl.github.io/scikit-feature/datasets.html> |
| St. Jude leukemia | Genomics | Pediatric leukemia: 6 subtypes (T-ALL, E2A-PBX1, BCR-ABL, TEL-AML1, MLL, hyperdiploid\>50) | 248 | TBD | 6 | Yeoh et al. (2002) | Monti et al. (2003) | CC, HC, SOM | CC + HC | RAND=0.948 | N/A | N/A |  |
| CBIR | Images | Content-based image retrieval | 1,545 | 183 | 8 | Dy et al. (1999) | Fern & Lin (2003) | KM, SpeCl, RP, FW, CAS, JC, CH, CSPA | CH + CSPA | NMI=0.341 | N/A | N/A |  |
| **COIL20** | Images | Columbia Object Image Library (20 objects, multi-view grayscale images) | 1,440 | 1,024 | 20 | Nene et al. (1996) | You et al. (2023) | KM + varied UFS | KM + NNSE | ACC=71.8, NMI=76.2 | N/A | N/A | <https://jundongl.github.io/scikit-feature/datasets.html> |
| JAFFE | Images | Facial expression images of Japanese female subjects, used for facial expression recognition and clustering tasks | 213 | 676 | 10 | Lyons et al. (1999) | You et al. (2023) | KM + varied UFS | KM + QS | ACC=82.4, NMI=87.4 | N/A | N/A |  |
| UMIST | Images | Face images from 20 individuals | 575 | 644 | 20 | Graham and Allinson (1998) | Guo et al. (2018) | KM + varied UFS | KM + DGUFS | ACC=57.1, NMI=74.4 | N/A | N/A |  |
| **warpAR10P** | Images | Subset of AR face images from 10 people, warping applied (multi-view grayscale) | 130 | 2,400 | 10 | Li, J., Wang, K., et al. | You et al. (2023) | KM + varied UFS | KM + NNSE | ACC=51.9, NMI=39.4 | N/A | N/A | <https://jundongl.github.io/scikit-feature/datasets.html> |
| **warpPIE10P** | Images | Subset of PIE (Pose, Illumination, and Expression) face images from 10 people, warping applied (multi-view grayscale) | 210 | 2,420 | 10 | Li, J., Wang, K., et al. | You et al. (2023) | KM + varied UFS | KM + NNSE | ACC=58.1, NMI=52.2 | IS + RP (GMM) + LCE | ACC=71.4, NMI=78 | <https://jundongl.github.io/scikit-feature/datasets.html> |
| ISOLET6 | Speech | Subset (6 classes) of the ISOLET dataset, spoken letter recognition, high-dimensional features | 1,440 | 617 | 6 | Blake & Merz (1998) | Fern & Lin (2008) | KM, SpeCl, RP, FW, CAS, JC, CH, CSPA | CSPA + CAS | NMI=0.851 | N/A | N/A |  |

Datasets in bold are included in `diceRplus`.

† Methods:  
- AffPr: Affinity propagation  
- CAS: Cluster and Select  
- CC: Consensus Clustering  
- CH: Convex Hull  
- CSPA: Cluster-based Similarity Partitioning Algorithm  
- DGUFS: Dependence Guided Unsupervised Feature Selection  
- FW: Feature Weighting  
- GMM: Gaussian Mixture Model  
- IS: Infinite Feature Selection - NNSE: Neural Networks embedded
Self-Expression  
- HC: Hierarchical Clustering  
- HDDC-B: High-dimensional data clustering (BIC method)  
- HDDC-C: High-dimensional data clustering (Cattel’s scree test)  
- JC: Joint Criterion  
- KM: K-means  
- RP: Random projection  
- RPEClu: Random Projection Ensemble Clustering - PAM: Partition around
Medoids  
- QS: Quick and Robust Feature Selection  
- SpeCl: Spectral Clustering  
- SOM: Self-Organizing Map  
- UFS: Unsupervised Feature Selection  
- Ward: Hierarchical clustering (Ward’s method)

### References

- Alizadeh, A.A., Eisen, M.B., Davis, R.E., Ma, C., Lossos, I.S.,
  Rosenwald, A., et al. (2000). Distinct types of diffuse large B-cell
  lymphoma identified by gene expression profiling. *Nature*,
  **403**(6769), 503–511.

- Anderlucci, L., Fortunato, F., Montanari, A. (2022). High-dimensional
  clustering via random projections. *Journal of Classification*,
  **39**(2), 191–216.

- Bhattacharjee, A., Richards, W.G., Staunton, J., Li, C., Monti, S.,
  Vasa, P., et al. (2001). Classification of human lung carcinomas by
  mRNA expression profiling reveals distinct adenocarcinomas
  sub-classes. *Proceedings of the National Academy of Sciences*,
  **98**(24), 13790–13795.

- Blake, C.L., Merz, C.J. (1998). UCI repository of machine learning
  databases. University of California, Irvine, CA.

- Dy, J., Brodley, C.E., Kak, A., Shyu, C., Broderick, L.S. (1999). The
  customized queries approach to CBIR using EM. Proceedings of IEEE
  Conference on Computer Vision and Pattern Recognition, 400–406.

- Fern, X.Z., Lin, W. (2008). Cluster ensemble selection. *Statistical
  Analysis and Data Mining*, **1**(3), 128–141.

- Fodor, S.P.A. (1997). Massively parallel genomics. *Science*,
  **277**(5324), 393.

- Graham, D. B., Allinson N. M. (1988). Characterising virtual
  eigensignatures for general purpose face recognition. *Face
  Recognition* pp. 446–456, 1998.

- Golub, T.R., Slonim, D.K., Tamayo, P., Huard, C., Gaasenbeek, M.,
  Mesirov, J.P., et al. (1999). Molecular classification of cancer:
  Class discovery and class prediction by gene expression monitoring.
  *Science*, **286**(5439), 531–537.

- Guo, J., Zhu, W. (2018). Dependence Guided Unsupervised Feature
  Selection. *Proceedings of the Thirty-Second AAAI Conference on
  Artificial Intelligence (AAAI-18)*, **32**, 2232–2239.

- Li, J., Wang, K., et al. (scikit-feature). Dataset available at:
  <https://jundongl.github.io/scikit-feature/datasets.html>

- Lyons, M.J., Budynek, J., Akamatsu, S. (1999). Automatic
  classification of single facial images. *IEEE Transactions on Pattern
  Analysis and Machine Intelligence*, **21**(12), 1357–1362.

- McElhinney, J., Downey, G., Fearn, T. (1999). Chemometric processing
  of visible and near infrared reflectance spectra for species
  identification in selected raw homogenised meats. *Journal of Near
  Infrared Spectroscopy*, **7**, 145–154.

- Monti, S., Tamayo, P., Mesirov, J.P., Golub, T.R. (2003). Consensus
  clustering: A resampling-based method for class discovery and
  visualization of gene expression microarray data. *Machine Learning*,
  **52**(1–2), 91–118.

- Nene, S.A., Nayar, S.K., Murase, H., Columbia Object Image Library
  (COIL-20), Technical Report, 1996.

- Pomeroy, S.L., Tamayo, P., Gaasenbeek, M., Sturla, L.M., Angelo, M.,
  McLaughlin, M.E., et al. (2002). Gene expression-based classification
  and outcome prediction of central nervous system embryonal tumors.
  *Nature*, **415**(6870), 436–442.

- Ramaswamy, S., Tamayo, P., Rifkin, R., Mukherjee, S., Yeang, C.-H.,
  Angelo, M., Ladd, C., et al. (2001). Multi-class cancer diagnosis
  using tumor gene expression signatures. *Proceedings of the National
  Academy of Sciences*, **98**(26), 15149–15154.

- Singh, D., Febbo, P.G., Ross, K., Jackson, D.G., Manola, J., Ladd, C.,
  et al. (2002). Gene expression correlates of clinical prostate cancer
  behavior. *Cancer Cell*, **1**(2), 203–209.

- Su, A.I., Cooke, M.P., Ching, K.A., Hakak, Y., Walker, J.R.,
  Wiltshire, T., Orth, A.P., Vega, R.G., et al. (2002). Large-scale
  analysis of the human and mouse transcriptomes. *Proceedings of the
  National Academy of Sciences*, **99**(7), 4465–4470.

- Yeoh, E. J., Ross, M.E., Shurtleff, S.A., Williams, W.K., Patel, D.,
  Mahfouz, R., Behm, F.G., et al. (2002). Classification, subtype
  discovery, and prediction of outcome in pediatric acute lymphoblastic
  leukemia by gene expression profiling. *Cancer Cell*, **1**(2),
  133–143.

- You, M., Yuan, A., He, D., Li, X. (2023). Unsupervised Feature
  Selection via Neural Networks and Self-Expression with Adaptive Graph
  Constraint. *Pattern Recognition*, **135**, 109173.
