# Amendment 1 to locked test v2 (written before any score, AUC or group comparison was computed)
Observation: the GSE255660 series matrix holds background-substituted signal intensities on a linear scale,
without between-array normalisation; spots below detection are blank. The locked protocol assumed normalised values.
Pre-specified handling, applied identically to all arrays and all miRNAs:
 1. Global normalisation (3D-Gene convention): each array is scaled so that the median of its detected signals equals 25.
 2. Undetected (blank) entries are set to half of the smallest normalised detected signal on that array.
 3. Values are log2-transformed.
Sensitivity analysis: the same scores computed without step 1 (raw intensities, steps 2-3 only).
Detection rate of each signature miRNA is reported per group.
Nothing else in locked test v2 is changed.
