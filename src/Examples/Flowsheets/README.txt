This folder contains example flowsheets using the IDAES platform.

The HDA example demonstrates hydrodealkylation of benzene to toluene via a standard set of unit operations. These operations include a purge and recycle for waste gases and a distillation unit for product recovery, requiring sequential decomposition to initialize the flowsheet and a Translator block to switch property packages used in different parts of the flowsheet.

The methanol synthesis example demonstrates simplified conversion of hydrogen to methanol via a standard set of unit operations. These operations include post-processing prior to a flash unit for product recovery, and implement a single VLE property package for the entire flowsheet.