# Benchmarking of the scientific software

**This project contains several steps**:
1. Preparation of the testing library
2. Running the programs under study on a test dataset and obtaining the results
3. Preparation of the resulting files
4. Statistical analysis
5. Visualization of the statistics

**Scenario of the project**:
1) Quality assessment [Q1 = internal data];
2) Benchmarking with other "downloadable, free for academic use and popular" software (IgBLAST, RosettaAntibody, ImmunediveRsity, ect.) [including Q2 = external data];
3) As our tool is applied just for annotation => to recommend a pypline for handling antibody sequences;
4) Testing the pypline

**Metrics used for assessment:**

*Error_rate* = probability of software to make incorrect predictions of ending+-starting amino acids for a testing data with respect to germline db

*Sensitivity(Recall)* = **TP/(TP+FN)**, defines how precise the current software can predict specific region (for ex. FR3-fragment) on the testing data with respect to germline db

*Specificity* = **TN/(TN+FP)**, defines how well current software can separate specific region from other regions

*Accuracy* = **(TP+TN)/(TP+TN+FP+FN)**, defines ability of current software to make a correct prediction of specific region among total number of predictions for the testing data with respect to germline db

*Precision* = **TP/(TP+FP)**

*F1-score* = **2TP/(2TP+FP+FN)**

To estimate performance and resource consumption: exec time, CPU time, max amount of RAM.

The **goal** is to end-up with ROC-curve summarizing all metrics to define more accurate and resource saving software.

**Other criteria:**
- ease of installation
- ease of performing
- user-guide
- ability to apply own germline db
- alignment algorithm
- etc

*This repo hosted by Dana Zotikova igCat Project*

## Some results:
![](/home/dana/Documents/dana/py_projects/igcat/images/Framework Region 1.png)
![](/home/dana/Documents/dana/py_projects/igcat/images/Framework Region 2.png)
![](/home/dana/Documents/dana/py_projects/igcat/images/Screenshot from 2022-09-27 20-31-13.png)
![](/home/dana/Documents/dana/py_projects/igcat/images/Screenshot from 2022-09-27 20-31-59.png)