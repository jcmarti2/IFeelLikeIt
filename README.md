# IFeelLikeIt
this repository contains the code for the final project of CEE 498 SIS.

authors:
  - Juan Carlos Martinez
  - Paul Gharzouzi
  - Jimmy Chang

note: a 'silly' name was chosen for this repository to avoid being found by classmates while the project was ongoing. 

note: to avoid future plagiarism, **this repository will be deleted after final grades are posted on May 21, 2016**.

copyright: the authors of this repository **do not** authorize reproduction of the files included to anyone other than the CEE 498 SIS staff for grading purposes. we will not hesitate reporting plagiarism if we become aware of it.

## directories

###data\
this directory contains the given data used for the IIM portion of the project
  - A.csv: matrix A given in the instructions sheet

###src\
this directory contains the source code used for the IIM and the Power System portions of the project
  - IIM.m: IIM class for holding the properties and functions needed to generate the report results
  - infrastructure_interdependence_analysis.m: script for instantiating the IIM class and generating the report results
  - power_system.m: power_system class holding the properties and functions needed to generate the report results
  - power_system_design.m: script for instantianting the power_system class and generating the report results
  
###report\
this directory contains the final pdf version of the report, the LaTex documents, and the figures
  - report.pdf: final version of the report
  - \latex: latex files to generate report
  - \figures: directory containing the figures included in the report

## instructions for generating results
the source included in this repository aims to answer the questions presented in the project instructions sheet. for this purpose, we have developed scripts that answer the questiosn required, one by one. these files are 
  - src\infrastructure_interdependence_analysis.m
  - src\power_system_desing.m
  
these two files are organized as follows: 1) user input section, 2) no user input section. 

the user input section contains a series of boolean variables assigned the valeus true or false. these are named after the questions to be answered, in the format q_X_YY, where q refers to 'question', X is the section number (2 or 4), and YY is the question number. the question numbering given in the instructions sheet is preserved.

to answer a specifc question, set the variable q_X_YY desired to true and run the code. please read any comments included in the files for more specific details.
