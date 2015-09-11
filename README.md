# SMC-DNA-Meta
================================================

Tools for participants in the ICGC-TCGA DREAM SMC-DNA Meta challenge

You will find here the tools and step-by-step instructions for uploading result files 
for participating in the The ICGC-TCGA DREAM Somatic Mutation Calling Meta-pipeline Challenge (SMC-DNA Meta)
(referred herein as The Challenge). You must first
[join the Challenge](https://www.synapse.org/#!Synapse:syn4258642/wiki/231962)
before you can submit.


# Submission steps

1. Obtain and install submission tool.
2. Create a private working project.
3. Run the validator on your submission file.
4. Upload and submit the entity to the evaluation


## Obtain and install submission tool
We provide a command line program to help validate and submit files to the contest.
First, you will need to install the program.
There are a few dependencies that are needed in order for the submission program to work. 
These package help to parse VCF files and contact the Synapse servers.

You can use the standard Python ‘pip’ installer to get these dependencies
```
pip install PyVCF
pip install synapseclient
```

or, to install in your home directory:
```
pip install --user PyVCF
```

Get a copy of the SMC-DNA Meta contestant toolkit: 
https://github.com/Sage-Bionetworks/SMC-DNA-Meta


```
git clone https://github.com/Sage-Bionetworks/SMC-DNA-Meta.git
```


##Create a private working project
If you don’t already have a personal working project, go to the front page of Synapse and 
use the ‘Create Project’ widget. Please make sure that this folder is not visible to the 
public. There is more [info available](https://www.synapse.org/#!ProjectsHome:0) if you 
want to find out about Synapse projects.


##Run the validator on your submission file


The ‘dream_vcf.py’ program will validate your submission format.

```
./dream_vcf.py validate <submitted_file>
```


####Submission format
For all submissions
* The file must be uncompressed (not ending in .gz)
* The file should include at least a 14-line header, indicating the list of pipelines associated with each 4 synthetic tumours and 10 real tumours:
```
  ##synthetic_1_Pipelines=X1,X2
  ##synthetic_2_Pipelines=X3,X4
  ##synthetic_3_Pipelines=X5,X6
  ##synthetic_4_Pipelines=X7,X8
  ##CPCG0100_Pipelines=X9,X10
  ##CPCG0183_Pipelines=X11,X12
  ##CPCG0184_Pipelines=X13,X14
  ##CPCG0196_Pipelines=X15,X16
  ##CPCG0235_Pipelines=X17,X18
  ##PCSI0023_Pipelines=X19,X20
  ##PCSI0044_Pipelines=X21,X22
  ##PCSI0046_Pipelines=X23,X24
  ##PCSI0048_Pipelines=X25,X26
  ##PCSI0072_Pipelines=X27,X28
```

* The file should have the following column names: CHROM, POS, Sample, Predicted, Probability (optional),
```
  #CHROM  POS Sample  Predicted Probability
```
  followed by the calls, specified in order by chromosome, position, and sample:
```
  1    12345    synthetic_1    1    0.758
  1    12345    synthetic_2    1    0.758
  1    12345    synthetic_3    1    0.758
  1    12345    synthetic_4    1    0.758
  1    12345    CPCG0100    1    0.758
  1    12345    CPCG0183    1    0.758
  1    12345    CPCG0184    1    0.758
  1    12345    CPCG0196    1    0.758
  1    12345    CPCG0235    1    0.758
  1    12345    PCSI0023    1    0.758
  1    12345    PCSI0044    1    0.758
  1    12345    PCSI0046    1    0.758
  1    12345    PCSI0048    1    0.758
  1    12345    PCSI0072    1    0.758
```

* No whitespaces in fields, fields must be tab-delimited (this is also in the spec, but 
it's a common reason for the parser to fail).
* The number of pipelines associated with each tumour list of pipelines MUST be the same.


##SynapseClient Login
You can either pass you name and password to the dream_submit program via the command 
line, or to can follow the instructions setting up an 
[authentication config file](https://www.synapse.org/#!Synapse:syn1768504/wiki/56068)

You can also visit [the settings page](https://www.synapse.org/#!Settings:0) to get an API 
key, then set the environmental variables:

```
export SYNAPSE_APIKEY=(really long string you copied from the setting page)
export SYNAPSE_EMAIL=(your email address)
```


##Upload and submit the entity to the evaluation
