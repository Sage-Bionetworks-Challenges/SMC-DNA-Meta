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
* The file should include at least a one-line header, indicating the number of pipelines used in the model: ##Pipelines=3
* The file should have the following column names: CHROM, POS, Sample, Predicted, Probability (optional)
* No whitespaces in fields, fields must be tab-delimited (this is also in the spec, but 
it's a common reason for the parser to fail).


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
