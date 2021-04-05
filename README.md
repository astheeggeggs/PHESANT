# Customised PHESANT - PHEnome Scan ANalysis Tool.

![plot](./docs/phesant.png)

Run a phenome scan and write the resultant phenotypes to disk. Much of the code is identical to [PHESANT](https://github.com/MRCIEU/PHESANT), and we refer users to that repository for running analyses and creating visualisations.

## General requirements

R with the R packages: data.table, optparse and MASS.

## Citing this project

Please cite:

Millard LAC, Davies NM, Gaunt TR, Davey Smith G, Tilling K. PHESANT: a tool for performing automated phenome scans in UK Biobank. bioRxiv (2017)

## Running a phenome scan

A phenome scan is run using `WAS/phenomeScan.r`. We have modified the original code so that now no tests are performed, we simply clean the phenotype information and provide a .tsv file which can be used for downstream analysis. In addition, we provide some simple scripts to summarise the resultant phenotype data.

The modifield custom PHESANT phenome scan processing pipeline is illustrated in the figure [here](pipeline.pdf). The original pipeline, together with extensive detailed descriptions can be found in the above paper.

The phenome scan is run with the following command:

```bash
cd WAS/

Rscript phenomeScan.r \
--phenofile=<phenotypesFilePath> \
--variablelistfile="../variable-info/outcome-info.tsv" \
--datacodingfile="../variable-info/data-coding-ordinal-info.csv" \
--resDir=<resultsDirectoryPath> \
--userId=<userIdFieldName>
```

### Required arguments

Arg | Description
-------|--------
phenofile 		| Comma separated file containing phenotypes. Each row is a participant, the first column contains the user id and the remaining columns are phenotypes. Where there are multiple columns for a phenotype these must be adjacent in the file. Specifically for a given field in Biobank the instances should be adjacent and within each instance the arrays should be adjacent. Each variable name is in the format 'x[varid]\_[instance]\_[array]' (we use the prefix 'x' so that the variable names are valid in R).
variablelistfile 	| Tab separated file containing information about each phenotype, that is used to process them (see below).
datacodingfile 		| Comma separated file containing information about data codings (see below).
resDir 			| Directory where you want the results to be stored.

### Optional arguments
Arg | Description
-------|--------
userId                  | User id column as in the traitofinterestfile and the phenofile (default: userId).
partIdx			| Subset of phenotypes you want to run (for parallelising).
numParts		| Number of subsets you are using (for parallelising).

The numParts and partIdx arguments are both used to parallelise the phenome scan. Note that this acts on the original phenotype file, so the resultant pieces will be of differing sizes as phenotypes are chosen for inclusion or split up according the rules outlined in the [pipeline](pipeline.pdf). E.g. setting numParts to 5 will divide the set of phenotypes into 5 (rough) parts and then partIdx can be used to call the phenome scan on a specific part (1-5).

#### Data coding file

Data codes define a set of values that can be assigned to a given field. A data code can be assigned to more than one variable, which is why we use
a separate file describing the necessary information for each data code. For example, there are several fields about diet that have data code [100009](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100009).

The data coding file should have the following columns:

1. dataCode - The ID of the data code.
2. ordinal - Whether the field is ordinal (value 1) or not (value 0). This field is only used for fields of the categorical (single) field type. Value -1 denotes this is not needed because the field is binary. 
3. ordering - Any needed corrections for the numeric ordering of a data codes specified by Biobank. This field is only used for data codes specified as ordinal in the ordinal column. For example, data code [100001](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100001) has values half, 1 and 2+ coded as 555, 1 and 200, respectively. We need the 'half' value to be less than the '1' value, so we change the order to '555|1|200'. NB: if this column is used then and any value is not included then this value is set to NA (i.e. this field can be used to remove and reorder values at the same time).
4. reassignments - Any value changes that are needed. For example, in data code [100662](http://biobank.ctsu.ox.ac.uk/showcase/coding.cgi?id=100662), the values
7 and 6 may be deemed equal (both representing 'never visited by friends/family' so we can set '7=6' to assign the value 6 to all participants with the value 7.
5. default_value - A default value assigned to all participants with no value for the field, but with a value for field stated in `default_value_related_field` column below. This is used where a category is not explicitly stated in the field but 
instead needs to be determined by looking at whether another field has a value. Typically, this occurs where there is no category for 'none' in a questionnaire field, because participants were told they did not have to mark 'none' but could instead leave it blank 
(see for example section 5.3 in the [24 hour diet questionnaire manual](http://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=118240)). Hence, we assume that if they completed the questionnaire and have not ticked a value, then the value is 'none'. See default value example below.
6. default_value_related_field - The field used to determine which participants are assigned the default value. All participants with a value in the field stated here, and with no value for a field with this data code, are assigned the default value stated in `default_value`.

##### Example of default value

In the data code information file we specify `default_value=0` and `default_value_related_field=20080` for data code 100006. 
Field [100200](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=100200), for example, has data code 100006. 
Therefore all participants with a value for field [20080](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=20080), but with no value in field 100200, are assigned value 0 for field 100200.
Intuitively, all participants who have answered the 24-hour recall diet questionnaire have a value in field 20080, and of these, we assume that those with no value for field 100200 have opted
for 'none' implicitly, by not ticking any option.

#### Variable information file

This file was initially the UK Biobank data dictionary, which can be downloaded from the UK Biobank website [here](http://biobank.ctsu.ox.ac.uk/~bbdatan/Data_Dictionary_Showcase.csv).
This data dictionary provides the following set of information about fields, used in this phenome scan tool:

1. ValueType column - the field type, either 'Integer', 'Continuous', 'Categorical single', 'Categorical multiple', or a few others we do not use.
2. Three Cat_ID and three Cat_Title columns - the three levels of the category hierarchy, that can be seen [here](http://biobank.ctsu.ox.ac.uk/showcase/label.cgi).
3. FieldID column - We use this to match the variable in our biobank data file to the correct row in this TSV file.
4. Field column -  The name of the field.

The variable information file also has the following columns that we have added, to provide additional information used in the phenome scan. For the parsed phenotype information used for GWAS in our group, we use the EXCLUDED column to manually curate the phenotypes of interest - see below.

1. TRAIT_OF_INTEREST - Specifies any field that represents the trait of interest (set this column to 'YES'). This is a marker so that after the phenome scan is run we can use these
results as validation only (e.g. a pheWAS of the BMI FTO SNP would expect the BMI phenotypes to show high in the results ranking), i.e. they do not contribute to the multiple testing burden. We have set this up for BMI, so have marked BMI/weight fields as the trait of interest - you will need to change this for your particular trait of interest. 
For categorical multiple fields you may want to mark the whole field as denoting the trait of interest (e.g all [cancers](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=20001)), or just a specific value (e.g. a particular type of cancer). To do the former set this column to YES, and to do the latter specify each particular value in this field separated by a bar, i.e. 'VALUE1|VALUE2'.
2. EXCLUDED - Phenotypes we apriori decide to exclude from the phenome scan. Any field with a value in this field is excluded, and we state a code that describes the reason we exclude a variable (for future reference). Our codes and reasons are as follows (of course for your phenome scan you can add others as you would like): 
 - YES-ACE: "Assessment center environment" variables that do not directly describe the participant.
 - YES-AGE: Age variables.
 - YES-ASSESSMENT-CENTRE: The assessment centre the participant attended.
 - YES-BIOBANK-SUGGESTED-VARIABLE: Variables not initially included in our data request, but that biobank suggested we receive.
 - YES-CAT-SIN-MUL-VAL: Fields that were 'Categorical single' types but had multiple values (several arrays in an instance). We do not deal with these currently so remove from the phenome scan.
 - YES-GENETIC: Genetic description variables.
 - YES-SENSITIVE: Variables not received from Biobank because they are sensitive so have more restricted access.
 - YES-SEX: Sex fields.
 - YES-NEALELAB: Removing these fields in addition to the above removals defines our manually curated phenotype collection (assumming access to the entire collection of phenotypes).
3. CAT_MULT_INDICATOR_FIELDS - every categorical multiple field must have a value in this column. 
The value describes which set of participants to include as the negative examples, when a binary variable is created from each value (see above cited paper for more information). 
The positive examples for a value `v` in this categorical multiple field are simply the people with this particular value. However the negative values can be determined in three ways:
 - ALL - Include all participants (any participant without value `v` is assigned `FALSE`, except those with a value denoting missingness (i.e. value is <0)). 
 - NO_NAN - Included only those who have at least one value for this field (assign `FALSE` to any participant with at least one value for this field, where these values do not include value `v`, and also do not include a value denoting missingness (i.e. value is <0).
 - fieldID - Include only those who have a value for another field, with ID `fieldID` (assign `FALSE` to any participant without value `v` and without a value denoting missingness (i.e. value is <0) and with a value in this other field with ID `fieldID`).
4. CAT_SINGLE_TO_CAT_MULT - Specifies fields that have the categorical single field type (as specified by UK Biobank) but that we actually want to treat as categorical multiple. State YES in this column to change this. To also convert the instances to arrays, state YES-INSTANCES. For example, field [20107](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=20107) has illnesses stored in 10 arrays (for each instance) so it makes sense to treat this as a categorical(multiple) field. Field [40011](http://biobank.ctsu.ox.ac.uk/showcase/field.cgi?id=40011) stores the histology of cancer tumours but is stored in 31 instances (rather than arrays) so we specify YES-INSTANCES so this field is treated as a categorical (multiple) field and the instances are treated as arrays.
5. DATA_CODING - The data coding IDs used to map a field to its data code in the data code information file described above. This is required for categorical (single) fields.

### Output
In the directory specified with the `resDir` argument, the following files will be created:

1. A .tsv file: A large file consisting of #individuals rows, and #PHESANT-phenotypes columns, which can be used for downstream analysis.
2. A .log file: One line for each Biobank field, providing information about the processing flow for this field.
3. Flow counts file: variable-flow-counts-all.txt - A set of counts denoting the number of variables reaching each point in the processing flow.

Where the phenome scan is run in parallel setup, then each parallel part will have one of each of the above files, with 'all' in each filename replaced with: [partIdx]-[numParts].

In addition, we also provide scripts to summarise the resultant phenotypes, taking the .tsv and .log files as input. These scripts use "variable-info/outcome-info.tsv" and "variable-info/data-coding-ordinal-info.csv" to create histograms and provide summary tables of the phenotype information.
