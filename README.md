# BJ-CNV
BJ-CNV is a modified version of Ginkgo for copy number call. The pipeline takes deduplicated bam file, calculates coverage, performs normalization, and make cnv calls. Pipeline generates various plots and creates a table with copy number information with genomic coordinates. The CNV Tertiary pipeline is dependent on the secondary BJ-DNA-QC pipeline. Once the BJ-DNA-QC pipeline completes successfully, the prerequisites are met to launch the BJ-CNV pipeline.

# Pipeline Overview
Following are the steps and tools that pipeline users to perform the analyses:

- Convert BAM input files into BED files using bamtobed command from bedtools
- Map the BED coverage files into prior computed bins using Ginkgo module
- Perform segmentation, ploidy estimation, and create CNV plots using the Ginkgo module
- Perform final CNV call using the Ginkgo module

# Running Locally

Following are instructions for running BJ-CNV in a local Ubuntu server

## Install Java

```
sudo apt-get install default-jdk

java -version
```

## Install AWS CLI

```
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

## Install Nextflow

```
wget -qO- https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

```

## Install Docker

```
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

```
## Resources Required

For running the pipeline, a typical dataset requires 4 CPU cores and 14 GB of memory. For larger datasets, you may need to increase the resources to 8 CPU cores. You can specify these resources in the command as follows:
```
--max_cpus 4 --max_memory 14.GB
```

## Test Pipeline Execution

All pipeline resources are publically available at `s3://bioskryb-public-data/pipeline_resources` users need not have to download this, and will be downloaded during nextflow run.

**Command**

example-

** csv input **

```
git clone https://github.com/BioSkryb/bj-cnv.git
cd bj-cnv
nextflow run main.nf --input_csv $PWD/tests/data/inputs/input.csv --max_cpus 4 --max_memory 14.GB
```

**Input Options**

The input for the pipeline can be passed via a input.csv with a meta data.

- **CSV Metadata Input**: The CSV file should have 2 columns: `biosampleName`and `bam`. 
The `biosampleName` column contains the name of the biosample, `bam` has the path to the input bam files. Ensure that the index files are located in the same directory as the BAM files. For example:

```
biosampleName,bam
ResolveOME-test-bam1,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/bam/ResolveOME-test-bam1.bam
ResolveOME-test-bam2,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/bam/ResolveOME-test-bam2.bam
ResolveOME-test-bam3,s3://bioskryb-public-data/pipeline_resources/dev-resources/local_test_files/bam/ResolveOME-test-bam3.bam
```

**Outputs**

The pipeline saves its output files in the designated "publish_dir" directory. For details: [BJ-CNV outputs](https://docs.basejumper.bioskryb.com/pipelines/tertiary/bj-cnv/1.1.4/docs/#output-files)

**command options**

```
    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config


        [required]
        --input_csv         FILE    Path to input csv file

        --genome            STR     Reference genome to use. Available options - GRCh38, GRCm39
                                    DEFAULT: GRCh38


        [optional]
        
        --genomes_base      STR     Path to the genomes
                                    DEFAULT: s3://bioskryb-shared-data

        --publish_dir       DIR     Path to run output directory
                                    DEFAULT: 
                                    
        --timestamp         STR     User can specify timestamp otherwise uses runtime generated timestamp
                                    
        --help              BOOL    Display help message
```
**Tool versions**

- `bedtools: 2.28.0`
- `Ginkgo: 0.0.2`

**nf-test**

The BioSkryb BJ-CNV nextflow pipeline run is tested using the nf-test framework.

Installation:

nf-test has the same requirements as Nextflow and can be used on POSIX compatible systems like Linux or OS X. You can install nf-test using the following command:
```
wget -qO- https://code.askimed.com/install/nf-test | bash
```
It will create the nf-test executable file in the current directory. Optionally, move the nf-test file to a directory accessible by your $PATH variable.

Usage:

```
nf-test test
```

The nf-test for this repository is saved at tests/ folder.

```
    test("BJ-CNV Test") {

        when {
            params {
                // define parameters here. Example: 
                publish_dir = "${outputDir}/results"
                timestamp = "test"
            }
        }

        then {
            assertAll(
                // Check if the workflow was successful
                { assert workflow.success },

                // Verify existence of the AllSample-GinkgoSegmentSummary text file
                {assert new File("${outputDir}/results_test/tertiary_analyses/cnv_ginkgo/AllSample-GinkgoSegmentSummary.txt").exists()},

                // Verify existence of the SegCopy.binsize_1000000 text file
                {assert new File("${outputDir}/results_test/tertiary_analyses/cnv_ginkgo/SegCopy.binsize_1000000.tsv").exists()},

                // Verify existence of the cnv bin size jpeg file
                {assert new File("${outputDir}/results_test/tertiary_analyses/cnv_ginkgo/cnv_binsize_1000000_AAC7WVGM5-NA-ResolveOME-WM-Frag-Exome-DNA-SC10_CN.jpeg").exists()},

                // Check for a match in the cnv1 text file
                {assert snapshot (path("${outputDir}/results_test/tertiary_analyses/cnv_ginkgo/AAC7WVGM5-NA-ResolveOME-WM-Frag-Exome-DNA-SC10_CNV1.tsv")).match("cnv1")}

            )
        }

    }
```
