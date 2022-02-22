# Setup for coding

## On Ucloud

Access to [Ucloud](cloud.sdu.dk) with your account and access to your own project, or a project you have been invited to.
Look for the `Jupyter` app, and choose the most recent version to run. 

Open your terminal and download the repository for this course:

`git clone https://github.com/hds-sandbox/scRNASeq_course.git`

Download the datasets into the data folder (a few GigaBytes each)

```
mkdir -p scRNAseq_course/Data/notebooks_data

cd scRNAseq_course/Data/notebooks_data

wget https://aarhusuniversitet-my.sharepoint.com/:u:/g/personal/au612681_uni_au_dk/Ef-p9YNvrrZJp0s5gqRJ02YBg19brRA0o6JQ49S3MR-lYA?email=jose.romero%40sund.ku.dk&e=o5MuQw

wget https://aarhusuniversitet-my.sharepoint.com/:u:/g/personal/au612681_uni_au_dk/Ecvn2TQ_BfVIlAFB4QLvJrEBaBGx1Rqcr_itY0U8jlS2gg?email=jose.romero%40sund.ku.dk&e=RhJTDL

wget https://aarhusuniversitet-my.sharepoint.com/:u:/g/personal/au612681_uni_au_dk/EcfQhugmdARDrp-3QVOXFXABEc1EdrGBZ3dkeuit9vvrsw?email=jose.romero%40sund.ku.dk&e=tSYM5U
```

Now, install the required packages by creating a conda environment:

```
cd ../../Environments/Python
conda env create --file scrna_twodays.yml
```

Activate you conda environment by

```
conda activate scrna_test_environment
```

Finally, install a couple of remaining packages

```
Rscript R_install.R
```