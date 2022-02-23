#!/usr/bin
eval "$(conda shell.bash hook)"
conda activate /work/sandbox_scRNA_testAndFeedback/scrna-environment
python -m ipykernel install --user --name scrna-environment --display-name "Python (scRNA)"
mkdir introduction_to_scrna_analysis
cd introduction_to_scrna_analysis
git init
git remote add origin https://github.com/hds-sandbox/scRNASeq_course.git
echo "Assignments/" >> .git/info/sparse-checkout
echo "Environments/" >> .git/info/sparse-checkout
echo "Notebooks/" >> .git/info/sparse-checkout
echo "Scripts/" >> .git/info/sparse-checkout
git pull --depth=1 origin main
