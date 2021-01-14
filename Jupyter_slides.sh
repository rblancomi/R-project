#Create a new environment (not mandatory but recommended)

conda create -y -n R_env

#activate the environment
conda activate R_env 

#Install R within Jupyter Notebook
conda install -c r r-essentials

#Install OncoSimulR
#v2.20.0 the same one that is available from RStudio
conda install -c bioconda bioconductor-oncosimulr

#Launch Jupyter notebook and create a new notebook with R
jupyter notebook

#Once you have the slides
#Use this to convert the notebook into slides that will open into the browser

jupyter nbconvert Slides.ipynb --to slides --post serve --SlidesExporter.reveal_theme=serif --SlidesExporter.reveal_scroll=True --SlidesExporter.reveal_transition=none
