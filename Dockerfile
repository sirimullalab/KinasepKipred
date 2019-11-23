FROM hassanmohsin/rdkit-ubuntu:py2
MAINTAINER Md Mahmudulla Hassan <mhassan@miners.utep.edu>

WORKDIR app/
COPY get_kinase_pki.py utility.py requirements.txt Random_forest_gridsearch_py27.mdl ./
COPY mayachemtools ./mayachemtools

RUN pip install -r requirements.txt

ENTRYPOINT ["python", "get_kinase_pki.py"]
