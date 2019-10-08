FROM hassanmohsin/rdkit-ubuntu:py2
MAINTAINER Md Mahmudulla Hassan <mhassan@miners.utep.edu>

WORKDIR app/
COPY app/ /app/ #app contains this repo

RUN pip install -r requirements.txt

ENTRYPOINT ["python", "get_kinase_pki.py"]
