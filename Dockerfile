FROM python:3.11
RUN pip3 install pandas
RUN pip3 install scipy
#RUN pip3 install leidenalg
RUN pip3 install 'scanpy[leiden]'
RUN pip3 install louvain
# remove this when we update all comp environments
RUN pip3 install anndata==0.7.8 

COPY . /app