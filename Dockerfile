FROM python:3
RUN pip3 install pandas
RUN pip3 install scipy
#RUN pip3 install leidenalg
RUN pip3 install 'scanpy[leiden]'

COPY . /app