FROM python:3.11-slim

RUN apt update && apt install -y gfortran 

RUN pip3 install pandas numpy

COPY . /app 
ENV PYTHONPATH="/app:$PYTHONPATH"
RUN mkdir -p /root/.config/autochem
COPY molecules.txt /root/.config/autochem/molecules.txt

ENTRYPOINT ["/app/bin/autochem"] 
