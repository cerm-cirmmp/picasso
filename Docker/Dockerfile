FROM docker.io/library/python:3.9.7

RUN pip install ortools==7.8.7959 pandas numpy scipy plotly matplotlib
RUN git clone https://github.com/cerm-cirmmp/picasso.git
RUN mkdir work
WORKDIR /work
